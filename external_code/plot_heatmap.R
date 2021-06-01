#! /usr/bin/env Rscript

library(RcppCNPy)
library(optparse)
library(gplots)
library(fastcluster)
library("RColorBrewer")

#####################
## FUNCTIONS
#####################

cluster_obj_to_groups <- function(matrix_transf, clust_obj, method, minProportionCluster) {
	max_clusters <- ceiling(ncol(matrix_transf)/2)
	min_clusters <- ceiling(ncol(matrix_transf)/10)
	# print("max");print(max_clusters)
	# print("min");print(min_clusters)
	if(method=="quantile") {
		quantValue <- quantile(matrix_transf, c(.2), na.rm = TRUE)
		groups <- cutree(clust_obj, h = quantValue)
	} else if (method == "min_height_increase") {
		inert.gain <- rev(clust_obj$height)
		intra <- rev(cumsum(rev(inert.gain)))
		quot <- intra[min_clusters:(max_clusters)]/intra[(min_clusters - 1):(max_clusters - 1)]
		nb.clust <- which.min(quot) + min_clusters -1
		groups <- cutree(clust_obj, k = nb.clust)
	} else if (method == "silhouette") {		
		# res <- NbClust::NbClust(diss=as.dist(matrix_transf), distance = NULL, min.nc=min_clusters, max.nc=max_clusters,  method = "ward.D2", index = "silhouette")
		# groups <- res$Best.partition
	} else if (method == "dynamic") {
		minClusterSize <- 2
		data_minClusterSize <- ceiling(ncol(matrix_transf) * minProportionCluster)
		minClusterSize <- max(c(minClusterSize, data_minClusterSize))
		# print(minClusterSize)
		groups <- dynamicTreeCut::cutreeDynamic(dendro = clust_obj, 
												distM = matrix_transf,
														deepSplit = 2, pamRespectsDendro = TRUE,
														minClusterSize = minClusterSize)
		names(groups) <- colnames(matrix_transf)
	} else {
		stop("Method not found")
	}
	return(groups)
}

get_matrix_subset_mean <- function(combo, matrix_transf, groups) {
  mean(matrix_transf[
		names(groups)[groups %in% combo[1]],
		names(groups)[groups %in% combo[2]]
      ]
  )
}

calc_sim_between_groups <- function(matrix_transf, groups) {
	unique_groups <- unique(groups)
	n_groups <- length(unique_groups)
	subset_means <- combn(unique_groups, m=2, 
		FUN=get_matrix_subset_mean, matrix_transf=matrix_transf, groups=groups)
	group_sim <- matrix(nrow = n_groups, ncol = n_groups, dimnames=list(unique_groups,unique_groups))
	group_sim[lower.tri(group_sim)] <- subset_means
	# print(1-group_sim)
	group_sim
}

vectdist <- function(vectA, vectB){
	# VectA and B must have same length. Exception not handled
	return(sqrt(sum((vectA - vectB)^2)))
}

toDistances <- function(vectors_matrix, rows = TRUE){
	if(!rows){
		vectors_matrix = t(vectors_matrix)
	}
	# Calc similitudes of rows
	numItems = nrow(vectors_matrix)
	Mdist = matrix(Inf,nrow = numItems, ncol = numItems)
	invisible(lapply(seq(numItems), function(i){
		if(i != numItems){
			invisible(lapply(seq(i+1, numItems), function(j){
				v = vectdist(vectors_matrix[i,],vectors_matrix[j,])
				Mdist[i,j] <<- v
				Mdist[j,i] <<- v
			}))
		}
		Mdist[i,i] <<- 0
	}))
	return(Mdist)
}



#####################
## OPTPARSE
#####################
option_list <- list(
	make_option(c("-d", "--data_file"), type="character",
		help="Tabulated file with information about each sample"),
	make_option(c("-H", "--header"), type="logical", default=FALSE, action="store_true",
		help="Indicates that file has header"),
	make_option(c("-P", "--pdf"), type="logical", default=FALSE, action="store_true",
		help="Indicates that file has header"),	
	make_option(c("-p", "--pairs"), type="logical", default=FALSE, action="store_true",
		help="Indicates if input file is a pairs-file instead a matrix"),
	make_option(c("-y", "--npy"), type="character", default=NULL,
		help="Indicates that input file is a numpy matrix and the given PATH is a file with axis labels"),
	make_option(c("-o", "--output"), type="character", default="output",
		help="Output figure file"),
	make_option(c("-s", "--same_sets"), type="logical", default=TRUE, action = "store_false",
		help="Flag to indicate that set A (rows) and B (columns) contain different items"),
	make_option(c("-c", "--collabel"), type="character", default="",
		help="Columns (x) graph label"),
	make_option(c("-r", "--rowlabel"), type="character", default="",
		help="Rows (y) graph label"),
	make_option(c("-m", "--matrix_transformation"), type="character", default=NULL,
		help="Matrix transformation parameter. Options: comp1, inverse, max."),
	make_option(c("-M","--minProportionCluster"), type="double", default=0.05,
		help="Minimun percentage of element per cluster."),
#	make_option(c("-C", "--min_clusters"), type="integer", default=NULL,
#		help="Minimum number of clusters to obtain"),
#	make_option(c("-M", "--max_clusters"), type="integer", default=NULL,
#		help="Maximum number of clusters to obtain"),
	make_option(c("-t", "--tree_cut_method"), type="character", default="quantile",
		help="Method to use to determine number of clusters in which to divide data"),
	make_option(c("-D", "--diagonal_replace"), type="numeric", default=NULL,
		help="Number to replace diagonal values in the input matrx")
)
opt <- parse_args(OptionParser(option_list=option_list))


################################################################
## MAIN
################################################################

######### LOAD
if(opt$pairs){ # Load pairs
	data_raw <- read.table(opt$data_file, sep="\t", header=opt$header, stringsAsFactors=FALSE)
	colnames(data_raw) <- c("SetA","SetB","Value")

	if(opt$same_sets){
		aux = unique(c(as.character(data_raw$SetA),as.character(data_raw$SetB)))
		rowSet = aux
		colSet = aux
	}else{
		rowSet = as.character(unique(data_raw$SetA))
		colSet = as.character(unique(data_raw$SetB))
	}

	data <- matrix(0, ncol = length(colSet), nrow = length(rowSet))
	colnames(data) <- colSet
	rownames(data) <- rowSet
	invisible(lapply(seq(nrow(data_raw)),function(i){
		data[as.character(data_raw$SetA[i]),as.character(data_raw$SetB[i])] <<- data_raw$Value[i]
	}))	
}else{ # Load matrix
	if(!is.null(opt$npy)){
		axis_labels <- read.table(opt$npy, header=FALSE, stringsAsFactors=FALSE)
		data <- npyLoad(opt$data_file)
		colnames(data) <- axis_labels$V1
		rownames(data) <- axis_labels$V1
	}else{
		data <- as.matrix(read.table(opt$data_file, sep="\t", header=opt$header, stringsAsFactors=FALSE, row.names= 1, check.names = FALSE))
	}
}

if(!is.null(opt$diagonal_replace)){
	diag(data) <- opt$diagonal_replace
}

######### NORMALIZE
if(is.null(opt$matrix_transformation)){
	matrix_transf <- data
} else if (opt$matrix_transformation == 'comp1'){
	matrix_transf <- 1 - data
} else if (opt$matrix_transformation == 'inverse'){
	matrix_transf <- 1 / data
} else if (opt$matrix_transformation == 'max'){
	matrix_transf <- max(data, na.rm = TRUE) - data	
} else {
	stop(paste('Invalid', opt$matrix_transformation, 'option value.'))
}


######### CLUSTERING
if(opt$same_sets){
	hr <- fastcluster::hclust(as.dist(matrix_transf), method="ward.D2")
	groups <- cluster_obj_to_groups(matrix_transf, hr, opt$tree_cut_method, minProportionCluster = opt$minProportionCluster)

	sim_between_groups <- calc_sim_between_groups(data, groups)
	distance_between_groups <- 1 - sim_between_groups
	groups_clustered <- fastcluster::hclust(as.dist(distance_between_groups), method="ward.D2")
	dendrogram_groups <- as.dendrogram(groups_clustered)

	# Plot dendrogram to check performance
	# png(file=file.path(opt$output, 'dendrogram_groups.png', sep=''))
	#    plot(dendrogram_groups)
	# dev.off()
	######### EXPORT
	write.table(groups, file=paste0(opt$output, '_clusters.txt'), sep="\t", quote=FALSE, col.names=FALSE, row.names= TRUE)
	save(dendrogram_groups, file=paste0(opt$output, '_dendrogram_groups.RData', sep=''))

}else{
	# Calc similitudes of rows
	mdistRows = toDistances(matrix_transf)
	mdistCols = toDistances(matrix_transf, FALSE)
	# Obtaing clustering
	quantValue_row <- quantile(mdistRows, c(.2), na.rm = TRUE)
	hr_row <- fastcluster::hclust(as.dist(mdistRows), method="ward.D2")
	groups_row <- cluster_obj_to_groups(mdistRows, hr, opt$tree_cut_method)

	quantValue_col <- quantile(mdistCols, c(.2), na.rm = TRUE)
	hr_col <- fastcluster::hclust(as.dist(mdistCols), method="ward.D2")
	groups_col <- cutree(hr_col, h = quantValue_col)
	######### EXPORT
	write.table(groups_row, file=paste(opt$output, '_clusters_rows.txt', sep=''), sep="\t", quote=FALSE, col.names=FALSE)
	write.table(groups_col, file=paste(opt$output, '_clusters_cols.txt', sep=''), sep="\t", quote=FALSE, col.names=FALSE)
}

######### RENDER
if(opt$pdf){
	pdf(paste0(opt$output, '_heatmap.pdf'), width = 1000, height = 1000, units = "px", res=175, pointsize = 8)
}else{
	png(paste0(opt$output, '_heatmap.png'), width = 1000, height = 1000, units = "px", res=175, pointsize = 8)
}
	if(opt$same_sets){
		group_colours <- colorRampPalette(brewer.pal(8, "Set1"))(length(unique(groups)))
		# print("cluster number:")
		# print(length(unique(groups)))
		# print(length(unique(group_colours)))
		group_colours_arranged <- group_colours[groups]
		heatmap.2(data, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hr), trace="none", col=brewer.pal(11,"RdBu"), dendrogram = c("row"), labRow = FALSE, labCol = FALSE,
					xlab = opt$collabel, ylab = opt$rowlabel, RowSideColors=group_colours_arranged)
	}else{
		heatmap.2(data, Rowv=as.dendrogram(hr_row), Colv=as.dendrogram(hr_col), trace="none", col=brewer.pal(11,"RdBu"), labRow = FALSE, labCol = FALSE,
					xlab = opt$collabel, ylab = opt$rowlabel)

	}
dev.off()
# save.image("test.RData")

