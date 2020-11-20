#! /usr/bin/env Rscript

library(RcppCNPy)
library(optparse)
library(gplots)
library("RColorBrewer")

#####################
## FUNCTIONS
#####################
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
		help="Matrix transformation parameter. Options: comp1, inverse, max.")

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
		data <- as.matrix(read.table(opt$data_file, sep="\t", header=opt$header, stringsAsFactors=FALSE, row.names="pat"))
	}
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
	quantValue <- quantile(matrix_transf, c(.2), na.rm = TRUE)
	hr <- hclust(as.dist(matrix_transf), method="ward.D2")
	groups <- cutree(hr, h = quantValue)
	######### EXPORT
	write.table(groups, file=paste(opt$output, '_clusters.txt', sep=''), sep="\t", quote=FALSE, col.names=FALSE)
}else{
	# Calc similitudes of rows
	mdistRows = toDistances(matrix_transf)
	mdistCols = toDistances(matrix_transf, FALSE)
	# Obtaing clustering
	quantValue_row <- quantile(mdistRows, c(.2), na.rm = TRUE)
	hr_row <- hclust(as.dist(mdistRows), method="ward.D2")
	groups_row <- cutree(hr_row, h = quantValue_row)
	quantValue_col <- quantile(mdistCols, c(.2), na.rm = TRUE)
	hr_col <- hclust(as.dist(mdistCols), method="ward.D2")
	groups_col <- cutree(hr_col, h = quantValue_col)
	######### EXPORT
	write.table(groups_row, file=paste(opt$output, '_clusters_rows.txt', sep=''), sep="\t", quote=FALSE, col.names=FALSE)
	write.table(groups_col, file=paste(opt$output, '_clusters_cols.txt', sep=''), sep="\t", quote=FALSE, col.names=FALSE)
}

######### RENDER
pdf(paste(opt$output, '_heatmap.pdf', sep=''))
	if(opt$same_sets){
		heatmap.2(data, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hr), trace="none", col=brewer.pal(11,"RdBu"), dendrogram = c("row"), labRow = FALSE, labCol = FALSE,
					xlab = opt$collabel, ylab = opt$rowlabel)
	}else{
		heatmap.2(data, Rowv=as.dendrogram(hr_row), Colv=as.dendrogram(hr_col), trace="none", col=brewer.pal(11,"RdBu"), labRow = FALSE, labCol = FALSE,
					xlab = opt$collabel, ylab = opt$rowlabel)

	}
dev.off()

