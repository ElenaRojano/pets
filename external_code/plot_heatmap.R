#! /usr/bin/env Rscript

library(optparse)
library(gplots)
library("RColorBrewer")

#####################
## OPTPARSE
#####################
option_list <- list(
	make_option(c("-d", "--data_file"), type="character",
		help="Tabulated file with information about each sample"),
	make_option(c("-o", "--output"), type="character", default="output",
		help="Output figure file"),
	make_option(c("-m", "--matrix_transformation"), type="character", default=NULL,
		help="Matrix transformation parameter. Options: comp1, inverse, max.")

)
opt <- parse_args(OptionParser(option_list=option_list))


################################################################
## MAIN
################################################################

data <- as.matrix(read.table(opt$data_file, sep="\t", header=TRUE, stringsAsFactors=FALSE, row.names="pat"))

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

quantValue <- quantile(matrix_transf, c(.2), na.rm = TRUE)

hr <- hclust(as.dist(matrix_transf), method="ward.D2")

pdf(paste(opt$output, '_heatmap.pdf', sep=''))
	heatmap.2(data, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hr), trace="none", col=brewer.pal(11,"RdBu"), dendrogram = c("row"), labRow = FALSE, labCol = FALSE)
dev.off()

groups <- cutree(hr, h = quantValue)
write.table(groups, file=paste(opt$output, '_clusters.txt', sep=''), sep="\t", quote=FALSE, col.names=FALSE)