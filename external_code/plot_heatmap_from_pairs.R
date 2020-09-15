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
	make_option(c("-H", "--header"), type="logical", default=FALSE, action="store_true",
		help="Indicates that file has header"),
	make_option(c("-m", "--matrix_transformation"), type="character", default=NULL,
		help="Matrix transformation parameter. Options: comp1, inverse, max.")

)
opt <- parse_args(OptionParser(option_list=option_list))


################################################################
## MAIN
################################################################

data_raw <- read.table(opt$data_file, sep="\t", header=opt$header, stringsAsFactors=FALSE)
colnames(data_raw) <- c("SetA","SetB","Value")

data <- matrix(0, ncol = length(unique(data_raw$SetB)), nrow = length(unique(data_raw$SetA)))
colnames(data) <- as.character(unique(data_raw$SetB))
rownames(data) <- as.character(unique(data_raw$SetA))
invisible(lapply(seq(nrow(data_raw)),function(i){
	data[as.character(data_raw$SetA[i]),as.character(data_raw$SetB[i])] <<- data_raw$Value[i]
}))

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