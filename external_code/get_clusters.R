#! /usr/bin/env Rscript
library(RcppCNPy)
library(optparse)
library(fastcluster)


option_list <- list(
	make_option(c("-d", "--data_file"), type="character",
		help="Tabulated file with information about each sample"),
	make_option(c("-y", "--npy"), type="character", default=NULL,
		help="Indicates that input file is a numpy matrix and the given PATH is a file with axis labels"),
	make_option(c("-o", "--output"), type="character", default="output",
		help="Output figure file")

)
opt <- parse_args(OptionParser(option_list=option_list))

if(!is.null(opt$npy)){
	x_axis_labels <- read.table(paste0(opt$npy, '_x.lst'), header=FALSE, stringsAsFactors=FALSE)
	y_axis_labels <- read.table(paste0(opt$npy, '_y.lst'), header=FALSE, stringsAsFactors=FALSE)
	matrix_data <- npyLoad(opt$data_file)
	colnames(matrix_data) <- x_axis_labels$V1
	rownames(matrix_data) <- y_axis_labels$V1
}else{
	matrix_data <- read.table(opt$data_file, sep="\t", quote = '')
}

d <- dist(matrix_data, method = "euclidean") # distance matrix
fit <- fastcluster::hclust(d, method="ward.D2")
fit$height <- round(fit$height, 6) 
groups <- cutree(fit, h=1.5) 
write.table(groups, file=opt$output, sep="\t", quote=FALSE, col.names=FALSE)

pdf(file.path(dirname(opt$output), 'figures.pdf'))
	plot(fit) # display dendogram
	rect.hclust(fit, h=1.5, border="red") # draw dendogram with red borders around the 5 clusters
dev.off()
