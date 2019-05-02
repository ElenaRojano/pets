#! /usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

file <- args[1]
output <- args[2]

# matrix_data <- read.table(file, sep="\t", header=TRUE, quote = '')
matrix_data <- read.table(file, sep="\t", quote = '')
d <- dist(matrix_data, method = "euclidean") # distance matrix
fit <- hclust(d, method="ward.D2")
fit$height <- round(fit$height, 6) 
groups <- cutree(fit, h=1.5) 
write.table(groups, file=file.path(output, 'cluster_asignation'), sep="\t", quote=FALSE, col.names=FALSE)

pdf(file.path(output, 'figures.pdf'))
	plot(fit) # display dendogram
	rect.hclust(fit, h=1.5, border="red") # draw dendogram with red borders around the 5 clusters
dev.off()
