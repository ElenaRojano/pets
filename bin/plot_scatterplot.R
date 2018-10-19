#! /usr/bin/env Rscript

library(ggplot2)
args <- commandArgs(trailingOnly = TRUE)

data <- read.table(args[1], header=TRUE)
output <- args[2]
x_axis <- args[3]
y_axis <- args[4]
density <- args[5]
x_tag <- args[6]
y_tag <- args[7]
size_tag <- args[8]
x_order <- unique(data[[x_axis]])
data[[x_axis]] <- factor(data[[x_axis]], levels = x_order)

pdf(file.path(output, 'scatterplot.pdf'))
	ggplot(data, aes(x=data[[x_axis]], y=data[[y_axis]])) + 
		geom_point(aes(size=data[[density]])) +
		xlab(x_tag) +
		ylab(y_tag) +
		labs(size = size_tag) +
		theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
