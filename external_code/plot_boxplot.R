#! /usr/bin/env Rscript

library(ggplot2)
args <- commandArgs(trailingOnly = TRUE)

data <- read.table(args[1], header=TRUE)
output <- args[2]
x_axis <- args[3]
y_axis <- args[4]
x_tag <- args[5]
y_tag <- args[6]
y_axis2 <- args[7]
y_tag2 <- args[8]

x_order <- unique(data[[x_axis]])
data[[x_axis]] <- factor(data[[x_axis]], levels = x_order)
pdf(file.path(output, 'boxplot.pdf'))
	gridExtra::grid.arrange ( 
		ggplot(data, aes(x=data[[x_axis]], y=data[[y_axis2]])) + 
			geom_boxplot() +
			xlab(x_tag) +
			ylab(y_tag2) +
			theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()),
		ggplot(data, aes(x=data[[x_axis]], y=data[[y_axis]])) + 
			geom_boxplot() +
			xlab(x_tag) +
			ylab(y_tag) +
			theme(axis.text.x = element_text(angle = 45, hjust = 1)),
		ncol = 1)

dev.off()
