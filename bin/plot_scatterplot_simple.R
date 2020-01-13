#! /usr/bin/env Rscript

library(ggplot2)
args <- commandArgs(trailingOnly = TRUE)

data <- read.table(args[1], header=TRUE)
output <- args[2]
x_axis <- args[3]
y_axis <- args[4]
x_tag <- args[5]
y_tag <- args[6]

pdf(file.path(output))
	p <- ggplot(data, aes(x=data[[x_axis]], y=data[[y_axis]])) + 
		geom_point(alpha = 1/10) +
		xlab(x_tag) +
		ylab(y_tag) +
		theme(axis.text.x = element_text(angle = 45, hjust = 1))
	ggExtra::ggMarginal(
	  p,
	  type = 'density',
	  margins = 'both',
	  size = 5,
	  colour = '#000000',
	  fill = '#A6A39E'
	)
dev.off()

