#! /usr/bin/env Rscript

library(optparse)

option_list <- list(
  optparse::make_option(c("-i", "--info"), type="character", default=NULL,
                        help="Input file - needs minimum 2 columns"),
  optparse::make_option(c("-o", "--output"), type="character", default="out.pdf",
                        help="Output file name"),
  optparse::make_option(c("-x", "--x_axis"), type="character", default=NULL,
                        help="column name for y-axis values."),
  optparse::make_option(c("-y", "--y_axis"), type="character", default=NULL,
                        help="column to use for y-axis values"),
  optparse::make_option("--x_tag", type="character", default="x-axis",
                        help="x-axis label"),
  optparse::make_option("--y_tag", type="character", default="y-axis",
                        help="y-axis label"),
  optparse::make_option("--x_lim", type="character", default="NA,NA",
                        help="x-axis limits (2 values, comma separated, NAs mean axes follow data)"),
  optparse::make_option("--y_lim", type="character", default="NA,NA",
                        help="y-axis limits (2 values, comma separated, NAs mean axes follow data)")
)
opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))

library(ggplot2)

info <- read.table(opt$info, header=TRUE)
x_lim <- as.numeric(strsplit(opt$x_lim, ",")[[1]])
y_lim <- as.numeric(strsplit(opt$y_lim, ",")[[1]])
x_lim_min <- x_lim[[1]]
x_lim_max <- x_lim[[2]]
y_lim_min <- y_lim[[1]]
y_lim_max <- y_lim[[2]]

pdf(file.path(opt$output))
	p <- ggplot(info, aes(x=.data[[opt$x_axis]], y=.data[[opt$y_axis]])) + 
		geom_point(alpha = 1/10) +
		xlim(x_lim_min, x_lim_max) +
		ylim(y_lim_min, y_lim_max) +
		xlab(opt$x_tag) +
		ylab(opt$y_tag) +
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

