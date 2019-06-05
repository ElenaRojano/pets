#! /usr/bin/env Rscript

library(ggplot2)
library(optparse)

#####################
## OPTPARSE
#####################
option_list <- list(
	make_option(c("-d", "--data_file"), type="character",
		help="Tabulated file with information about each sample"),
	make_option(c("-o", "--output"), type="character", default="results",
		help="Output figure file"),
	make_option(c("-x", "--x_column"), type="character", 
		help="Name of column to be used for X dimension"),
	make_option(c("-y", "--y_column"), type="character", 
		help="Name of column to be used for Y dimension"),
	make_option(c("-s", "--set_column"), type="character", default="",
        	help="Name of column to be used on set groups"),
	make_option(c("-L", "--no_legend"), action="store_true", default=FALSE,
	       	help="Remove legend"),
	make_option(c("-l", "--legend_title"), type="character", default="Association methods",
	       	help="Title for legend"),
	make_option(c("-c", "--colours"), type="character", default="",
		help="Define which color is asigned to each data series. List colours comma separated."),
	make_option(c("-m", "--set_geom"), type="character", default="line", 
		help="Choose the type of graphical representation, using points or lines"),
	make_option(c("-e", "--establish_limits"), action="store_true", default=FALSE,
 	 	help="Allow establishing limits for X and Y axis. If true, please set x_limit and y_limit"),
	make_option(c("-X", "--x_limit"), type="integer", default=0,
	 	help="Set x axis limit"),
	make_option(c("-Y", "--y_limit"), type="integer", default=1,
		help="Set y axis limit")

)
opt <- parse_args(OptionParser(option_list=option_list))


################################################################
## MAIN
################################################################

data <- read.table(opt$data_file, sep="\t", header=TRUE)

pdf(paste(opt$output, '.pdf', sep=""))
	if(opt$set_column != ""){
		obj <- ggplot(data, aes(x=data[[opt$x_column]], y=data[[opt$y_column]], color=data[[opt$set_column]]))
	}else{
		obj <- ggplot(data, aes(x=data[[opt$x_column]], y=data[[opt$y_column]]))
	}
        if(opt$colours != ""){
                colours <- unlist(strsplit(opt$colours, ','))
                obj <- obj + scale_color_manual(values=c(colours))
        }

	if(opt$set_geom == 'point'){ 	
 		obj <- obj + geom_point()
 	}else if(opt$set_geom == 'line'){
 		obj <- obj + geom_line()
	}

	obj <- obj + xlim(0, 1)
	obj <- obj + ylim(0, 1)
	
	if (opt$establish_limits){
	 	obj <- obj + xlim(opt$x_limit, 1)
		obj <- obj + ylim(0, opt$y_limit)
	 }
	obj <- obj + xlab(opt$x_column)
	obj <- obj + ylab(opt$y_column)
	obj <- obj + labs(fill = opt$legend_title)
	if(opt$no_legend){
		obj <- obj + guides(color=FALSE)	
	}
	obj
dev.off()

