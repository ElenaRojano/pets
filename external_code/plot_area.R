#! /usr/bin/env Rscript
# x,y graph

library(ggplot2)
library(optparse)

################################################################
# OPTPARSE
################################################################
option_list <- list(
	make_option(c("-d", "--data_file"), type="character",
		help="Tabulated file with information about each sample"),
	make_option(c("-o", "--output"), type="character", default="results",
		help="Output figure file"),
	make_option(c("-x", "--x_values"), type="character", 
		help="Name of column with values to be plotted"),
	make_option(c("-y", "--y_values"), type="character", 
		help="Name of column with values to be plotted"),
	make_option(c("-f", "--density_values"), type="character", 
		help="Name of column to be used as density values"),
	make_option(c("-H", "--header"), action="store_false", default=TRUE,
        help="The input table not have header line"),
	 make_option(c("-X", "--x_title"), type="character", 
	 	help="Name of column to be used for bars titles"), 
	make_option(c("-Y", "--y_title"), type="character", 
	 	help="Title of y axis"),
	make_option(c("-F", "--output_format"), type="character", default="pdf", 
	 	help="pdf or jpeg file output format"),
    make_option(c("-m", "--maxs_file"), type="character", default="",
        help="Tabulated file maximum of each sample"),
    make_option(c("-t", "--graph_title"), type="character", default="",
        help="Title of the graph")	

)
opt <- parse_args(OptionParser(option_list=option_list))


################################################################
## MAIN
################################################################

data <- read.table(opt$data_file, sep="\t", header=opt$header)
if (opt$output_format == "pdf"){
	pdf(paste(opt$output, '.pdf', sep=""))
}else if(opt$output_format == "jpeg"){
	jpeg(paste(opt$output, '.jpeg', sep=""))
}	
	goodChrOrder <- c(1:22,"X","Y")
	data$V1 <- factor(data$V1,levels=goodChrOrder)

	maxs <- c()
	if(opt$maxs_file != ""){
		maxs <- read.table(opt$maxs_file, sep="\t", header=FALSE)
		#print(maxs)
	}
	#ggplot(data=data, aes(x=data[[opt$x_values]], y=data[[opt$y_values]] ))  +
	obj <- ggplot(data=data, aes(x=V2, y=V3 ))  
	#geom_area(aes(fill=data[[opt$density_values]], )) +
	obj <- obj + geom_area(aes(fill=V1, ))
	obj <- obj + facet_wrap(~ V1, ncol=2, strip.position = "right" )
	if(length(maxs) > 0){
		obj <- obj + geom_vline(data = maxs, aes(xintercept = V2)) 
	}
	obj <- obj + xlab(opt$x_title)
	obj <- obj + ylab(opt$y_title)
	obj <- obj + theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
	obj <- obj + guides(fill=FALSE)
	#obj <- obj + labs(title = opt$graph_title)
	obj <- obj + ggtitle(label = opt$graph_title)
	obj
dev.off()
