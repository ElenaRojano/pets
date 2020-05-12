#! /usr/bin/env Rscript

library(optparse)

#####################
## OPTPARSE
#####################
option_list <- list(
	make_option(c("-d", "--data_file"), type="character",
		help="Tabulated file with information about each sample"),
	make_option(c("-o", "--output"), type="character", default="results",
		help="Output figure file"),
	make_option(c("-e", "--external_score"), type="double", default=NULL,
		help="Use external score"),
    make_option(c("-s", "--set_column"), type="character", default="",
       	help="Name of column to be converted to Z-scores")

)
opt <- parse_args(OptionParser(option_list=option_list))


################################################################
## MAIN
################################################################

data <- read.table(opt$data_file, sep="\t", header=FALSE)
raw_data <- data[[opt$set_column]]
if(!is.null(opt$external_score)){
	raw_data <- c(opt$external_score, raw_data)
}
z_scores = scale(raw_data, center=TRUE, scale=TRUE)
 if(!is.null(opt$external_score)){
	external_score2z_score <- z_scores[1]
	cat("ExtZScore\t", external_score2z_score, "\n")
	z_scores <- z_scores[-1] #remove external score
}

data[[opt$set_column]] <- z_scores

write.table(data, file=opt$output, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
