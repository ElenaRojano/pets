#! /usr/bin/env Rscript

library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

file <- args[1]

values <- args[2]

#x_axis_limit <- strtoi(args[3])
x_axis_limit <- as.numeric(args[3])

groups <- args[4]
x_axis_limit_min <- as.numeric(args[5])
#categories <- args[3]

#xtitle <- args[4]

#ytitle <- args[5]

data <- read.table(file, header = TRUE , sep="\t")

pdf('out.pdf')
	#ggplot(data, aes(x=data[[values]], colour=data[[categories]], fill=data[[categories]])) + 
		#geom_histogram(binwidth=.5, position="dodge") +
		#geom_histogram(position="dodge") +
		#xlab(xtitle) + 
		#ylab('Count') +
                if(is.na(groups)){
                        obj <- ggplot(data, aes(x=data[[values]])) 
                        obj <- obj + geom_density() 
                }else{
                        obj <- ggplot(data, aes(x=data[[values]], fill=data[[groups]]))
                        obj <- obj + geom_density(alpha=.3) 
                }
		obj <- obj + theme(legend.title=element_blank())
		if(!is.na(x_axis_limit)){
			xmin <- 0
			if(!is.na(x_axis_limit_min)){
				xmin <-x_axis_limit_min
			}
			obj <- obj + xlim(xmin, x_axis_limit) 
		}
		obj
dev.off()
