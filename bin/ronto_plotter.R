#!/usr/bin/env Rscript

#' @author Fernando Moreno Jabato <jabato(at)uma(dot)es>

############################################################
##                      SETUP PROGRAM                     ##
############################################################
# Loading libraries
suppressPackageStartupMessages(require(optparse)) 
suppressPackageStartupMessages(require(ggplot2)) 
suppressPackageStartupMessages(require(ontologyIndex)) 

# Prepare command line input 
option_list <- list(
  make_option(c("-i", "--input_file"), type="character", default=NULL,
    help="Input file to be loaded. Format must be: 1) Term, 2) Value [color], 3) (OPTIONAL) Second value [size]"),
  make_option(c("-o", "--output_file"), type="character", default=NULL,
    help="Output graph file name"),
  make_option(c("-O", "--ontology"), type="character", default=NULL,
    help="Ontology OBO file to be loaded"),
  make_option(c("-H", "--header"), type="logical", default=FALSE, action = "store_true",
    help="Flag to indicate that input file has header"),
  make_option(c("-v", "--verbose"), type="logical", default=FALSE, action = "store_true",
    help="Input file to be loaded")
)
opt <- parse_args(OptionParser(option_list=option_list))


############################################################
##                        FUNCTIONS                       ##
############################################################

sparseLevelItems <- function(levels){
  ulvls <- unique(levels)
  base_radians <- rep(0,length(levels))
  invisible(lapply(ulvls,function(lvl){
    items <- which(levels == lvl)
    values = seq(length(items)) * ((2*pi)/length(items))
    base_radians[items] <<- values 
  }))
  return(base_radians)
}


############################################################
##                           MAIN                         ##
############################################################

if(opt$verbose) message("Loading input file") # Verbose point

# Load data
df <- read.table(file = opt$input_file, sep = "\t", stringsAsFactors = FALSE, header = opt$header) 
cnames <- c("Term","Color","Size")
colnames(df) <- cnames[1:ncol(df)]

if(opt$verbose) message("Loading ontology") # Verbose point

# Load ontology
onto <- get_ontology(file = opt$ontology)

if(opt$verbose) message("Populating with all available terms") # Verbose point

# Populate input df with all allowed terms
df$Type <- rep("Raw",nrow(df))
invisible(lapply(setdiff(onto$id,unique(df$Term)),function(id){
  df <<- rbind(df, as.list(c(id, rep(0, ncol(df)-2), "Added")))
}))

if(opt$verbose) message("Calculating levels") # Verbose point

# Obtain levels
levels <- lapply(onto$id,function(id){length(get_ancestors(onto,id)) - 1})
names(levels) <- onto$id

# Include plotting data
df$Level <- unlist(lapply(df$Term,function(term){levels[term]}))
df$Radians <- sparseLevelItems(df$Level)
df$LevelD <- factor(as.character(df$Level), levels = as.character(c(0,seq(max(df$Level)))))
df$Color <- as.numeric(df$Color)

# Sort by code to be (always) the same
df <- df[order(df$Term),]

if(opt$verbose) message("Generating plot") # Verbose point

if("Size" %in% colnames(df)){
  aes_aux <- aes(x = LevelD, y = Radians, color = Color, size = Size)
}else{
  aes_aux <- aes(x = LevelD, y = Radians, color = Color)
}

pp = ggplot() + 
      geom_point(data = df[df$Type == "Added",],mapping = aes(x = LevelD, y = Radians), alpha = 0.3, size = 0.7) + 
      geom_point(data = df[df$Type == "Raw",],mapping = aes_aux) + 
      scale_x_discrete("LevelD") + 
      ylim(c(0,2*pi)) +
      coord_polar(theta = "y")

if(opt$verbose) message("Rendering plot ...") # Verbose point

siz <- 1000
png(paste0(opt$output_file,".png"), height = siz, width = siz)
plot(pp)
dev.off()

if(opt$verbose) message("Program finish") # Verbose point

