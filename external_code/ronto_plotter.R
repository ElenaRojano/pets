#!/usr/bin/env Rscript

#' @author Fernando Moreno Jabato <jabato(at)uma(dot)es>

############################################################
##                      SETUP PROGRAM                     ##
############################################################
# Loading libraries
suppressPackageStartupMessages(require(optparse)) 
suppressPackageStartupMessages(require(ggplot2)) 
suppressPackageStartupMessages(require(ontologyIndex)) 
suppressPackageStartupMessages(require(pbapply)) 

# Prepare command line input 
option_list <- list(
  make_option(c("-i", "--input_file"), type="character", default=NULL,
    help="Input file to be loaded. Format must be: 1) Term, 2) Value [color], 3) (OPTIONAL) Second value [size]"),
  make_option(c("-o", "--output_file"), type="character", default=NULL,
    help="Output graph file name"),
  make_option(c("-O", "--ontology"), type="character", default=NULL,
    help="Ontology OBO file to be loaded"),
  make_option(c("-w", "--white_list"), type="character", default=NULL,
    help="White list file of terms to be plotted (plus input file terms). If not specified, all ontology terms will be plotted"),
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


findPaths <- function(ontology, verbose = FALSE){
  calculated_paths <- list()
  appfun <- ifelse(verbose,pblapply,lapply)
  invisible(appfun(ontology$id,function(id){
    if(!id %in% names(calculated_paths)){
      paths <- findPath(id,ontology)
      calculated_paths[id] <<- paths
      # Add already expanded terms
    }
  }))
  return(calculated_paths)
}


findShortestPath <- function(term,ontology,calculated_paths = NULL){
  # Check
  if(!is.null(calculated_paths)){
    if(term %in% names(calculated_paths)){
      return(min(unlist(lapply(calculated_paths[term],length)))-1)
    }
  }

  paths <- findPath(term,ontology)
  return(min(unlist(lapply(paths,length)))-1)
}

findPath <- function(term,ontology,calculated_paths = NULL){
  # Check
  if(!is.null(calculated_paths)){
    if(term %in% names(calculated_paths)){
      return(calculated_paths[term])
    }
  }
  # Take parentals
  parents <- ontology$parents[term]
  # Expand if possible
  if(length(parents) > 0){
    expansion <- lapply(parents,function(p){findPath(p,ontology)})
  }else{
    expansion <- list() 
  }
  # Concat expansions
  expansion_updated <- lapply(expansion,function(ex){append(term,ex)})
  # Return
  return(expansion_updated)
}


############################################################
##                           MAIN                         ##
############################################################

apply_fun <- lapply
if(opt$verbose){ # Verbose point
  message("Loading input file")
  apply_fun <- pblapply
  pboptions(type="txt")
} 

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
topopulate <- onto$id
if(!is.null(opt$white_list)){ # Populate with white list
  wl <- read.table(file = opt$white_list, sep = "\t", stringsAsFactors = FALSE, header = FALSE)[,1]
  wl <- unique(c(wl,df$Term))
  topopulate <- unique(unlist(lapply(wl,function(item){get_ancestors(onto,item)})))
}
invisible(apply_fun(setdiff(topopulate,unique(df$Term)),function(id){
  df <<- rbind(df, as.list(c(id, rep(0, ncol(df)-2), "Added")))
}))  

if(opt$verbose) message("Calculating levels ...") # Verbose point

# Obtain levels
# paths <- findPaths(onto,verbose = opt$verbose)
# levels <- apply_fun(onto$id,function(id){findShortestPath(id,onto,paths)})
levels <- apply_fun(onto$id,function(id){length(get_ancestors(onto,id)) - 1})
names(levels) <- onto$id

# Sort by code to be (always) the same
df <- df[order(df$Term),]

# Include plotting data
df$Level <- unlist(lapply(df$Term,function(term){levels[term]}))
df$Radians <- sparseLevelItems(df$Level)
df$LevelD <- factor(as.character(df$Level), levels = as.character(c(0,seq(max(df$Level)))))
df$Color <- as.numeric(df$Color)

if(opt$verbose) message("Generating plot") # Verbose point

if("Size" %in% colnames(df)){
  aes_aux <- aes(x = LevelD, y = Radians, color = Color, size = Size)
}else{
  aes_aux <- aes(x = LevelD, y = Radians, color = Color)
}

pp = ggplot() + 
      geom_point(data = df[df$Type == "Added",],mapping = aes(x = LevelD, y = Radians), alpha = 0.3, size = 0.7) + 
      geom_point(data = df[df$Type == "Raw",],mapping = aes_aux, alpha = 0.6) + 
      scale_colour_gradient(low = "orange", high = "red") + 
      scale_x_discrete("LevelD") + 
      ylim(c(0,2*pi)) +
      coord_polar(theta = "y")

if(opt$verbose) message("Rendering plot ...") # Verbose point

siz <- 1000
png(paste0(opt$output_file,".png"), height = siz, width = siz)
plot(pp)
dev.off()

if(opt$verbose) message("Program finish") # Verbose point

