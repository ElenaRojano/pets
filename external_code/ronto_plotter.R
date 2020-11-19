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
  make_option(c("-e", "--expand_wl"), type="logical", default=FALSE, action = "store_true",
    help="Flag to indicate that white list must be expanded using input terms and their parentals"),
  make_option(c("-v", "--verbose"), type="logical", default=FALSE, action = "store_true",
    help="Activate verbose mode")
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


calc_paths <- function(ontology, ids = NULL){
  # Initialize
  visited_terms <- c()
  paths <- list()
  if(is.null(ids)){
    ids <- as.vector(ontology$id)
  }else{
    ids <- as.vector(ids)
  }
  env <- environment()
  for(id in ids){
    if(!id %in% visited_terms){ # Only if it's not visited
      paths[[id]] <- list(ID = id, NPaths = 0, ShortestPath = 0, Paths = list())
      direct_ancestors <- ontology$parents[[id]]
      if(length(direct_ancestors) <= 0){
        paths[[id]]$Paths <- append(paths[[id]]$Paths, list(c(id)))
      }else{
        invisible(lapply(direct_ancestors,function(anc){
          if(!anc %in% visited_terms){
            calc_path(anc,ontology,env)
          }
          paths[[id]]$Paths <<- append(paths[[id]]$Paths, lapply(paths[[anc]]$Paths,function(path){c(id,unlist(path))}))
          visited_terms <<- c(visited_terms,anc)
        }))
      }
      visited_terms <- c(visited_terms,id)
    }
    # Calc metadata
    paths[[id]]$NPaths <- length(paths[[id]]$Paths)
    paths[[id]]$ShortestPath <- min(unlist(lapply(paths[[id]]$Paths,function(path){length(unlist(path))})))
  }
  return(paths)
}


calc_path <- function(id, ontology, env){
  if(!id %in% env$visited_terms){
    env$paths[[id]] <- list(ID = id, NPaths = 0, ShortestPath = 0, Paths = list())
    direct_ancestors <- ontology$parents[[id]]
    if(length(direct_ancestors) <= 0){
      env$paths[[id]]$Paths <- append(env$paths[[id]]$Paths, list(c(id)))
    }else{
      invisible(lapply(direct_ancestors,function(anc){
        if(!anc %in% env$visited_terms){
          calc_path(anc,ontology,env)
        }
        env$paths[[id]]$Paths <- append(env$paths[[id]]$Paths, lapply(env$paths[[anc]]$Paths,function(path){c(id,unlist(path))}))
        env$visited_terms <- c(env$visited_terms,anc)
      }))
    }
    env$visited_terms <- c(env$visited_terms,id)
  }
}

############################################################
##                           MAIN                         ##
############################################################

apply_fun <- lapply
if(opt$verbose){ # Verbose point
  message("Loading input file")
  apply_fun <- pblapply
  pboptions(type="txt")
  ltimes <- list()
  ltimes$init <- Sys.time()
} 

# Load data
df <- read.table(file = opt$input_file, sep = "\t", stringsAsFactors = FALSE, header = opt$header) 
cnames <- c("Term","Color","Size")
colnames(df) <- cnames[1:ncol(df)]

if(opt$verbose){
  ltimes$load_init <- Sys.time()
  message("Loading ontology ...") # Verbose point
}

# Load ontology
onto <- get_ontology(file = opt$ontology)

if(opt$verbose){
  ltimes$populate_init <- Sys.time()
  message(paste0("\tElapsed : ",round(as.numeric(ltimes$populate_init) - as.numeric(ltimes$load_init),2)," s"))
  message("Populating with all available terms ...") # Verbose point
}

# Populate input df with all allowed terms
df$Type <- rep("Raw",nrow(df))
topopulate <- onto$id
if(!is.null(opt$white_list)){ # Populate with white list
  wl <- read.table(file = opt$white_list, sep = "\t", stringsAsFactors = FALSE, header = FALSE)[,1]
  if(opt$expand_wl) wl <- unique(c(wl,df$Term))
  topopulate <- unique(unlist(lapply(wl,function(item){get_ancestors(onto,item)})))
}

aux_df <- as.data.frame(do.call(rbind,apply_fun(setdiff(topopulate,unique(df$Term)),function(id){
  return(t(as.data.frame(c(id, rep(0, ncol(df)-2), "Added"))))
})))
rownames(aux_df) <- NULL
colnames(aux_df) <- colnames(df)
df <- rbind(df, aux_df)

if(opt$verbose){
  ltimes$levels_init <- Sys.time()
  message(paste0("\tElapsed : ",round(as.numeric(ltimes$levels_init) - as.numeric(ltimes$populate_init),2)," s"))
  message("Calculating levels ...") # Verbose point
}

# Obtain levels
uniqterms <- unique(df$Term)
paths <- calc_paths(onto, uniqterms)
levels <- lapply(paths,function(path){path$ShortestPath})

# Sort by code to be (always) the same
df <- df[order(df$Term),]

# Include plotting data
df$Level <- unlist(lapply(df$Term,function(term){levels[term]}))
df$Radians <- sparseLevelItems(df$Level)
df$LevelD <- factor(as.character(df$Level), levels = as.character(c(0,seq(max(df$Level)))))
df$Color <- as.numeric(df$Color)

if(opt$verbose){
  ltimes$plot_init <- Sys.time()
  message(paste0("\tElapsed : ",round(as.numeric(ltimes$plot_init) - as.numeric(ltimes$levels_init),2)," s"))
  message("Generating plot") # Verbose point
}

if("Size" %in% colnames(df)){
  aes_aux <- aes(x = LevelD, y = Radians, color = Color, size = Size)
}else{
  aes_aux <- aes(x = LevelD, y = Radians, color = Color)
}

pp = ggplot() + 
      geom_point(data = df[df$Type == "Added",],mapping = aes(x = LevelD, y = Radians), alpha = 0.3, size = 0.8) + 
      geom_point(data = df[df$Type == "Raw",],mapping = aes_aux, alpha = 0.7) + 
      scale_colour_gradient(low = "orange", high = "red") + 
      scale_x_discrete("LevelD") + 
      ylim(c(0,2*pi)) +
      coord_polar(theta = "y")

if(opt$verbose){
  ltimes$render_plot <- Sys.time()
  message(paste0("\tElapsed : ",round(as.numeric(ltimes$render_plot) - as.numeric(ltimes$plot_init),2)," s"))
  message("Rendering plot ...") # Verbose point
}

siz <- 1000
png(paste0(opt$output_file,".png"), height = siz, width = siz)
plot(pp)
dev.off()

if(opt$verbose){
  ltimes$end <- Sys.time()
  message(paste0("\tElapsed : ",round(as.numeric(ltimes$end) - as.numeric(ltimes$render_plot),2)," s"))
  message(paste0("Program finished. Elapsed : ",round(as.numeric(ltimes$end) - as.numeric(ltimes$init),2)," s"))
}

