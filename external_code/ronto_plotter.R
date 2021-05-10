#!/usr/bin/env Rscript
#' @author SysBioLab

############################################################
##                        FUNCTIONS                       ##
############################################################

sparseLevelItems <- function(term_levels){
  u_levels <- unique(term_levels)
  base_radians <- rep(0,length(term_levels))
  for (level in u_levels) {
    terms_in_level <- which(term_levels == level)
    radial_pos <- seq(length(terms_in_level)) * 
                   ((2*pi)/length(terms_in_level))
    base_radians[terms_in_level] <- radial_pos 
  }
  return(base_radians)
}

prepare_data <- function(profile_terms, onto, white_list, expand_wl){
  profile_terms$Type <- "Profile"
  ontTerms2plot <- onto$id[!onto$obsolete]
  if(!is.null(white_list)){ # Populate with white list
    whitelist <- read.table(file = white_list, sep = "\t", stringsAsFactors = FALSE, header = FALSE)[,1]
    if(expand_wl){
      whitelist <- unique(c(whitelist, profile_terms[[1]]))
    }
    ontTerms2plot <- ontologyIndex::get_ancestors(onto, whitelist)
  }
  ref_terms <- setdiff(ontTerms2plot, unique(profile_terms[[1]]))
  ref_terms <- data.frame(Term= ref_terms)
  ref_terms[,1 + seq(ncol(profile_terms) - 2)] <- 0
  ref_terms$Type <- "Reference"
  colnames(ref_terms) <- colnames(profile_terms)
  terms_info <- rbind(profile_terms, ref_terms)
  return(terms_info)
}

calc_all_paths <- function(ids, ontology){
    env <- environment()
    paths <- list()
    calc_path(ids = ids, ontology = ontology, env = env)
    return(paths)
}

calc_path <- function(ids, ontology, env){
  for (id in ids) {
    if (is.null(env$paths[[id]])) {
      env$paths[[id]] <- list()
      direct_ancestors <- ontology$parents[[id]]
      if (length(direct_ancestors) == 0) {
        env$paths[[id]] <- c(env$paths[[id]], list(c(id)))
      } else {
        for (anc in direct_ancestors) {
          calc_path(ids = anc, ontology = ontology,env = env)
          env$paths[[id]] <- c(env$paths[[id]], 
                               lapply(env$paths[[anc]], function(path){
                                c(id,unlist(path))}))
        }
      }
    }
  }
}


############################################################
##                    OPTIONS PARSING                     ##
############################################################

option_list <- list(
  optparse::make_option(c("-i", "--input_file"), type="character", default=NULL,
    help="Input file to be loaded. Format must be: 1) Term, 2) Value [color], 3) (OPTIONAL) Second value [size]"),
  optparse::make_option(c("-o", "--output_file"), type="character", default=NULL,
    help="Output graph file name"),
  optparse::make_option(c("-O", "--ontology"), type="character", default=NULL,
    help="Ontology OBO file to be loaded"),
  optparse::make_option(c("-w", "--white_list"), type="character", default=NULL,
    help="White list file of terms to be plotted (plus input file terms). If not specified, all ontology terms will be plotted"),
  optparse::make_option(c("-H", "--header"), type="logical", default=FALSE, action = "store_true",
    help="Flag to indicate that input file has header"),
  optparse::make_option(c("-e", "--expand_wl"), type="logical", default=FALSE, action = "store_true",
    help="Flag to indicate that white list must be expanded using input terms and their parentals"),
  optparse::make_option(c("-v", "--verbose"), type="logical", default=FALSE, action = "store_true",
    help="Activate verbose mode")
)
opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))

############################################################
##                           MAIN                         ##
############################################################

time_ref <- Sys.time()

# Load data
profile_terms <- read.table(file = opt$input_file, sep = "\t", stringsAsFactors = FALSE, header = opt$header) 
if (!opt$header){
  cnames <- c("Term","Color","Size")
  colnames(profile_terms) <- cnames[1:ncol(profile_terms)]
} 

# Load ontology
onto <- ontologyIndex::get_ontology(file = opt$ontology)

# Populate input df with all allowed terms

terms_info <- prepare_data(profile_terms, onto, white_list = opt$white_list, expand_wl = opt$expand_wl)

# Obtain levels
paths <- calc_all_paths(ids = terms_info[[1]], ontology = onto)
levels <- lapply(paths, function(id){min(lengths(id))})

# Sort by code 
terms_info <- terms_info[order(terms_info[[1]]),]

# Include plotting data
terms_info$Level <- unlist(lapply(terms_info[[1]],function(term){levels[term]}))
terms_info$Radians <- sparseLevelItems(terms_info$Level)
terms_info$LevelD <- factor(as.character(terms_info$Level), levels = as.character(c(0,seq(max(terms_info$Level)))))

aes_profiles <- ggplot2::aes_string(color = colnames(terms_info)[2])
if(length(colnames(terms_info)) == 4){
  aes_profiles <- c(aes_profiles, ggplot2::aes_string(size = colnames(terms_info)[3]))
}

pp <- ggplot2::ggplot(mapping = ggplot2::aes(x = LevelD, y = Radians)) + 
      ggplot2::geom_point(data = terms_info[terms_info$Type == "Reference",], alpha = 0.3, size = 0.8) + 
      ggplot2::geom_point(data = terms_info[terms_info$Type == "Profile",],mapping = aes_profiles, alpha = 0.7) + 
      ggplot2::scale_colour_gradient(low = "orange", high = "red") + 
      ggplot2::scale_x_discrete("LevelD") + 
      ggplot2::ylim(c(0,2*pi)) +
      ggplot2::coord_polar(theta = "y")


png(paste0(opt$output_file,".png"), height = 1000, width = 1000)
  plot(pp)
dev.off()


