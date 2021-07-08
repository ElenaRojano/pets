#! /usr/bin/env Rscript

######################################################################################################
####################### LIBRARY INSTALLING SCRIPT for PETS #################################
######################################################################################################
print("Installing libraries from CRAN")
packages_list <-c("optparse","RcppCNPy","ggplot2","fastcluster","dplyr","gplots","RColorBrewer","tidyr","data.table","gridExtra", "dynamicTreeCut", "ggExtra", "ontologyIndex", "magrittr")
installed <- library()$results[,1]
packages_list <- setdiff(packages_list, installed)
install.packages(packages_list, repos='https://cloud.r-project.org')

