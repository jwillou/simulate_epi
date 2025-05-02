#Set working directory, import packages, source functions
setwd(paste(directory,"/source/", sep = ''))    # set temp working directory 

#import packages
library(dplyr)
library(hierfstat)

#source functions
source(paste(getwd(), "/CalcFst.R", sep = ''))
source(paste(getwd(), "/Simulate.R", sep = ''))
source(paste(getwd(), "/Theory_simple.R", sep = ''))
source(paste(getwd(), "/Theory_migration.R", sep = ''))
source(paste(getwd(), "/Theory_full.R", sep = ''))
source(paste(getwd(), "/Theory_WP.R", sep = ''))
