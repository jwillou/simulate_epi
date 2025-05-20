#Set working directory, import packages, source functions
setwd(paste(directory,"/source/", sep = ''))    # set temp working directory 

#import packages
library(dplyr)
library(hierfstat)
library(tidyr)

#source functions
source(paste(getwd(), "/CalcFstEpi.R", sep = ''))
source(paste(getwd(), "/CalcFstGen.R", sep = ''))
source(paste(getwd(), "/Simulate.R", sep = ''))
source(paste(getwd(), "/Theory_simple.R", sep = ''))
source(paste(getwd(), "/Theory_migration.R", sep = ''))
source(paste(getwd(), "/Theory_full.R", sep = ''))
source(paste(getwd(), "/Theory_WP.R", sep = ''))
source(paste(getwd(), "/Theory_adjusted_migration.R", sep = ''))
source(paste(getwd(), "/Theory_accumulative_migration.R", sep = ''))