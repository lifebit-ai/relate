#!/usr/bin/env Rscript

# @Author: daniel 
# @Date: 2019-11-11T10:44:09+00:00
# @Email:daniel.rhodes@genomicsengland.co.uk
# @Project: Aggregate_QC 
# @Last modified by:daniel  

options(bitmap='cairo')

#Project setup
if(grepl('^hpc-|^phpgridz', Sys.info()['nodename'])){
  myrepo = getOption("repos")
  myrepo["CRAN"] = "https://artifactory.gel.zone/artifactory/cran"
  options(repos = myrepo)
  rm(myrepo)
  pack_load <- function(x) {
    if (!require(x, character.only = TRUE)) {
      cat('Install package in research environment')
    }
    else{   
      require(x, character.only = TRUE, 
                    lib.loc = '~/public_data_resources/Rpackages/3.6.1/' )
    }
  }
}


if(grepl('daniel', Sys.info()['effective_user'])){
  pack_load <- function(x, install_only = FALSE) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      if (!install_only) {
        require(x, character.only = TRUE)
      }
    }
  }
}  

lapply(c("data.table",
         "magrittr",
         "ggplot2",
         'dplyr',
         'reshape2',
         'tidyr',
         'stringr'),
       pack_load)