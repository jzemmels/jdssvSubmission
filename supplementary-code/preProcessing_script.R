## This file is used to pre-process the cartridge case scans in the DFSC-scans
## repo. Note that we assume that the file is saved in the home directory of the
## DFSC-scans repo (https://github.com/heike/DFSC-scans). You may need to change
## the path to wherever you've saved the repo.

# install.packages("cmcR")
library(cmcR)
# install.packages("x3ptools")
library(x3ptools)
# install.packages("tidyverse")
library(tidyverse)

source("scanNames.R")

if(!dir.exists("sample_400_preprocessed")){
  dir.create("sample_400_preprocessed")
}

# Define pre-processing function to be applied to each scan
preprocessScans <- function(file){

  # group number defines the sub-directory the scan is in
  groupNumber <- paste0("G",str_sub(file,2,4))

  # could be the first time we're processing a scan from a specific group
  if(!dir.exists(pastse0("sample_400_preprocessed/",groupNumber))){
    dir.create(pastse0("sample_400_preprocessed/",groupNumber))
  }

  # read it from the sample_400 file, then perform the necessary pre-processing
  # below
  ret <- x3ptools::read_x3p(paste0("sample_400/",groupNumber,"/",file))

  if(!is.null(ret$mask)){

    maskValues <- unique(c(ret$mask))

    # some mask values have an extra "FF" at the end
    # that isn't part of the hexidecimal color ID
    for(color in maskValues){

      ret$mask[ret$mask == color] <-
        str_sub(color,1,7)

    }
    maskCounts <- table(c(ret$mask))

    # assume that the least common color inthe mask
    # is the one to keep
    maskColorKeep <- names(maskCounts)[which.min(maskCounts)]

    # replace non-mask values with NA
    ret$surface.matrix[t(as.matrix(ret$mask)) != maskColorKeep] <- NA

    ret <- ret %>%
      cmcR:::preProcess_cropWS(robust = FALSE,
                               croppingProp = .5)

  }

  # apply automatic pre-processing
  ret <- ret %>%
    cmcR::preProcess_removeTrend(statistic = "quantile",
                                 tau = .5,
                                 method = "fn") %>%
    cmcR::preProcess_gaussFilter() %>%
    cmcR::preProcess_removeTrend(statistic = "quantile",
                                 tau = .5,
                                 method = "fn") %>%
    cmcR::preProcess_gaussFilter() %>%
    cmcR::preProcess_erode(region = "interior",
                           morphRadius = 12) %>%
    cmcR::preProcess_erode(region = "exterior",
                           morphRadius = 12)

  # cropping whitespace will make the dimensions of the surface matrix fit to
  # the non-missing observations in the scan
  ret <- ret %>%
    cmcR:::preProcess_cropWS(robust = FALSE,
                             croppingProp = .5)

  ret$mask <- NULL

  x3ptools::x3p_write(ret,file = paste0("sample_400_preprocessed/",groupNumber,"/",file))

}

# Use future package if you want to parallelize the processing
# install.packages("future")

# future:::ClusterRegistry("stop") # this function will kill any currently running parallel jobs
future::plan(future::multisession(workers = round(future::availableCores()/2)))

data.frame(file = scanNames) %>%
  # skip pre-processing scans that have already been pre-processed
  filter(!(file %in% list.files("sample_400_preprocessed",recursive = TRUE))) %>%
  furrr::future_pwalk(~ {

    # print(..1)
    preprocessScans_safe(..1)

  })
