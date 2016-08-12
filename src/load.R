rm(list = ls()[!ls() %in% 'prj'])
### ----------------------------------------------------------------------------
### Setup project structure
## Project Path
# prj <- '/home/edisz/Documents/Uni/Projects/PHD/4BFG/Project' # EDiLD
# prj <- 'C:\\Users\\Edi\\Documents\\BFG\\monitoring_ger' # EDiLD2
if (!exists("prj")) {
  stop("You need to create a object 'prj' that points to the top folder, 
       e.g. prj <- '/home/edisz/Documents/Uni/Projects/PHD/4BFG/Project'!")
}


## set paths
srcdir <- file.path(prj, "src")     # source code
datadir <- file.path(prj, "data")	# raw data
cachedir <- file.path(prj, "cache")	# intermediate data


### ----------------------------------------------------------------------------
### load packages
# data handling
#! use developmental version of data.table
#! CRAN version should also work...
# install.packages('devtools')
library(devtools)
# install_github("Rdatatable/data.table")
library(data.table)
library(xlsx)    # to read .xls files
# Note this needs Java JDK (http://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html)
library(reshape2)
library(plyr)
library(bit64)
library(Hmisc)
library(stringr)
# install.packages("devtools")
# devtools::install_github("hadley/readxl")
library(readxl)
library(stringi)
library(stringdist)

# DB
library(RPostgreSQL)

# spatial
library(sp)
library(rgdal)
library(maptools)
library(rgeos)
library(raster)
library(ggmap)

# plotting
library(ggplot2)
library(gridExtra)
library(ggExtra)
library(cowplot)
library(viridis)

# misc
library(pbapply)
library(esmisc) # https://github.com/EDiLD/esmisc
library(webchem)
library(xtable)

# modelling
library(vegan)
library(mgcv)
library(gamlss)
library(gamlss.cens)
library(gamlss.dist)
library(gamlss.add)


### ----------------------------------------------------------------------------
### Source defined functions
source(file.path(srcdir, 'functions.R'))


### ----------------------------------------------------------------------------
# connection details to database
source(file.path(srcdir, 'credentials.R'))




### ----------------------------------------------------------------------------
### switches
# should file be written to cache (as .csv)?
tocache <- TRUE  
options(scipen = 999)

# object to check if load.R was run
ld <- TRUE

# keep these objects
load_obj <- ls()
load_obj <- c(load_obj, 'load_obj')

# use utf8 encoding
options(encoding = 'utf8')
options(stringsAsFactors = FALSE)
# check OS
os <- Sys.info()['sysname']
# use english time
if (os == 'windows')
  Sys.setlocale("LC_TIME", "English")

