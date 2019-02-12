## ---------------------------------------------------------
## Test functions for IBM Hydroperiod modelling
## 
## By: J. Burkhart and B. Ousterhout
## ---------------------------------------------------------

## Load Packages:
## --------------
  if(!require(actuar)) install.packages('actuar', dependencies=T); library('actuar')
  if(!require(raster)) install.packages('raster', dependencies=T); library('raster')
  if(!require(landscapeR)) install.packages('landscapeR', dependencies=T); library('landscapeR')
  if(!require(rgdal)) install.packages('rgdal', dependencies=T); library('rgdal')
## --------------


## Set Directories:
## ----------------
  setwd("C:/Users/Jacob/Git_Repository/IBM_Hydro/")
  input.dir <- paste0(getwd(), "/Data/")
  output.dir <- paste0(getwd(), "/Output/")
  plot.dir <- paste0(output.dir, "Plots/")
## ----------------


## Initialize Data Files: 
## ----------------------
  gen.df <- read.csv(paste0(input.dir, "WEP_Modularity_DataFile.csv"), header=T, stringsAsFactors = F)
  
## ----------------------


## Execute IBM:
## ------------
   ## Breed Function
  
   ## Dipserse Function
  
   ## Age/Grow/Die
  
  
## ------------