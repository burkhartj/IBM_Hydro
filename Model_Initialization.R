##################################
# To Do: Consolidate variables
#        Replace variables in disp model trunc call with flags to this file
##################################


## ---------------------------------------------------------
## Test functions for IBM Hydroperiod modelling
## 
## By: J. Burkhart and B. Ousterhout
## ---------------------------------------------------------

## Load Packages:
## --------------
  if(!require(actuar)) install.packages('actuar', dependencies=T); library('actuar')
  if(!require(adegenet)) install.packages('adegenet'); library("adegenet")
  if(!require(ggplot2)) install.packages('ggplot2'); library("ggplot2")
  if(!require(hierfstat)) install.packages('hierfstat'); library("hierfstat")
  if(!require(landscapeR)) install.packages('landscapeR', dependencies=T); library('landscapeR')
  if(!require(magicfor)) devtools::install_github("hoxo-m/magicfor"); library("magicfor")
  if(!require(mmod)) install.packages('mmod'); library("mmod")
  if(!require(pegas)) install.packages('pegas'); library("pegas")
  if(!require(raster)) install.packages('raster', dependencies=T); library('raster')
  if(!require(rgdal)) install.packages('rgdal', dependencies=T); library('rgdal')
  if(!require(reshape2)) install.packages('reshape2'); library("reshape2")
  if(!require(truncdist)) install.packages('truncdist'); library('truncdist')
  if(!require(plyr)) install.packages('plyr'); library("plyr")

## --------------


## Set Directories:
## ----------------
  setwd("C:/Users/jjbxb3/Git_Repository/IBM_Hydro/")
  input.dir <- paste0(getwd(), "/Data/")
  output.dir <- paste0(getwd(), "/Output/")
  plot.dir <- paste0(output.dir, "Plots/")
## ----------------


## Initialize Data Files: 
## ----------------------
  gen.df <- read.csv(paste0(input.dir, "WEP_Modularity_DataFile.csv"), header=T, stringsAsFactors = F)
  

## Set model parameters: 
## ----------------------

  ## Landscape size:
  ##----------------
  
  landscape.x.rows <- 100
  landscape.y.columns <- 100
  landscape.x.coor.min <- 0
  landscape.x.coor.max <- 3000
  landscape.y.coor.min <- 0
  landscape.y.coor.max <- 3000
  # n.patch <- ls[[1]]@ncols * ls[[1]]@nrows      ## Used to calculate landscape carrying capacity later 
  terrestrial.k <- 180                 ## Terrestrial carrying capacity
  patch.K.mult <- 180             ## Carrying capacity multiplier for each 30x30 m grid cell
  
  ##----------------
  
    
  ## Pond info:
  ##-------------

  n.ponds <- 1                        ## number of ponds to create
  pond.size.mean <- 1                  ## mean size of pond (in 30x30 m grid cells)
  pond.size.sd <- 0                    ## st devation of pond
  pond.hydro.class <- 1                ## number of pond classes (used for hydroperiod later?)
  pond.K.mult <- 2.25             ## Carrying capacity multiplier. Based off a Semlitsch paper 
  
  min.hydro <- 1                   ## minimum hydroperiod class; if < this threshold, no reproduction in pond
  n.hydro.class <- 3               ## number of hydroperiod classes (used in random draws currently, will be used for deriving parameters later rather than hard coding like below)
  hydro.0.mu <- 42.8
  hydro.0.sd <- 5
  hydro.1.mu <- 76.5
  hydro.1.sd <- 10
  hydro.2.mu <- 437.5
  hydro.2.sd <- 50
  hydro.3.mu <- 3003.5
  hydro.3.sd <- 100
  ##-------------

  ##Dispersal:
  ##-------------
  
  sensing.matrix <- matrix(c(NA,1, 1, 1,NA,  ## Can sense 2 rings of cells around current location (opts)
                             1, 1, 1, 1, 1,
                             1, 1, 0, 1, 1, 
                             1, 1, 1, 1, 1,
                             NA,1, 1, 1,NA),
                           ncol=5, byrow=TRUE)
  
  max.disp.dist <- 5000                ## maximum dispersal distance (in meters) for creating a buffer 
                                       ## (quick and dirty fix for world wrapping)
  
  new.move <- 30                       ## Step length
  philo.rate <- 0.90               ## rate of philopatry
  disp.shape <- 1.5                ## shape param for rllogis dispersal dist draws
  disp.scale <- 30                 ## scale param for rllogis dispersal dist draws
  disp.beta1 <- 0.45               ## shape1 for rbeta dispersal dist
  disp.beta2 <- 1.00               ## shape2 for rbeta dispersal dist
  ##--------------
  
  ## Breeding:
  ##-------------
  n.inds <- 100              ## number of individuals to initialize model with
  n.gens <- 100              ## number of generations to iterate over 
  ##--------------
  
  ## Demographic info:
  ##--------------
  
  min.SVL.F <- 48            ## minimum SVL for sexual maturity - females ---- ARBITRARY RIGHT NOW 
  min.SVL.M <- 45            ## minimum SVL for sexual maturity - males ---- ARBITRARY RIGHT NOW 
  SVL.0.mu <- 26.00          ## SVL params based on Taylor and Scott equations, MU's are data driven, SD's ARE ARBITRARY (sizes are for AMOP, maybe use Winters thesis [skeletal chronology] for AMAN size classes?), 
  SVL.0.sd <- 1.5
  SVL.1.mu <- 39.61
  SVL.1.sd <- 0.5
  SVL.2.mu <- 44.47
  SVL.2.sd <- 0.25
  SVL.3.mu <- 48.47 
  SVL.3.sd <- 0.10
  SVL.4p.mu <- 48.97
  SVL.4p.sd <- 0.01
  
  min.age.F <- 2                   ## minimum age for reproduction - females
  min.age.M <- 1                   ## minumum age for reproduction - males
  max.age <- 15                    ## maximum age for all adults (based on winters thesis)
  
  mort.prob.mu <- 0.69             ## probability of ADULT mortality - MU for rnorm draw
  mort.prob.sd <- 0.10             ## probability of ADULT mortality - SD for rnorm draw
  ##--------------

## ----------------------  

## Execute IBM:
## ------------
   ## Initialize landscape Function 
  
   ## Breed Function
  
   ## Dipserse Function
  
   ## Age/Grow/Die
  
  
## ------------