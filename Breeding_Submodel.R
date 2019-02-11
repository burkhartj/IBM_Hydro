## ---------------------------------------------------------
## Test functions for IBM Hydroperiod modelling
## 
## By: J. Burkhart and B. Ousterhout
## ---------------------------------------------------------

## Load Packages:
## --------------
  if(!require(actuar)) install.packages('actuar'); library('actuar')
  if(!require(adegenet)) install.packages('adegenet'); library("adegenet")
  if(!require(pegas)) install.packages('pegas'); library("pegas")
  if(!require(hierfstat)) install.packages('hierfstat'); library("hierfstat")
  if(!require(mmod)) install.packages('mmod'); library("mmod")
  if(!require(reshape2)) install.packages('reshape2'); library("reshape2")
  if(!require(ggplot2)) install.packages('ggplot2'); library("ggplot2")
## --------------

## Import Data Files: 
## ------------------
  gen.df <- read.csv(paste0(input.dir, "WEP_Modularity_DataFile.csv"), header=T, stringsAsFactors = F)        ## import genetic data (move to the model_initialization.R script later)
    gen.df <- gen.df[,-((dim(gen.df)[2]-2):dim(gen.df)[2])]                                                   ## trim up the loci that were ommitted
## ------------------

## Initialize Models:
## ------------------
  ## Input Parameters to automate changes:
    n.inds <- 100               ## number of individuals to create across all ponds
    n.ponds <- 5               ## number of initial ponds to create
    n.patch <- ls[[1]]@ncols * ls[[1]]@nrows            ## total number of patches ---- TEMPORARY, DELETE WHEN LANDSCAPE UPDATE WORKS
    n.gens <- 100              ## number of generations to iterate over 
    
    pond.K.mult <- 2.25             ## Carrying capacity multiplier. Based off a Semlitsch paper 
    patch.K.mult <- 180             ## Carrying capacity multiplier for each 30x30 m grid cell
    
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
    max.age <- 15                    ## maximum age for all adults
    
    max.disp.dist <- 5000            ## maximum possible dispersal distance (used for landscape border buffer)
    philo.rate <- 0.90               ## rate of philopatry
    disp.shape <- 1.5                ## shape param for rllogis dispersal dist draws
    disp.scale <- 30                 ## scale param for rllogis dispersal dist draws
    disp.beta1 <- 0.45               ## shape1 for rbeta dispersal dist
    disp.beta2 <- 1.00               ## shape2 for rbeta dispersal dist
    mort.prob.mu <- 0.69             ## probability of ADULT mortality - MU for rnorm draw
    mort.prob.sd <- 0.10             ## probability of ADULT mortality - SD for rnorm draw
    
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
    # hydro.4.mu <- 
    # hydro.4.sd <-   
    
    temp.terrestrial.K <- n.patch * patch.K.mult           ## DELETE LATER. WILL BE IRRELEVANT ONCE THE SPATIAL STUFF IS INCORPORATED
    
    
  ## Initialize Pond Data Frame:
  ponds <- data.frame(Pond.ID = seq(1:n.ponds),
                      Pond.X = numeric(n.ponds), 
                      Pond.Y = numeric(n.ponds),
                      Hydroperiod = sample(0:n.hydro.class, n.ponds, T), 
                      Pond.Area = numeric(n.ponds), 
                      Pond.K = numeric(n.ponds), 
                      N.inds = numeric(n.ponds), 
                      Mort.Pond = numeric(n.ponds)
                      )
  
  for(i in 1:dim(ponds)[1]){
    ponds$Pond.Area[i] <- ifelse(ponds$Hydroperiod[i] == 0, yes=round(rnorm(1, hydro.0.mu, hydro.0.sd), 3), 
                                 no=ifelse(ponds$Hydroperiod[i] == 1, yes=round(rnorm(1, hydro.1.mu, hydro.1.sd), 3),
                                           no=ifelse(ponds$Hydroperiod[i] == 2, yes=round(rnorm(1, hydro.2.mu, hydro.2.sd), 3), 
                                                     no=ifelse(ponds$Hydroperiod[i] == 3, yes=round(rnorm(1, hydro.3.mu, hydro.3.sd), 3), 
                                                               no=NA))))
  }
  
  ponds$Pond.K <- round(pond.K.mult * ponds$Pond.Area)
  
  ponds$Pond.X <- as.data.frame(rasterToPoints(pond.r, function(x){x == 1}))$x
  ponds$Pond.Y <- as.data.frame(rasterToPoints(pond.r, function(x){x == 1}))$y
          
  ## Create Data Frame of Individuals:  inds-own [age, sex, svl, loci, my-land-home (Nat.Patch), my-natal-pond (Nat.Pond),  my-breed-pond (Breed.Pond), can-breed (Rep.Active), years-since-breed (IBI)]
  inds <- data.frame(Sex = sample(c("M", "F"), n.inds, replace=T), 
                     Age = sample(0:round(max.age / 2), n.inds, T), 
                     SVL = numeric(n.inds), 
                     Rep.Active = NA, 
                     IBI = numeric(n.inds), 
                     Bred = 0, 
                     Generation = 0, 
                     Nat.Pond = sample(1:n.ponds, n.inds, T), 
                     Patch.X = numeric(n.inds), 
                     Patch.Y = numeric(n.inds),
                     Disp.Prob = round(rbeta(n.inds, disp.beta1, disp.beta2), 3), 
                     Breed.Pond = numeric(n.inds),
                     Init.Angle = numeric(n.inds), 
                     Disp.Dist = numeric(n.inds), 
                     Mort.Prob = numeric(n.inds), 
                     LocA=numeric(n.inds), LocB=numeric(n.inds), 
                     LocC=numeric(n.inds), LocD=numeric(n.inds), 
                     LocE=numeric(n.inds), LocF=numeric(n.inds),
                     LocG=numeric(n.inds), LocH=numeric(n.inds),
                     LocI=numeric(n.inds), LocJ=numeric(n.inds),
                     LocK=numeric(n.inds), LocL=numeric(n.inds),
                     LocM=numeric(n.inds), LocN=numeric(n.inds),
                     LocO=numeric(n.inds), LocP=numeric(n.inds),
                     LocQ=numeric(n.inds), LocR=numeric(n.inds),
                     LocS=numeric(n.inds), LocT=numeric(n.inds),
                     LocU=numeric(n.inds)
                     )
  
     inds$Breed.Pond <- ifelse(inds$Disp.Prob > philo.rate, yes=sample(unique(inds$Nat.Pond), 1, T), no=inds$Nat.Pond)
     
     ponds$N.inds <- as.data.frame(table(inds$Nat.Pond))$Freq
     
     
    ## Add Genetic Data 
    for(i in 1:dim(inds)[1]){
      ## Calculate SVL:              ## automate these steps later
      inds$SVL[i] <- ifelse(inds$Age[i] == 0, yes=round(rnorm(1, SVL.0.mu, SVL.0.sd), 2),
                         no=ifelse(inds$Age[i] == 1, yes=round(rnorm(1, SVL.1.mu, SVL.1.sd), 2), 
                                   no=ifelse(inds$Age[i] == 2, yes=round(rnorm(1, SVL.2.mu, SVL.2.sd), 2), 
                                             no=ifelse(inds$Age[i] == 3, yes=round(rnorm(1, SVL.3.mu, SVL.3.sd), 2), 
                                                       no=ifelse(inds$Age[i] >= 4, yes=round(rnorm(1, SVL.4p.mu, SVL.4p.sd) + (inds$Age[i] - 4) * 0.5, 2), no=NA)))))
     
      ## Initialize the Reproductive States:
      inds$Rep.Active[i] <- ifelse((inds$Sex[i]=="F" & inds$Age[i]<min.age.F) | (inds$Sex[i]=="F" & inds$SVL[i]<min.SVL.F), 
                                yes = F, 
                                no = ifelse(inds$Sex[i]=="F" & inds$Age[i]>=min.age.F & inds$SVL[i]>=min.SVL.F, 
                                         yes=sample(x = c(T, F), size = 1, replace=T),
                                         no = ifelse((inds$Sex[i]=="M" & inds$Age[i]<min.age.M) | (inds$Sex[i]=="M" & inds$SVL[i]<min.SVL.M), 
                                                  yes = F,
                                                  no = ifelse((inds$Sex[i]=="M" & inds$Age[i]>=min.age.M & inds$SVL[i]>=min.SVL.M),
                                                           yes = T, 
                                                           no = NA))))
      
      ## Initialize Genetic Data:
          ## Try to automate this later. Can loop over columns if you tell where the columns start and where they end...  
      inds$LocA[i] <- paste0(sample(unlist(strsplit(gen.df$Aa_37, ":")), 1), ":", sample(unlist(strsplit(gen.df$Aa_37, ":")), 1))
      inds$LocB[i] <- paste0(sample(unlist(strsplit(gen.df$Aa_50, ":")), 1), ":", sample(unlist(strsplit(gen.df$Aa_50, ":")), 1))
      inds$LocC[i] <- paste0(sample(unlist(strsplit(gen.df$ac300, ":")), 1), ":", sample(unlist(strsplit(gen.df$ac300, ":")), 1))
      inds$LocD[i] <- paste0(sample(unlist(strsplit(gen.df$Aa_25, ":")), 1), ":", sample(unlist(strsplit(gen.df$Aa_25, ":")), 1))
      inds$LocE[i] <- paste0(sample(unlist(strsplit(gen.df$Aa_258, ":")), 1), ":", sample(unlist(strsplit(gen.df$Aa_258, ":")), 1))
      inds$LocF[i] <- paste0(sample(unlist(strsplit(gen.df$Aa_85, ":")), 1), ":", sample(unlist(strsplit(gen.df$Aa_85, ":")), 1))
      inds$LocG[i] <- paste0(sample(unlist(strsplit(gen.df$Aa_44, ":")), 1), ":", sample(unlist(strsplit(gen.df$Aa_44, ":")), 1))
      inds$LocH[i] <- paste0(sample(unlist(strsplit(gen.df$Aa_39, ":")), 1), ":", sample(unlist(strsplit(gen.df$Aa_39, ":")), 1))
      inds$LocI[i] <- paste0(sample(unlist(strsplit(gen.df$Aa_40, ":")), 1), ":", sample(unlist(strsplit(gen.df$Aa_40, ":")), 1))
      inds$LocJ[i] <- paste0(sample(unlist(strsplit(gen.df$Aa_312, ":")), 1), ":", sample(unlist(strsplit(gen.df$Aa_312, ":")), 1))
      inds$LocK[i] <- paste0(sample(unlist(strsplit(gen.df$aj_346, ":")), 1), ":", sample(unlist(strsplit(gen.df$aj_346, ":")), 1))
      inds$LocL[i] <- paste0(sample(unlist(strsplit(gen.df$Aa_31, ":")), 1), ":", sample(unlist(strsplit(gen.df$Aa_31, ":")), 1))
      inds$LocM[i] <- paste0(sample(unlist(strsplit(gen.df$Aa_311, ":")), 1), ":", sample(unlist(strsplit(gen.df$Aa_311, ":")), 1))
      inds$LocN[i] <- paste0(sample(unlist(strsplit(gen.df$Aa_36, ":")), 1), ":", sample(unlist(strsplit(gen.df$Aa_36, ":")), 1))
      inds$LocO[i] <- paste0(sample(unlist(strsplit(gen.df$Aa_21, ":")), 1), ":", sample(unlist(strsplit(gen.df$Aa_21, ":")), 1))
      inds$LocP[i] <- paste0(sample(unlist(strsplit(gen.df$Aa_27, ":")), 1), ":", sample(unlist(strsplit(gen.df$Aa_27, ":")), 1))
      inds$LocQ[i] <- paste0(sample(unlist(strsplit(gen.df$Aa_28, ":")), 1), ":", sample(unlist(strsplit(gen.df$Aa_28, ":")), 1))
      inds$LocR[i] <- paste0(sample(unlist(strsplit(gen.df$Aa_86, ":")), 1), ":", sample(unlist(strsplit(gen.df$Aa_86, ":")), 1))
      inds$LocS[i] <- paste0(sample(unlist(strsplit(gen.df$Aa_153, ":")), 1), ":", sample(unlist(strsplit(gen.df$Aa_153, ":")), 1))
      inds$LocT[i] <- paste0(sample(unlist(strsplit(gen.df$Aa_84, ":")), 1), ":", sample(unlist(strsplit(gen.df$Aa_84, ":")), 1))
      inds$LocU[i] <- paste0(sample(unlist(strsplit(gen.df$Aa_314, ":")), 1), ":", sample(unlist(strsplit(gen.df$Aa_314, ":")), 1))
      
      
      ## Intialize Dispersal Data: 
      inds$Init.Angle[i] <- round(runif(1, 1, 360)) 
      inds$Disp.Dist[i] <- rllogis(1, shape=disp.shape, scale=disp.scale)
      inds$Patch.X[i] <- ponds[ponds$Pond.ID %in% inds$Nat.Pond[i], "Pond.X"] + inds$Disp.Dist[i] * sin((inds$Init.Angle[i])*(pi/180))
      inds$Patch.Y[i] <- ponds[ponds$Pond.ID %in% inds$Nat.Pond[i], "Pond.Y"] + inds$Disp.Dist[i] * cos((inds$Init.Angle[i])*(pi/180))

      
        repeat{
          Patch.K <- extract(ls[[2]], cbind(inds$Patch.X[i], inds$Patch.Y[i]))    ## extract patch 
          patch.inds <- dim(subset(inds, Patch.X == inds$Patch.X[i] & Patch.Y == inds$Patch.Y[i]))[1]   ## calculate # inds on patch Resident surface and update that. 
          
          print(paste0("Check Terr K for ind #", i, " ----- Patch K = ", Patch.K, " ----- Patch.inds = ", patch.inds))
          
          if(patch.inds > floor(Patch.K / 2) | is.na(Patch.K)==T) {
            # print(paste0("Old Angle = ", round(inds$Init.Angle[i], 3), 
            #              " ----- Old Disp Dist = ", round(inds$Disp.Dist[i], 3), 
            #              " ----- Old XY-Coor = (", floor(inds$Patch.X[i]), ", ",
            #              floor(inds$Patch.Y[i]), ")"))
            
            inds$Init.Angle[i] <- round(runif(1, 1, 360)) 
            inds$Disp.Dist[i] <- rllogis(1, shape=disp.shape, scale=disp.scale)
            inds$Patch.X[i] <- ponds[ponds$Pond.ID %in% inds$Nat.Pond[i], "Pond.X"] + 
                                  inds$Disp.Dist[i] * sin((inds$Init.Angle[i])*(pi/180))
            inds$Patch.Y[i] <- ponds[ponds$Pond.ID %in% inds$Nat.Pond[i], "Pond.Y"] + 
                                  inds$Disp.Dist[i] * cos((inds$Init.Angle[i])*(pi/180))
          }
            else {
              # print("Patch not full, next ind")
              break 
            }      ## end else statement
        }      ## end repeat loop
    }      ## end for loop
   
  

## ------------------  

## Breeding Migration:
## -------------------
# plot(x=seq(from=0, to=max.age+1, length.out = max.age+1), 
#      y=seq(from=0, to=max(ponds$Pond.K), length.out = max.age+1), 
#      pch=16, col="white", xlab="Age (years)", ylab="# Individuals")

inds0 <- inds    

start.time <- Sys.time()       
for(g in 1:n.gens){    
  ## Set Up Breeding Functionality
  p.sub <- subset(ponds, Hydroperiod >= min.hydro)   ## Subset out the pond with non-suitable hydroperiods
  pond.list <- unique(p.sub$Pond.ID)                 ## Create a vector of unique ID's to iterate over.  CAN BE DELETED IN THE FUTURE; however, 
  off.df <- inds[-c(1:dim(inds)[1]), ]               ## create empty data frame to populate with offspring within the loop
  
  
  ## Pond Test Condition
  
  ## Iterate through suitable breeding ponds
  for(p in 1:length(pond.list)){
    off.pond <- inds[-c(1:dim(inds)[1]), ]             ## create empty data frame for offspring within each pond
    rep.feme <- subset(inds, Sex == "F" & Rep.Active == T & Age >= min.age.F & 
                         Breed.Pond == pond.list[p] & SVL >= min.SVL.F)
    rep.male <- subset(inds, Sex == "M" & Rep.Active == T & Age >= min.age.M & 
                         Breed.Pond == pond.list[p] & SVL >= min.SVL.M)
    
    ## Iterate through reproductive females at each pond
    if(dim(rep.feme)[1] > 0 & dim(rep.male)[1] > 0){             ## only execute if there is at least 1 reproductive male and 1 reproductive female
      
      inds$Bred <- ifelse(((inds$Sex=="F" & inds$Rep.Active==T & inds$Breed.Pond==pond.list[p] &
                              inds$SVL>=min.SVL.F & inds$Age>=min.age.F) |
                           (inds$Sex=="M" & inds$Rep.Active==T & inds$Breed.Pond==pond.list[p] & 
                              inds$SVL>=min.SVL.M & inds$Age>=min.age.M)),
                          yes = 1, no = inds$Bred)
      
      for(f in 1:dim(rep.feme)[1]){   
        n.off <- round((8.47 * rep.feme$SVL[f] - 380) * abs(rnorm(1, 0.355, 0.11)))      ## Size based fecundity * survival to age 1
        n.off <- ifelse(n.off < 0, yes=0, no=n.off)
      
        male.df <- rep.male[sample(1:dim(rep.male)[1], n.off, replace=T), ]
        
        if(n.off > 0){  
          num.off <- as.data.frame(matrix(nrow = n.off, ncol=dim(inds)[2]))                ## create empty data frame to populate with an individual female's offspring during breeding loop
          colnames(num.off) <- colnames(inds)
          
          ## Demographics
          num.off$Sex <- sample(c("F", "M"), n.off, T)
          num.off$Age <- rep(0, times=n.off)
          num.off$SVL <- round(rnorm(n.off, SVL.0.mu, SVL.0.sd), 2)            ## random normal draw for now. will update to density dependent
          num.off$Rep.Active <- rep(F, times=n.off)                            ## Set reproductively active to "FALSE"
          num.off$IBI <- rep(0, times=n.off)                                   ## Set interbreeding interval to 0
          num.off$Nat.Pond <- rep(pond.list[p], times=n.off)                   ## Set natal pond to mother's breeding pond
          num.off$Patch.X <- numeric(n.off)                                   ## patch x-coordinate 
          num.off$Patch.Y <- numeric(n.off)                                   ## patch y-coordinate
          num.off$Disp.Prob <- round(rbeta(n.off, disp.beta1, disp.beta2), 3)  ## create a vector of dispersal probability
          num.off$Breed.Pond <- ifelse(num.off$Disp.Prob > philo.rate, yes=sample(unique(ponds$Pond.ID), 1, T), no=num.off$Nat.Pond)
          num.off$Init.Angle <- round(runif(n.off, 1, 360))                    ## choose initial movement angle
          num.off$Disp.Dist <- rllogis(1, shape=disp.shape, scale=disp.scale)  ## Peterman et al. 2015 dispersal distance for A. annulatum? 1,693m (95% CI 1,645-1,740 m for AMAN)
          num.off$Mort.Prob <- numeric(n.off)                                  ## empty vector for later imposing mortality 
          num.off$Bred <- rep(0, times=n.off)                                  ## vector of zeros for breeding
          num.off$Generation <- rep(g, times=n.off)                            ## vector to store generation of origin 
          
          ## Genetics --
            ## System mimics polyandry. Randomly selects one of two females alleles, equal to number of offspring. 
            ## Randomly selects on of males alleles at each locus, each offspring has single paternity.
<<<<<<< HEAD
            num.off$LocA <- paste0(unlist(strsplit(rep.feme$LocA[f], ":"))[sample(1:2, n.off, T)], ":", unlist(strsplit(male.df$LocA, ":"))[sample(1:2, n.off, T) + seq(0, n.off*2-1, by=2)])
            num.off$LocB <- paste0(unlist(strsplit(rep.feme$LocB[f], ":"))[sample(1:2, n.off, T)], ":", unlist(strsplit(male.df$LocB, ":"))[sample(1:2, n.off, T) + seq(0, n.off*2-1, by=2)])
            num.off$LocC <- paste0(unlist(strsplit(rep.feme$LocC[f], ":"))[sample(1:2, n.off, T)], ":", unlist(strsplit(male.df$LocC, ":"))[sample(1:2, n.off, T) + seq(0, n.off*2-1, by=2)])
            num.off$LocD <- paste0(unlist(strsplit(rep.feme$LocD[f], ":"))[sample(1:2, n.off, T)], ":", unlist(strsplit(male.df$LocD, ":"))[sample(1:2, n.off, T) + seq(0, n.off*2-1, by=2)])
            num.off$LocE <- paste0(unlist(strsplit(rep.feme$LocE[f], ":"))[sample(1:2, n.off, T)], ":", unlist(strsplit(male.df$LocE, ":"))[sample(1:2, n.off, T) + seq(0, n.off*2-1, by=2)])
            num.off$LocF <- paste0(unlist(strsplit(rep.feme$LocF[f], ":"))[sample(1:2, n.off, T)], ":", unlist(strsplit(male.df$LocF, ":"))[sample(1:2, n.off, T) + seq(0, n.off*2-1, by=2)])
            num.off$LocG <- paste0(unlist(strsplit(rep.feme$LocG[f], ":"))[sample(1:2, n.off, T)], ":", unlist(strsplit(male.df$LocG, ":"))[sample(1:2, n.off, T) + seq(0, n.off*2-1, by=2)])
            num.off$LocH <- paste0(unlist(strsplit(rep.feme$LocH[f], ":"))[sample(1:2, n.off, T)], ":", unlist(strsplit(male.df$LocH, ":"))[sample(1:2, n.off, T) + seq(0, n.off*2-1, by=2)])
            num.off$LocI <- paste0(unlist(strsplit(rep.feme$LocI[f], ":"))[sample(1:2, n.off, T)], ":", unlist(strsplit(male.df$LocI, ":"))[sample(1:2, n.off, T) + seq(0, n.off*2-1, by=2)])
            num.off$LocJ <- paste0(unlist(strsplit(rep.feme$LocJ[f], ":"))[sample(1:2, n.off, T)], ":", unlist(strsplit(male.df$LocJ, ":"))[sample(1:2, n.off, T) + seq(0, n.off*2-1, by=2)])
            num.off$LocK <- paste0(unlist(strsplit(rep.feme$LocK[f], ":"))[sample(1:2, n.off, T)], ":", unlist(strsplit(male.df$LocK, ":"))[sample(1:2, n.off, T) + seq(0, n.off*2-1, by=2)])
            num.off$LocL <- paste0(unlist(strsplit(rep.feme$LocL[f], ":"))[sample(1:2, n.off, T)], ":", unlist(strsplit(male.df$LocL, ":"))[sample(1:2, n.off, T) + seq(0, n.off*2-1, by=2)])
            num.off$LocM <- paste0(unlist(strsplit(rep.feme$LocM[f], ":"))[sample(1:2, n.off, T)], ":", unlist(strsplit(male.df$LocM, ":"))[sample(1:2, n.off, T) + seq(0, n.off*2-1, by=2)])
            num.off$LocN <- paste0(unlist(strsplit(rep.feme$LocN[f], ":"))[sample(1:2, n.off, T)], ":", unlist(strsplit(male.df$LocN, ":"))[sample(1:2, n.off, T) + seq(0, n.off*2-1, by=2)])
            num.off$LocO <- paste0(unlist(strsplit(rep.feme$LocO[f], ":"))[sample(1:2, n.off, T)], ":", unlist(strsplit(male.df$LocO, ":"))[sample(1:2, n.off, T) + seq(0, n.off*2-1, by=2)])
            num.off$LocP <- paste0(unlist(strsplit(rep.feme$LocP[f], ":"))[sample(1:2, n.off, T)], ":", unlist(strsplit(male.df$LocP, ":"))[sample(1:2, n.off, T) + seq(0, n.off*2-1, by=2)])
            num.off$LocQ <- paste0(unlist(strsplit(rep.feme$LocQ[f], ":"))[sample(1:2, n.off, T)], ":", unlist(strsplit(male.df$LocQ, ":"))[sample(1:2, n.off, T) + seq(0, n.off*2-1, by=2)])
            num.off$LocR <- paste0(unlist(strsplit(rep.feme$LocR[f], ":"))[sample(1:2, n.off, T)], ":", unlist(strsplit(male.df$LocR, ":"))[sample(1:2, n.off, T) + seq(0, n.off*2-1, by=2)])
            num.off$LocS <- paste0(unlist(strsplit(rep.feme$LocS[f], ":"))[sample(1:2, n.off, T)], ":", unlist(strsplit(male.df$LocS, ":"))[sample(1:2, n.off, T) + seq(0, n.off*2-1, by=2)])
            num.off$LocT <- paste0(unlist(strsplit(rep.feme$LocT[f], ":"))[sample(1:2, n.off, T)], ":", unlist(strsplit(male.df$LocT, ":"))[sample(1:2, n.off, T) + seq(0, n.off*2-1, by=2)])
            num.off$LocU <- paste0(unlist(strsplit(rep.feme$LocU[f], ":"))[sample(1:2, n.off, T)], ":", unlist(strsplit(male.df$LocU, ":"))[sample(1:2, n.off, T) + seq(0, n.off*2-1, by=2)])
=======
            num.off$LocA <- paste0(unlist(strsplit(rep.feme$LocA[f], ":"))[sample(1:2, n.off, T)], ":",  
                                   unlist(strsplit(male.df$LocA, ":"))[sample(1:2, n.off, T) + seq(0, n.off*2-1, by=2)])
            num.off$LocB <- paste0(unlist(strsplit(rep.feme$LocB[f], ":"))[sample(1:2, n.off, T)], ":",  
                                   unlist(strsplit(male.df$LocB, ":"))[sample(1:2, n.off, T) + seq(0, n.off*2-1, by=2)])
            num.off$LocC <- paste0(unlist(strsplit(rep.feme$LocC[f], ":"))[sample(1:2, n.off, T)], ":",  
                                   unlist(strsplit(male.df$LocC, ":"))[sample(1:2, n.off, T) + seq(0, n.off*2-1, by=2)])
            num.off$LocD <- paste0(unlist(strsplit(rep.feme$LocD[f], ":"))[sample(1:2, n.off, T)], ":",  
                                   unlist(strsplit(male.df$LocD, ":"))[sample(1:2, n.off, T) + seq(0, n.off*2-1, by=2)])
            num.off$LocE <- paste0(unlist(strsplit(rep.feme$LocE[f], ":"))[sample(1:2, n.off, T)], ":",  
                                   unlist(strsplit(male.df$LocE, ":"))[sample(1:2, n.off, T) + seq(0, n.off*2-1, by=2)])
            num.off$LocF <- paste0(unlist(strsplit(rep.feme$LocF[f], ":"))[sample(1:2, n.off, T)], ":",  
                                   unlist(strsplit(male.df$LocF, ":"))[sample(1:2, n.off, T) + seq(0, n.off*2-1, by=2)])
            num.off$LocG <- paste0(unlist(strsplit(rep.feme$LocG[f], ":"))[sample(1:2, n.off, T)], ":",  
                                   unlist(strsplit(male.df$LocG, ":"))[sample(1:2, n.off, T) + seq(0, n.off*2-1, by=2)])
            num.off$LocH <- paste0(unlist(strsplit(rep.feme$LocH[f], ":"))[sample(1:2, n.off, T)], ":",  
                                   unlist(strsplit(male.df$LocH, ":"))[sample(1:2, n.off, T) + seq(0, n.off*2-1, by=2)])
            num.off$LocI <- paste0(unlist(strsplit(rep.feme$LocI[f], ":"))[sample(1:2, n.off, T)], ":",  
                                   unlist(strsplit(male.df$LocI, ":"))[sample(1:2, n.off, T) + seq(0, n.off*2-1, by=2)])
            num.off$LocJ <- paste0(unlist(strsplit(rep.feme$LocJ[f], ":"))[sample(1:2, n.off, T)], ":",  
                                   unlist(strsplit(male.df$LocJ, ":"))[sample(1:2, n.off, T) + seq(0, n.off*2-1, by=2)])
            num.off$LocK <- paste0(unlist(strsplit(rep.feme$LocK[f], ":"))[sample(1:2, n.off, T)], ":",  
                                   unlist(strsplit(male.df$LocK, ":"))[sample(1:2, n.off, T) + seq(0, n.off*2-1, by=2)])
            num.off$LocL <- paste0(unlist(strsplit(rep.feme$LocL[f], ":"))[sample(1:2, n.off, T)], ":",  
                                   unlist(strsplit(male.df$LocL, ":"))[sample(1:2, n.off, T) + seq(0, n.off*2-1, by=2)])
            num.off$LocM <- paste0(unlist(strsplit(rep.feme$LocM[f], ":"))[sample(1:2, n.off, T)], ":",  
                                   unlist(strsplit(male.df$LocM, ":"))[sample(1:2, n.off, T) + seq(0, n.off*2-1, by=2)])
            num.off$LocN <- paste0(unlist(strsplit(rep.feme$LocN[f], ":"))[sample(1:2, n.off, T)], ":",  
                                   unlist(strsplit(male.df$LocN, ":"))[sample(1:2, n.off, T) + seq(0, n.off*2-1, by=2)])
            num.off$LocO <- paste0(unlist(strsplit(rep.feme$LocO[f], ":"))[sample(1:2, n.off, T)], ":",  
                                   unlist(strsplit(male.df$LocO, ":"))[sample(1:2, n.off, T) + seq(0, n.off*2-1, by=2)])
            num.off$LocP <- paste0(unlist(strsplit(rep.feme$LocP[f], ":"))[sample(1:2, n.off, T)], ":",  
                                   unlist(strsplit(male.df$LocP, ":"))[sample(1:2, n.off, T) + seq(0, n.off*2-1, by=2)])
            num.off$LocQ <- paste0(unlist(strsplit(rep.feme$LocQ[f], ":"))[sample(1:2, n.off, T)], ":",  
                                   unlist(strsplit(male.df$LocQ, ":"))[sample(1:2, n.off, T) + seq(0, n.off*2-1, by=2)])
            num.off$LocR <- paste0(unlist(strsplit(rep.feme$LocR[f], ":"))[sample(1:2, n.off, T)], ":",  
                                   unlist(strsplit(male.df$LocR, ":"))[sample(1:2, n.off, T) + seq(0, n.off*2-1, by=2)])
            num.off$LocS <- paste0(unlist(strsplit(rep.feme$LocS[f], ":"))[sample(1:2, n.off, T)], ":",  
                                   unlist(strsplit(male.df$LocS, ":"))[sample(1:2, n.off, T) + seq(0, n.off*2-1, by=2)])
            num.off$LocT <- paste0(unlist(strsplit(rep.feme$LocT[f], ":"))[sample(1:2, n.off, T)], ":",  
                                   unlist(strsplit(male.df$LocT, ":"))[sample(1:2, n.off, T) + seq(0, n.off*2-1, by=2)])
            num.off$LocU <- paste0(unlist(strsplit(rep.feme$LocU[f], ":"))[sample(1:2, n.off, T)], ":",  
                                   unlist(strsplit(male.df$LocU, ":"))[sample(1:2, n.off, T) + seq(0, n.off*2-1, by=2)])
>>>>>>> eb5ba32de58bf98a729efbb6b97bbb5a1f7f9bcc
        }
        
        ## Update output data frame
        off.pond <- rbind(off.pond, num.off)
        
      }  ## end FEMALES iterations
      
      ## Impose LARVAL carrying capacity catch.  
      if(ponds[ponds$Pond.ID %in% pond.list[p], "Pond.K"] < dim(off.pond)[1]){
         off.pond <- off.pond[-sample(x = 1:dim(off.pond)[1], replace=F,
                                      size = (dim(off.pond)[1]) - ponds[ponds$Pond.ID %in% pond.list[p], "Pond.K"]), ]
      }
      
      ## Create the offspring data set... might be able to pull straight into the 'inds' data frame... 
      off.df <- rbind(off.df, off.pond)
      
    } ## end IF statement
  
  } ## end POND iterations
  
  ## Join new and old data frames 
  inds <- rbind(inds, off.df)
  
  print(table(inds$Breed.Pond))
  
## Grow, Mature, Death:
## --------------------
  ## Update demographic information 
  inds$Age <- inds$Age + 1
  inds$SVL <- ifelse(inds$Age == 1, yes=0.785 * inds$SVL + 19.9,  
                     no=ifelse(inds$Age == 2, yes=0.937 * inds$SVL + 7.36, 
                               no=ifelse(inds$Age == 3, yes=inds$SVL + 2.5, 
                                         no=ifelse(inds$Age == 4, yes=inds$SVL + 1.5, 
                                                   no=ifelse(inds$Age >=5, yes=inds$SVL + (inds$Age - 4) * 0.5, no=NA)))))
  
  inds$IBI <- ifelse(inds$Bred == 1, yes = 0, no = inds$IBI + 1)
  inds$Rep.Active <- ifelse((inds$Bred==1 & inds$Sex=="F"),
                            yes=F,
                            no=ifelse((inds$Sex=="F" & inds$Age>=min.age.F & inds$SVL>=min.SVL.F), 
                                      yes = sample(c(T, F), dim(subset(inds, Sex=="F" & Age>=min.age.F & SVL>=min.SVL.F)[1]), replace=T), 
                                      no = ifelse((inds$Sex=="M" & inds$SVL>=min.SVL.M & inds$Age>=min.age.M), 
                                                   yes=T, no=inds$Rep.Active)))
  
  inds$Bred <- 0                                ## reset breeding counter
  
  
  ## Impose Age based mortality
  inds <- inds[which(inds$Age <= max.age), ]                                              ## keep only individuals less than the hard upper age cut off
  inds$Mort.Prob <- runif(dim(inds)[1], 0, 1) > rnorm(dim(inds)[1], mort.prob.mu, mort.prob.sd)      ## create T/F vector for imposing mortality
  inds <- inds[which(inds$Mort.Prob == F), ]                                         ## keep only individuals that pass the random mortality catch
  
  ## Impose TEMPORARY terrestrial carrying capacity
  if(temp.terrestrial.K < dim(inds)[1]) {
    inds <- inds[-sample(x = 1:dim(inds)[1], replace=F,
                         size = (dim(inds)[1]) - temp.terrestrial.K), ]
  }
  
    
  ## Update count of individuals per pond
  # find.ponds <- p.test$Pond.ID %in% as.numeric(as.character(unique(as.data.frame(table(inds$Breed.Pond))$Var1)))
  ponds[ponds$Pond.ID %in% as.numeric(as.character(unique(as.data.frame(table(inds$Breed.Pond))$Var1))), "N.inds"] <- as.data.frame(table(inds$Breed.Pond))$Freq
  ponds[!ponds$Pond.ID %in% as.numeric(as.character(unique(as.data.frame(table(inds$Breed.Pond))$Var1))), "N.inds"] <- 0
  print(ponds) 


  # print(lines(x=density(inds$Age)$x, 
  #            y=density(inds$Age)$y * length(inds$Age), type="l"))
  print(paste0("Generation = ", g, " --- Total # Inds = ", dim(inds)[1], "; Terrestrial K = ", temp.terrestrial.K))
  print(round(Sys.time() - start.time, 2))
} ## end generations loop
print(round(Sys.time() - start.time, 2))         ## end timer

  
## --------------------
  

## Save Genetic Data:
## ------------------
  ## Manipulate data? Export Genepop from gstudio?
    g.df <- inds
    g.df$Ind.ID <- paste0("N", g.df$Nat.Pond, "B", g.df$Breed.Pond, "-", g.df$Sex, "-", rownames(g.df))   
    g.df <- gen.df[,c("Ind.ID", "Breed.Pond", 
                      "LocA", "LocB", "LocC", "LocD", "LocE", "LocF", "LocG", "LocH", "LocI", "LocJ", "LocK", 
                      "LocL", "LocM", "LocN", "LocO", "LocP", "LocQ", "LocR", "LocS", "LocT", "LocU")]

    gen.data <- df2genind(X=g.df[,c("LocA", "LocB", "LocC", "LocD", "LocE", "LocF", "LocG", "LocH", "LocI", "LocJ", "LocK", 
                                    "LocL", "LocM", "LocN", "LocO", "LocP", "LocQ", "LocR", "LocS", "LocT", "LocU")],
                          ind.names=g.df$Ind.ID,
                          loc.names=c("LocA", "LocB", "LocC", "LocD", "LocE", "LocF", "LocG", "LocH", "LocI", "LocJ", "LocK", 
                                      "LocL", "LocM", "LocN", "LocO", "LocP", "LocQ", "LocR", "LocS", "LocT", "LocU"),
                          type="codom",
                          # strata=gen.df$Breed.Pond,
                          pop=g.df$Breed.Pond,
                          sep=":",
                          ncode=3)
    
    summary(gen.data)
    bartlett.test(list(summary(gen.data)$Hexp, summary(gen.data)$Hobs))
    
    gen.hier <- genind2hierfstat(gen.data, pop = gen.data@pop)
    gen.hier$pop <- paste0("P-", gen.hier$pop)
    basicstat <- basic.stats(gen.hier, diploid = TRUE)
    names(basicstat)
    
    boot.ppfis(gen.hier) 
    indpca(gen.hier)
    plot(indpca(gen.hier), cex = 0.7)
    
    
    setPop(gen.data) <- ~pop
    gen_diff <- diff_stats(gen.data)
    
    
    
    
    
      ## Manipulate data? Export Genepop from gstudio?
      gen.data <- inds
      gen.data$ind <- paste0("NP", gen.data$Nat.Pond, "-BP", gen.data$Breed.Pond, "-", gen.data$Sex, "-", row.names(gen.data))
      gen.data$pop <- paste0("Pond-", gen.data$Breed.Pond)
      gen.data <- gen.data[ , c("ind", "pop", "Patch.X", "Patch.Y", "Sex", "Age", 
                                "LocA", "LocB", "LocC", "LocD", "LocE", "LocF", "LocG", "LocH", "LocI", "LocJ", "LocK",
                                "LocL", "LocM", "LocN", "LocO", "LocP", "LocQ", "LocR", "LocS", "LocT", "LocU")]
    
      gen.data2 <- df2genind(X = gen.data, sep=":", ploidy=2, pop=gen.data$pop, ind.names=gen.data$ind)
      g.out <- popgenreport(gen.data2, mk.counts = T, mk.fst = T, mk.allele.dist = T, mk.pdf = F)
      #write.csv()  
## ------------------
    
    
## Plot Demographic Data:
## ----------------------
  par(mfrow=c(2,3))
   
    hist(inds0$Age, col=rgb(0,0,1,0.5), xlim=c(0,max.age+1),
         main="Age Distribution - Initial", xlab="Age (years)")
    plot(inds0$SVL ~ jitter(inds0$Age, 1), col=ifelse(inds0$Sex=="M", "blue", "red"), pch=16, cex=0.75, 
         main="Age x SVL x Sex - Inital", xlab="Age (years)", ylab="SVL (mm)")
    plot(inds0$SVL ~ jitter(inds0$Age, 1), col=ifelse(inds0$Rep.Active==T, "orange", "black"), pch=16, cex=0.75, 
         main="Age x SVL x Rep. Active - Initial", 
         xlab="Age (years)", ylab="SVL (mm)")
    
    hist(inds$Age, col=rgb(1,0,0,0.5), xlim=c(0,max.age+1), 
         main="Age Distribution - Final", xlab="Age (years)")
    plot(inds$SVL ~ jitter(inds$Age, 1), col=ifelse(inds$Sex=="M", "blue", "red"), pch=16, cex=0.75, 
         main="Age x SVL x Sex - Final", xlab="Age (years)", ylab="SVL (mm)")
    plot(inds$SVL ~ jitter(inds$Age, 1), col=ifelse(inds$Rep.Active==T, "orange", "black"), pch=16, cex=0.75,
         main="Age x SVL x Rep. Active - Final", 
         xlab="Age (years)", ylab="SVL (mm)")
    
## ----------------------