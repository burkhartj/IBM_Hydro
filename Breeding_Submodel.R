## ---------------------------------------------------------
## Test functions for IBM Hydroperiod modelling
## 
## By: J. Burkhart and B. Ousterhout
## ---------------------------------------------------------


## Build Pond and Terrestrial Landscape:
## -------------------------------------
  source("Landscape_Submodel.R")
## -------------------------------------


# ## Initialize Models:  ---> Moved to Model_Initialization.R code 
# ## ------------------
#   ## Input Parameters to automate changes:
#     n.inds <- 100               ## number of individuals to create across all ponds
#     n.ponds <- 5               ## number of initial ponds to create
#      
#     n.gens <- 100              ## number of generations to iterate over 
#     
#     pond.K.mult <- 2.25             ## Carrying capacity multiplier. Based off a Semlitsch paper 
#     patch.K.mult <- 180             ## Carrying capacity multiplier for each 30x30 m grid cell
#     
#     min.SVL.F <- 48            ## minimum SVL for sexual maturity - females ---- ARBITRARY RIGHT NOW 
#     min.SVL.M <- 45            ## minimum SVL for sexual maturity - males ---- ARBITRARY RIGHT NOW 
#     SVL.0.mu <- 26.00          ## SVL params based on Taylor and Scott equations, MU's are data driven, SD's ARE ARBITRARY (sizes are for AMOP, maybe use Winters thesis [skeletal chronology] for AMAN size classes?), 
#     SVL.0.sd <- 1.5
#     SVL.1.mu <- 39.61
#     SVL.1.sd <- 0.5
#     SVL.2.mu <- 44.47
#     SVL.2.sd <- 0.25
#     SVL.3.mu <- 48.47 
#     SVL.3.sd <- 0.10
#     SVL.4p.mu <- 48.97
#     SVL.4p.sd <- 0.01
# 
#     min.age.F <- 2                   ## minimum age for reproduction - females
#     min.age.M <- 1                   ## minumum age for reproduction - males
#     max.age <- 15                    ## maximum age for all adults (based on winters thesis)
#     
#     max.disp.dist <- 5000            ## maximum possible dispersal distance (used for landscape border buffer)
#     philo.rate <- 0.90               ## rate of philopatry
#     disp.shape <- 1.5                ## shape param for rllogis dispersal dist draws
#     disp.scale <- 30                 ## scale param for rllogis dispersal dist draws
#     disp.beta1 <- 0.45               ## shape1 for rbeta dispersal dist
#     disp.beta2 <- 1.00               ## shape2 for rbeta dispersal dist
#     mort.prob.mu <- 0.69             ## probability of ADULT mortality - MU for rnorm draw
#     mort.prob.sd <- 0.10             ## probability of ADULT mortality - SD for rnorm draw
#     
#     min.hydro <- 1                   ## minimum hydroperiod class; if < this threshold, no reproduction in pond
#     n.hydro.class <- 3               ## number of hydroperiod classes (used in random draws currently, will be used for deriving parameters later rather than hard coding like below)
#     hydro.0.mu <- 42.8
#     hydro.0.sd <- 5
#     hydro.1.mu <- 76.5
#     hydro.1.sd <- 10
#     hydro.2.mu <- 437.5
#     hydro.2.sd <- 50
#     hydro.3.mu <- 3003.5
#     hydro.3.sd <- 100
#     # hydro.4.mu <- 
#     # hydro.4.sd <-   
# 
# ### --------------------------------
  ## Calculate number of patches in landscape
  n.patch <- ls[[1]]@ncols * ls[[1]]@nrows         ## total number of patches    

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
  
  ## Calculate Pond area based on Hydroperiod classes based upon the Semlitsch diversity/disturbance paper
  for(i in 1:dim(ponds)[1]){
    ponds$Pond.Area[i] <- ifelse(ponds$Hydroperiod[i] == 0, yes=round(rnorm(1, hydro.0.mu, hydro.0.sd), 3), 
                                 no=ifelse(ponds$Hydroperiod[i] == 1, yes=round(rnorm(1, hydro.1.mu, hydro.1.sd), 3),
                                           no=ifelse(ponds$Hydroperiod[i] == 2, yes=round(rnorm(1, hydro.2.mu, hydro.2.sd), 3), 
                                                     no=ifelse(ponds$Hydroperiod[i] == 3, yes=round(rnorm(1, hydro.3.mu, hydro.3.sd), 3), 
                                                               no=NA))))
  }
  
  ponds$Pond.K <- round(pond.K.mult * ponds$Pond.Area)                             ## calculate pond K based upon 
  
  ponds$Pond.X <- as.data.frame(rasterToPoints(ls[[1]], function(x){x == 1}))$x     ## extract x-coordinates from pond raster (pond.r)
  ponds$Pond.Y <- as.data.frame(rasterToPoints(ls[[1]], function(x){x == 1}))$y     ## extract y-coordinates from pond raster (pond.r)
          
  
  ## Create Data Frame of Individual Demographic Information:
  inds <- data.frame(Ind.ID = character(n.inds),
                     Sex = sample(c("M", "F"), n.inds, replace=T), 
                     Age = sample(0:round(max.age / 2), n.inds, T),      
                     SVL = numeric(n.inds),   
                     Rep.Active = NA,                                                 ## counter to assess whether reproductively active
                     IBI = numeric(n.inds),                                           ## inter-breeding interval 
                     Bred = 0, 
                     Generation = 0, 
                     Nat.Pond = sample(1:n.ponds, n.inds, T),                         ## Natal pond
                     Patch.X = numeric(n.inds), 
                     Patch.Y = numeric(n.inds),
                     Disp.Prob = round(rbeta(n.inds, disp.beta1, disp.beta2), 3),      
                     Breed.Pond = numeric(n.inds),                                    ## Breeding pond
                     Init.Angle = numeric(n.inds),                                    ## Angle (in degrees) for which individuals leave pond
                     dist.max = numeric(n.inds),                                      ## max dispersal distance
                     Mort.Prob = numeric(n.inds),                                     ## probability of mortality
                     LocA=numeric(n.inds), LocB=numeric(n.inds),                      ## genetic info, add number loci equal to input data
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
  
     inds$Breed.Pond <- inds$Nat.Pond   #ifelse(inds$Disp.Prob > philo.rate, yes=sample(unique(inds$Nat.Pond), 1, T), no=inds$Nat.Pond)
     inds$Ind.ID <- paste0("N", inds$Nat.Pond, "-B", inds$Breed.Pond, "-G", inds$Generation, "-", rownames(inds), "-", inds$Sex)
     
     ponds$N.inds <- as.data.frame(table(inds$Breed.Pond))$Freq        ## calculate the number of individuals assigned to each pond
     
     
    ## Calculate SVL, Reproductive Activity Status, Genetic Data, and Assign to terrestrial patches 
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
      ## Grab two random alleles from each loci per individual per locus
        inds$LocA[i] <- paste0(sample(unlist(strsplit(gen.df$Aa_37, ":")), 1), ":", sample(unlist(strsplit(gen.df$Aa_37, ":")), 1))
        ## -----
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
      ## -----
      
      ## Assign individuals to a terrestrial patch 
      inds$Init.Angle[i] <- round(runif(1, 1, 360))                                                                                     ## draw dispersal angle
      inds$dist.max[i] <- rtrunc(n=1, spec="llogis", a=min.disp.dist, b=max.disp.dist, shape=disp.shape, scale=disp.scale)              ## draw a maximum dispersal distance
      inds$Patch.X[i] <- ponds[ponds$Pond.ID %in% inds$Nat.Pond[i], "Pond.X"] + inds$dist.max[i] * sin((inds$Init.Angle[i])*(pi/180))   ## move individual from pond in x-coor
      inds$Patch.Y[i] <- ponds[ponds$Pond.ID %in% inds$Nat.Pond[i], "Pond.Y"] + inds$dist.max[i] * cos((inds$Init.Angle[i])*(pi/180))   ## move individual from pond in y-coor

        ## Check that number of individuals on patch is < K/2; 
          ## if yes -> assign to terrestrial.resident.r; if no -> loop until find unoccupied patch
        repeat{
          Patch.K <- extract(ls[[2]], cbind(inds$Patch.X[i], inds$Patch.Y[i]))                            ## extract patch K from 'terrestrial.k.r' layer
          patch.inds <- ls[[3]][cellFromXY(ls, cbind(inds$Patch.X[i],inds$Patch.Y[i]))]                   ## pull # of inds from 'terrestrial.resident.r' layer
                        #terrestrial.resident.r[cellFromXY(ls, cbind(inds$Patch.X[i],inds$Patch.Y[i]))]
                            #subset(inds, Patch.X == inds$Patch.X[i] & Patch.Y == inds$Patch.Y[i]))[1]   ## calculate # inds on patch Resident surface and update that. 
          
          print(paste0("Check Terr K for ind #", i, " ----- Patch K = ", Patch.K, " ----- Patch.inds = ", patch.inds))  ## report progress
          
          if(patch.inds > floor(Patch.K / 2) | is.na(Patch.K)==T) {   ## check if patch has less than half carrying capacity, if not -> redraw patch
            inds$Init.Angle[i] <- round(runif(1, 1, 360))             
            inds$dist.max[i] <- rtrunc(n=1, spec="llogis", a=min.disp.dist, b=max.disp.dist, shape=disp.shape, scale=disp.scale)
            inds$Patch.X[i] <- ponds[ponds$Pond.ID %in% inds$Nat.Pond[i], "Pond.X"] + inds$dist.max[i] * sin((inds$Init.Angle[i])*(pi/180))
            inds$Patch.Y[i] <- ponds[ponds$Pond.ID %in% inds$Nat.Pond[i], "Pond.Y"] + inds$dist.max[i] * cos((inds$Init.Angle[i])*(pi/180))
          }
            else { 
              # print(i, " - passed")
              ls[[3]][cellFromXY(ls, cbind(inds$Patch.X[i],inds$Patch.Y[i]))] <- ls[[3]][cellFromXY(ls, cbind(inds$Patch.X[i], inds$Patch.Y[i]))]+1
              # terrestrial.resident.r[cellFromXY(ls, cbind(inds$Patch.X[i],inds$Patch.Y[i]))] <- terrestrial.resident.r[cellFromXY(ls, cbind(inds$Patch.X[i], inds$Patch.Y[i]))]+1
              break
            }      ## end else statement
        }      ## end repeat loop
    }      ## end for loop

    plot(ls[[3]], main="Initial Resident Layer (pre loop)", #terrestrial.resident.r, 
         xlim=c(-690,3690), ylim=c(-690,3690))#, 
         # breaks=c(0,1,10,25,50,100,180), 
         # col=c("black", "grey", "tan", "yellow", "orange", "red", "purple")) 
    points(ponds$Pond.X, ponds$Pond.Y, 
           cex=as.numeric(ponds$Hydroperiod)+1, col="magenta") 
    
    table(as.matrix(ls[[3]]))
## ------------------  

     
## Breeding Migration:
## -------------------
inds0 <- inds
ponds$Generation <- 0                             ## add generation counter
pond.output <- ponds[-c(1:dim(ponds)[1]), ]  
gen.output <- inds

start.time <- Sys.time()       
g <- 0
repeat{   
  g <- g + 1                                           ## set up generation counter
  ## Set Up Breeding Functionality
  p.sub <- subset(ponds, Hydroperiod >= min.hydro)   ## Subset out the pond with non-suitable hydroperiods
  pond.list <- unique(p.sub$Pond.ID)                 ## Create a vector of unique ID's to iterate over.  CAN BE DELETED IN THE FUTURE; however, 
  off.df <- inds[-c(1:dim(inds)[1]), ]               ## create empty data frame to populate with offspring within the loop
  

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
      
      for(f in 1:dim(rep.feme)[1]){                                                     ## loop over reproductively active females
        n.off <- round((8.47 * rep.feme$SVL[f] - 380) * abs(rnorm(1, 0.355, 0.11)))     ## Size based fecundity * survival to age 1
        n.off <- ifelse(n.off < 0, yes=0, no=n.off)                                     ## adjust number of offspring if negative value (happens with the -380, occassionally)
      
        male.df <- rep.male[sample(1:dim(rep.male)[1], n.off, replace=T), ]             ## data frame of males chosen to breed with a female (sampled with replacement) 
        
        if(n.off > 0){                                                                  ## check that there are offspring to create, if so -> breed, if not -> next reproductive female
          num.off <- as.data.frame(matrix(nrow = n.off, ncol=dim(inds)[2]))             ## create empty data frame to populate w/ female's offspring during breeding loop
          colnames(num.off) <- colnames(inds)                                           ## assign column names from "ind" data frame, in the same order (essentially for later joining)
          
          ## Demographics
          num.off$Ind.ID <- character(n.off)                                   ## create empty individual ID vector
          num.off$Sex <- sample(c("F", "M"), n.off, T)                         ## randomly assign sex (50:50 draw)
          num.off$Age <- rep(0, times=n.off)                                   ## initialize age
          num.off$SVL <- round(rnorm(n.off, SVL.0.mu, SVL.0.sd), 2)            ## SVL at metamorphosis; ---- UPDATE LATER: random normal draw for now. will update to density dependent
          num.off$Rep.Active <- rep(F, times=n.off)                            ## Set reproductively active to "FALSE"
          num.off$IBI <- rep(0, times=n.off)                                   ## Set interbreeding interval to 0
          num.off$Nat.Pond <- rep(pond.list[p], times=n.off)                   ## Set natal pond to mother's breeding pond
          num.off$Patch.X <- as.numeric(rep(NA, n.off))                                 ## patch x-coordinate; populate with dispersal submodel later
          num.off$Patch.Y <- as.numeric(rep(NA, n.off))                                    ## patch y-coordinate; populate with dispersal submodel later
          num.off$Disp.Prob <- round(rbeta(n.off, disp.beta1, disp.beta2), 3)  ## create a vector of dispersal probability
          num.off$Breed.Pond <- as.numeric(rep(NA, n.off))      ## assing new breeding pond --- DOUBLE CHECK THAT DISPERSAL SUBMODEL DOESN'T OVERRIDE THIS, if so, make an empty numeric vector
          num.off$Init.Angle <- round(runif(n.off, 1, 360))                    ## choose initial movement angle
          num.off$dist.max <- rtrunc(n=1, spec="llogis", a=min.disp.dist, b=max.disp.dist,
                                     shape=disp.shape, scale=disp.scale)       ## Peterman et al. 2015 dispersal distance for A. annulatum? 1,693m (95% CI 1,645-1,740 m for AMAN)
          num.off$Mort.Prob <- numeric(n.off)                                  ## empty vector for later imposing mortality 
          num.off$Bred <- rep(0, times=n.off)                                  ## vector of zeros for breeding
          num.off$Generation <- rep(g, times=n.off)                            ## vector to store generation of origin 
          
          ## Genetics --
            ## System mimics polyandry. Randomly selects one of two females alleles, equal to number of offspring. 
            ## Randomly selects one of males alleles at each locus, each offspring has single paternity.
            ## -----------------
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
            ## -----------------
        }
        
        ## Update output data frame
        off.pond <- rbind(off.pond, num.off)           ## save female offspring data to a temporary pond data frame
        
      }  ## end FEMALES iterations
      
      ## Impose LARVAL carrying capacity catch.  
      # if(ponds[ponds$Pond.ID %in% pond.list[p], "Pond.K"] < dim(off.pond)[1]){
      #    off.pond <- off.pond[-sample(x = 1:dim(off.pond)[1], replace=F,
      #                                 size = (dim(off.pond)[1]) - ponds[ponds$Pond.ID %in% pond.list[p], "Pond.K"]), ]
      # }
       
      ## Update output data set... might be able to pull straight into the 'inds' data frame... 
      off.df <- rbind(off.df, off.pond)             ## save pond data to a temporary generation offspring data frame 
      # print(paste0("number of offspring = ", dim(off.df)[1]))
    } ## end IF statement
  
  } ## end POND iterations
  
  ##########################
  # Dispersal submodel
  ##########################
  if(dim(off.df)[1] > 0) {
    #Parameters to save from model runs
    run.steps <- 0                   ## create empty vector for storing the number of movment steps

    total.dist <- NULL               ## stores the total distance traveled for each animal
    ind.moves <- NULL                ## stores the total number of moves for each animal
    die <- NULL
    success <- NULL
    new.move <- 30              ## DELETE LATER. Step length in the model init file

    ## Simulation Code
    repeat{                ## Loop over all offspring
      run.steps <- run.steps + 1
      new.x <- ponds[ponds$Pond.ID %in% off.df$Nat.Pond[run.steps], "Pond.X"] #Extract individuals pond x coord
      new.y <- ponds[ponds$Pond.ID %in% off.df$Nat.Pond[run.steps], "Pond.Y"] #Extract individuals pond y coord
      x <- NULL                    ## create empty vector for storing x-coordinates
      y <- NULL                    ## create empty vector for storing y-coordinates
      dist <- NULL                 ## create empty vector for storing distances from home
      die.disp <- 0
      success.disp <- 0
      disp.to.adj <- 0
      # new.ang <- sample(1:360, replace = TRUE, size = 1)     ## randomly select an initial angle
      # dist.max <- rtrunc(n = 1, spec = "llogis", a = 30, b = 2000, shape = 1.5, scale = 30) #Max dispersal distance ----- Changed to call from inds file
      ang <- NULL                  ## create empty vector for storing movement angles?
      dist.traveled <- 0             ## create object for storing the individual animal movement information

      repeat{              ## Loop over number of movment steps (run.steps) for each animal
        new.ang <- rnorm(n = 1, mean = off.df$Init.Angle[run.steps], sd = 10)     ## randomly select an angle
        new.x <- new.x + new.move*sin(new.ang*(pi/180))      ## calculate new x-coordinate using radians
        new.y <- new.y + new.move*cos(new.ang*(pi/180))      ## calculate new y-coordinate using radians
        dist.traveled <- sqrt(((new.x-ponds[ponds$Pond.ID %in% off.df$Nat.Pond[run.steps], "Pond.X"])^2) +  ## calculate distance from new position to pond
                              ((new.y-ponds[ponds$Pond.ID %in% off.df$Nat.Pond[run.steps], "Pond.Y"])^2))

        #If move max distance, die and move on to next individual
        if (dist.traveled >= off.df$dist.max[run.steps]) {
          die.disp <- 1
          break}

        #If find available home (Terrestrial.Resident < Terrestrial.k), stop and move on to next individual
        if (extract(ls[[1]], cbind(new.x,new.y)) == 0 &       #If not on pond cell
            extract(ls[[3]], cbind(new.x,new.y)) < extract(ls[[2]], cbind(new.x,new.y))){ #And k is greater than number of patch occupants
            success.disp <- 1
            # disp.cell <- cellFromXY(ls, cbind(new.x,new.y))
            off.df$Patch.X[run.steps] <- new.x
            off.df$Patch.Y[run.steps] <- new.y
              # terrestrial.resident.r <- raster(ls, layer = 3)
            ls[[3]][cellFromXY(ls, cbind(new.x,new.y))] <- ls[[3]][cellFromXY(ls, cbind(new.x,new.y))]+1    # terrestrial.resident.r[cellFromXY(ls, cbind(new.x,new.y))] <- terrestrial.resident.r[cellFromXY(ls, cbind(new.x,new.y))]+1
            # ls <- stack(pond.r, terrestrial.k.r, terrestrial.resident.r, dist.r)
          break}

        #Check neighborhood for available home (Terrestrial.Resident < Terrestrial.k), stop and move on to next individual
        adj.cells <- adjacent(x = ls[[3]], #terrestrial.resident.r,
                              cells = cellFromXY(ls[[3]], cbind(new.x,new.y)),   #terrestrial.resident.r
                              directions = sensing.matrix, pairs = FALSE, sorted=TRUE)
        available <- extract(x = ls[[2]], y = adj.cells) - extract(x = ls[[3]], y = adj.cells)   # extract(x = terrestrial.k.r, y = adj.cells) - extract(x = terrestrial.resident.r, y = adj.cells)

        #Assign salamander to first available
        i <- 0
        repeat{
          i <- i + 1
          if(i > length(available)){
            break}

          if(available[i] > 0){
            disp.to.adj <- 1
            # disp.cell <- adj.cells[i]
            # terrestrial.resident.r <- raster(ls, layer = 3)
            off.df$Patch.X[run.steps] <- xFromCell(ls[[3]], adj.cells[i])  # xFromCell(terrestrial.resident.r, adj.cells[i])
            off.df$Patch.Y[run.steps] <- yFromCell(ls[[3]], adj.cells[i])  # yFromCell(terrestrial.resident.r, adj.cells[i])
            ls[[3]][adj.cells[i]] <- ls[[3]][adj.cells[i]]+1 #terrestrial.resident.r[adj.cells[i]] <- terrestrial.resident.r[adj.cells[i]]+1
            # ls <- stack(pond.r, terrestrial.k.r, terrestrial.resident.r, dist.r)              ## double check if we need to keep re-making these
            break}
        }

        if(disp.to.adj == 1){
          success.disp <- 1
          break}

        dist.traveled <- dist.traveled + new.move    ## calculate the total distance traveled for the animal
      }
      
      
      ## NEED TO FIGURE OUT SOMETHING HERE...RIGHT NOW, THIS IS BREAKING THE LOOP AND ALL INDIVIDUALS 
      ## ARE GETTING ASSIGNED TO A BREEDING POND, REGARDLESS OF WHETHER THEY FIND A TERRESTRIAL PATCH
      ## TRIED AN 'IF STATEMENT' CURRENTLY, WILL SEE IF THAT FIXES THINGS
      if(success.disp == 1){
        d <- pointDistance(p1 =cbind(ponds$Pond.X, ponds$Pond.Y),
                           p2 = cbind(off.df[run.steps,'Patch.X'], off.df[run.steps, 'Patch.Y']),
                           lonlat = FALSE)
        off.df$Breed.Pond[run.steps] <- ponds$Pond.ID[which.min(d)]
      }
      
      total.dist <- c(total.dist, dist.traveled)                      ## store total distance traveled for each animal
      die <- c(die, die.disp)
      success <- c(success, success.disp)

      if (run.steps == dim(off.df)[1]) {break} ##n.ind) {break}         ## exit repeat loop when number of run steps equals number of offspring produced across all ponds
    }
  
    off.df <- subset(off.df, is.na(off.df$Breed.Pond) == F)                      ## POTENTIALLY DELETE THIS LATER. DEPENDS ON HOW THE IF STATEMENT EXECUTES ABOVE
    off.df$Ind.ID <- paste0("N", off.df$Nat.Pond, "-B", off.df$Breed.Pond, "-G", g, "-", rownames(off.df), "-", off.df$Sex)
    # print(plyr::count(die))
    # print(plyr::count(success))
    # print(plot(terrestrial.resident.r))

  }  ## end off spring catch
   # plot(ls[[3]], main=paste0("Post breeding, pre mortality, gen = ", g), #terrestrial.resident.r, 
   #               xlim=c(-690,3690), ylim=c(-690,3690))

  ## Join new and old data frames
  inds <- rbind(inds, off.df)                 ## join surviving offspring data into the full data set
  
  # print(table(inds$Breed.Pond))               ## report the number of individuals that bred
  
## Grow, Mature, Death:
## --------------------
  ## Update demographic information 
  inds$Age <- inds$Age + 1                                                             ## age everybody by one year
  inds$SVL <- ifelse(inds$Age == 1, yes=0.785 * inds$SVL + 19.9,  
                     no=ifelse(inds$Age == 2, yes=0.937 * inds$SVL + 7.36, 
                               no=ifelse(inds$Age == 3, yes=inds$SVL + 2.5, 
                                         no=ifelse(inds$Age == 4, yes=inds$SVL + 1.5, 
                                                   no=ifelse(inds$Age >=5, yes=inds$SVL + (inds$Age - 4) * 0.5, no=NA)))))   ## grow individuals based on Taylor and Scott equations
  
  inds$IBI <- ifelse(inds$Bred == 1, yes = 0, no = inds$IBI + 1)              ## calculate inter-breeding interval
  inds$Generation <- g
  gen.output <- rbind(gen.output, inds)
  
  inds$Rep.Active <- ifelse((inds$Bred==1 & inds$Sex=="F"),                   ## update reproductive activity status
                            yes=F,
                            no=ifelse((inds$Sex=="F" & inds$Age>=min.age.F & inds$SVL>=min.SVL.F), 
                                      yes = sample(c(T, F), dim(subset(inds, Sex=="F" & Age>=min.age.F & SVL>=min.SVL.F)[1]), replace=T), 
                                      no = ifelse((inds$Sex=="M" & inds$SVL>=min.SVL.M & inds$Age>=min.age.M), 
                                                   yes=T, no=inds$Rep.Active)))

  inds$Bred <- 0                                ## reset breeding counter
  
  
  ## Impose Age based mortality
  inds <- inds[which(inds$Age <= max.age), ]                                                        ## keep only individuals less than the hard upper age cut off
  inds$Mort.Prob <- runif(dim(inds)[1], 0, 1) > rnorm(dim(inds)[1], mort.prob.mu, mort.prob.sd)     ## create T/F vector for imposing mortality
  inds <- inds[which(inds$Mort.Prob == F), ]                                                        ## keep only individuals that pass the random mortality catch
  
  ## Update terrestrial resident layers
  count.cell <- count(cellFromXY(ls[[3]], cbind(inds$Patch.X, inds$Patch.Y)))      ## count the frequency of individuals 
  ls[[3]][] <- 0                                                                   ## reset occupancy to zero 
  ls[[3]][count.cell$x] <- count.cell$freq                                         ## update residency with cell counts from 'inds' dataframe
   
  
  ## Update count of individuals per pond
  # find.ponds <- p.test$Pond.ID %in% as.numeric(as.character(unique(as.data.frame(table(inds$Breed.Pond))$Var1)))
  ponds[ponds$Pond.ID %in% as.numeric(as.character(unique(as.data.frame(table(inds$Breed.Pond))$Var1))), "N.inds"] <- as.data.frame(table(inds$Breed.Pond))$Freq      ## calculate data for occupied ponds
  ponds[!ponds$Pond.ID %in% as.numeric(as.character(unique(as.data.frame(table(inds$Breed.Pond))$Var1))), "N.inds"] <- 0                                              ## set unoccupied pond counts to zero 
  ponds$Generation <- g
  # print(ponds)
  
  pond.output <- rbind(pond.output, ponds)             ## print pond counts
  
  print(paste0("Generation = ", g, " --- Total # Inds = ", dim(inds)[1], 
               "; Terrestrial K = ", terrestrial.k * n.patch))              ## report progress and pop sizes
  print(round(Sys.time() - start.time, 2))                                  ## report time progress
  
  ## Temporary tracker for terrestiral.resident.r surface behavior
  # print(plot(ls[[3]], main=paste0("generation - " , g), xlim=c(-690,3690), ylim=c(-690,3690)))
 
  if (g >= n.gens | dim(inds)[1] == 0) {break}   ## end loop if the number of individuals = 0 or the max number of generations is reached. 
} ## end generations loop 
print(round(Sys.time() - start.time, 2))         ## end timer

# gen.magic <- magic_result()                        ## stores all the individual data in a pseudo array, can access by calling gen.magic[[generation number]]
## --------------------
  

## Save Genetic Data:
# ## ------------------
  ## Manipulate data? Export Genepop from gstudio?
    g.magic <- subset(gen.output, Generation == n.gens)                    ## extract generation number from uber array
    # g.magic$Ind.ID <- paste0(g.magic$Ind.ID, "-", g.magic$Sex)        ## create an individual ID metric
    g.df <- g.magic[ , c("Ind.ID", "Breed.Pond",
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
   
    plot(indpca(gen.hier), col=rainbow(40)[as.numeric(gen.hier$pop)], cex = 0.7)
   
    # setPop(gen.data) <- ~pop
    # gen_diff <- diff_stats(gen.data)





      ## Manipulate data? Export Genepop from gstudio?
      # gen.data <- inds
      # # gen.data$ind <- paste0("NP", gen.data$Nat.Pond, "-BP", gen.data$Breed.Pond, "-", gen.data$Sex, "-", row.names(gen.data))
      # # gen.data$pop <- paste0("Pond-", gen.data$Breed.Pond)
      # # gen.data <- gen.data[ , c("ind", "pop", "Patch.X", "Patch.Y", "Sex", "Age",
      # #                           "LocA", "LocB", "LocC", "LocD", "LocE", "LocF", "LocG", "LocH", "LocI", "LocJ", "LocK",
      # #                           "LocL", "LocM", "LocN", "LocO", "LocP", "LocQ", "LocR", "LocS", "LocT", "LocU")]
      # 
      # gen.data2 <- df2genind(X = gen.data, sep=":", ploidy=2, pop=gen.data$pop, ind.names=gen.data$ind)
      # g.out <- popgenreport(gen.data, mk.counts = T, mk.fst = T, mk.allele.dist = T, mk.pdf = F)
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
  
  ## Plot population sizes across time    
    pond.output$hydro_class <- as.factor(pond.output$Hydroperiod)
    
    p <- ggplot(pond.output, aes(x = Generation, y = N.inds, group = Pond.ID, colour = hydro_class)) +
          geom_line(aes(linetype = hydro_class), size=1, show.legend = T) +
          scale_fill_brewer(palette = "Blues")+ #rainbow(length(unique(pond.output$hydro.fact))))+
          labs(x = "Generation", y = "Number of Individuals") + 
          ggtitle(expression(atop(bold("Population Size by Hydroperiod Class"), 
                                  atop("20 ponds, 4 Hydroperiod Classes, 200 Generations")))) +
          theme_classic()+
          theme(panel.grid.major.y = element_blank(), 
                panel.grid.major.x = element_blank(),
                panel.grid.minor = element_blank(),
                plot.title = element_text(size = rel(1.5), face = "bold", vjust=1.5),
                axis.title.y = element_text(face = "bold"), axis.title.x = element_text(face = "bold"))
    p
    
## ----------------------