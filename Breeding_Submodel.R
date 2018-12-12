## ---------------------------------------------------------
## Test functions for IBM Hydroperiod modelling
## 
## By: J. Burkhart and B. Ousterhout
## ---------------------------------------------------------

## Load Packages:
## --------------
  if(!require(actuar)) install.packages('actuar'); library('actuar')
## --------------

## Initialize Models:
## ------------------
  ## Input Parameters to automate changes:
    n.inds <- 50               ## number of individuals to create across all ponds
    n.ponds <- 4               ## number of initial ponds to create
    n.patch <- 144             ## total number of patches ---- TEMPORARY, DELETE WHEN LANDSCAPE UPDATE WORKS
    n.gens <- 200              ## number of generations to iterate over 
    
    K.mult <- 2.25             ## Carrying capacity multiplier. Based off a Semlitsch paper 
    
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
    max.age <- 12                    ## maximum age for all adults
    
    max.disp.dist <- 5000            ## maximum possible dispersal distance (used for landscape border buffer)
    philo.rate <- 0.90               ## rate of philopatry
    disp.shape <- 1.5                ## shape param for rllogis dispersal dist draws
    disp.scale <- 30                 ## scale param for rllogis dispersal dist draws
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
    
    temp.terrestrial.K <- n.patch * 180            ## DELETE LATER. WILL BE IRRELEVANT ONCE THE SPATIAL STUFF IS INCORPORATED
  
  ## Create Data Frame of Patches (eventually pull pond coordinates from the landscape submodel)
    # Make an empty landscape
    m <- matrix(0, sqrt(n.patch), sqrt(n.patch))                #Set landscape size
    r <- raster(m, xmn = 0, xmx = sqrt(n.patch) * 30, ymn = 0 , ymx = sqrt(n.patch) * 30)   #Set raster min and max y and x coords
    
    ## Make border landscape
    bm <- matrix(0, (dim(m)[1]+ceiling(max.disp.dist/30)+1), (dim(m)[2]+ceiling(max.disp.dist/30)+1))
    br <- raster(bm, 
                 xmn = extent(pond.r)@xmin - (ceiling(max.disp.dist/2) + 30),
                 xmx = extent(pond.r)@xmax + (ceiling(max.disp.dist/2) + 30), 
                 ymn = extent(pond.r)@ymin - (ceiling(max.disp.dist/2) + 30), 
                 ymx = extent(pond.r)@ymax + (ceiling(max.disp.dist/2) + 30))   #Set raster min and max y and x coords
    res(br) <- 30
    
    p2 <- pond.r + br
    # Add in pond patches. Patches are allowed to be contiguous. Can add in pond starting position
    # and split background into multiple classes (eg. make ponds first, add in forest / fields).
    pond.num <- n.ponds
    pond.size.mean <- 1
    pond.size.sd <- 0
    pond.class <- 1
    pond.r <- makeClass(r,                      #Raster name
                        pond.num,               #Number of ponds
                        rnorm(pond.num, pond.size.mean, pond.size.sd),      #Pond size
                        val = pond.class)       #Pond class
    
    pond.r <- pond.r + br
    plot(pond.r, col = c("green", "blue"))
    
    #Add layer to raster with terrestrial carrying capacity
    terrestrial.k <-180    #Assuming 180 salamanders per 30 x 30 m cell
    
    terrestrial.k.r <- r + 1
    terrestrial.k.r <- terrestrial.k.r - pond.r
    terrestrial.k.r <- terrestrial.k.r * terrestrial.k
    
  plot(pond.r, col=c("green", "blue"))  
  plot(terrestrial.k.r, col=c("light blue", "dark green"))  

  
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
  
  ponds$Pond.K <- round(K.mult * ponds$Pond.Area)
  
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
                     LocA=numeric(n.inds), 
                     LocB=numeric(n.inds), 
                     LocC=numeric(n.inds), 
                     LocD=numeric(n.inds)
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
      inds$LocA[i] <- paste0(sample(seq(100, 130, by=2), 1), ":", sample(seq(100, 130, by=2), 1))
      inds$LocB[i] <- paste0(sample(seq(151, 171, by=2), 1), ":", sample(seq(151, 171, by=2), 1))
      inds$LocC[i] <- paste0(sample(seq(200, 230, by=2), 1), ":", sample(seq(200, 230, by=2), 1))
      inds$LocD[i] <- paste0(sample(seq(241, 261, by=2), 1), ":", sample(seq(241, 261, by=2), 1))
      
      ## Intialize Dispersal Data: 
      inds$Init.Angle[i] <- round(runif(1, 1, 360)) 
      inds$Disp.Dist[i] <- rllogis(1, shape=disp.shape, scale=disp.scale)
      inds$Patch.X[i] <- ponds[ponds$Pond.ID %in% inds$Nat.Pond[i], "Pond.X"] + inds$Disp.Dist[i] * sin((inds$Init.Angle[i])*(pi/180))
      inds$Patch.Y[i] <- ponds[ponds$Pond.ID %in% inds$Nat.Pond[i], "Pond.Y"] + inds$Disp.Dist[i] * cos((inds$Init.Angle[i])*(pi/180))

      
        # repeat{
          ## extract patch K, patch coords
          ## patch.inds <- dim(subset(inds, Patch.X == inds$Patch.X[i] & Patch.Y[i]))[1]
          ## calculate # inds on patch
          # if(patch.inds <= patch.K / 2) {break}
            # inds$Init.Angle[i] <- round(runif(1, 1, 360)) 
            # inds$Disp.Dist[i] <- rllogis(1, shape=disp.shape, scale=disp.scale)
            # patch.K <- 
            
        # }
        ## Pick intial angle and dispersal distance. Check the patch K. 
        ## If less than K/2, assign individual here. If not, re-draw angle and distance.  
    }
   
  

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
          num.off$Patch.X <- numeric(n.inds)                                   ## patch x-coordinate 
          num.off$Patch.Y <- numeric(n.inds)                                   ## patch y-coordinate
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
            num.off$LocA <- paste0(unlist(strsplit(rep.feme$LocA[f], ":"))[sample(1:2, n.off, T)], ":",  
                                   unlist(strsplit(male.df$LocA, ":"))[sample(1:2, n.off, T) + seq(0, n.off*2-1, by=2)])
            num.off$LocB <- paste0(unlist(strsplit(rep.feme$LocB[f], ":"))[sample(1:2, n.off, T)], ":",  
                                   unlist(strsplit(male.df$LocB, ":"))[sample(1:2, n.off, T) + seq(0, n.off*2-1, by=2)])
            num.off$LocC <- paste0(unlist(strsplit(rep.feme$LocC[f], ":"))[sample(1:2, n.off, T)], ":",  
                                   unlist(strsplit(male.df$LocC, ":"))[sample(1:2, n.off, T) + seq(0, n.off*2-1, by=2)])
            num.off$LocD <- paste0(unlist(strsplit(rep.feme$LocD[f], ":"))[sample(1:2, n.off, T)], ":",  
                                   unlist(strsplit(male.df$LocD, ":"))[sample(1:2, n.off, T) + seq(0, n.off*2-1, by=2)])
          
         ## Dispersal: 
            
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
  # write.csv()  
## ------------------
    
    
## Plot Demographic Data:
## ----------------------
  par(mfrow=c(2,3))
   
    hist(inds0$Age, col=rgb(0,0,1,0.5), xlim=c(0,max.age+1),
         main="Age Dist. - Final", xlab="Age (years)")
    plot(inds0$Age, inds0$SVL, col=as.factor(inds0$Sex), pch=17, 
         main="Age x SVL x Sex - Inital", xlab="Age (years)", ylab="SVL (mm)")
    plot(inds0$Age, inds0$SVL, col=as.factor(inds0$Rep.Active), pch=17,
         main="Age x SVL x Rep. Active - Initial", 
         xlab="Age (years)", ylab="SVL (mm)")
    
    hist(inds$Age, col=rgb(1,0,0,0.5), xlim=c(0,max.age+1), 
         main="Age Dist. - Final", xlab="Age (years)")
    plot(inds$Age, inds$SVL, col=as.factor(inds$Sex), pch=16, 
         main="Age x SVL x Sex - Final", xlab="Age (years)", ylab="SVL (mm)")
    plot(inds$Age, inds$SVL, col=as.factor(inds$Rep.Active), pch=16,
         main="Age x SVL x Rep. Active - Final", 
         xlab="Age (years)", ylab="SVL (mm)")
    
## ----------------------