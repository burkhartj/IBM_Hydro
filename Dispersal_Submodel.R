#setwd("C:/Users/semlitschlab/Box Sync/Mizzou/Research/Manuel_Lizards/")           ## set working directory to where files are located and where you would like data output to
library(actuar)
#moves <- read.table("lizard_moves.txt", h = F)    ## data file holding the movement data for Anolis lizards
run.steps <- NULL   
angle <- 1:360        ## seq(from = 1, to = 360, by =8) 
homex <- 0            ## x-coordinate for home location
homey <- 0            ## y-coordinate for home location
startx <- 0           ## x-coordiante of displacement location
starty <- 40          ## y-coordinate of displacement location 
no.success <- NULL    ## empty vector to track the number of successes 

## Parameters to Change
disp.x <- startx             ## distance lizard displaced in the x-dimension
disp.y <- starty             ## distance lizard displaced in the y-dimension
home.radius <- 100            ## radius in which salamander can detect pond (100 m)
n.moves <- rllogis(1500              ## number of movements allowed per individual
n.lizards <- 20              ## number of lizards to test per trial
n.trials <- 10000            ## number of trials to run the lizard for
final.dist <- NULL           ## stores the final distance for each individual lizard
total.dist <- NULL           ## stores the total distance traveled for each lizard
lizard.moves <- NULL         ## stores the total number of moves for each lizard     

## Simulation Code
start_time <- Sys.time()         ## Start timer for running code 
repeat{                          ## Loop over number of trials
  success <- NULL                ## create empty vector for success
  run.steps <- NULL              ## create empty vector for storing the number of movment steps made during trial
  # plot.lizard <- sample(x = 1:n.lizards, 1)   ## select a random lizard to plot
  
  repeat{                ## Loop over number of lizards for each trial
    new.x <- 0                   ## create starting x-coordinate of displacement?
    new.y <- disp.y                  ## create starting y-coordinate of displacement?
    x <- NULL                    ## create empty vector for storing x-coordinates
    y <- NULL                    ## create empty vector for storing y-coordinates 
    dist <- NULL                 ## create empty vector for storing distances from home
    ang <- NULL                  ## create empty vector for storing movement angles? 
    dist.travel <- 0             ## create object for storing the individual lizard movement information
    fin.dist <- NULL             ## create object for storing the final distance per individual (fixes bug in code)

    repeat{              ## Loop over number of movment steps (run.steps) for each lizard 
      new.ang <- sample(angle, size = 1)                      ## randomly select an angle
      new.move <- sample(moves$V1, size = 1)                  ## randomly select a movement distance
      rad <- new.ang*(pi/180)                                 ## convert angle from degrees to radians
      new.x <- new.x + new.move*sin(rad)                      ## calculate new x-coordinate 
      new.y <- new.y + new.move*cos(rad)                      ## calculate new y-coordinate
      dist.home <- sqrt((abs(new.x))^2 + (abs(new.y))^2)      ## calculate distance from new position to home
      
      x <- c(x, new.x)                                        ## create vector of x-coordinates for each movement step 
      y <- c(y, new.y)                                        ## create vector of y-coordinates for each movement step
      dist <- c(dist, dist.home)                              ## create vector of distances from home for each movement step
      ang <- c(ang, new.ang)                                  ## create vector of angles for each movement step
     
      #If reach home radius, move on to next individual
      if (dist.home <= home.radius | length(dist) == n.moves+1) {break}     ## exit repeat loop if distance is less than 20 meters or if greater than 1000 movement steps
      
      dist.travel <- dist.travel + new.move    ## calculate the total distance traveled for the lizard
      fin.dist <- dist.home                    ## if lizard has not homed or timed out, update final distance          
    }
    
    dum.dist <- ifelse(length(dist) == n.moves+1, yes=fin.dist, no=dist.home)
    final.dist <- c(final.dist, dum.dist)                    ## store final distance to home for each lizard
    total.dist <- c(total.dist, dist.travel)                 ## store total distance traveled for each lizard
    run.steps <- c(run.steps, length(dist))                  ## create vector of movement steps made to reach home
    
    ## plot lizard locations
      # if(length(dist) < n.moves){
      #   tiff(paste0("Homing-Success_Plot_Trial_", length(no.success)+1,"_Lizard_", length(run.steps), ".tiff"), res=250, width=12.35, height=12.35, units="cm", compression=c("lzw"))
      #   plot(x=x, y=y, type="l", ylim=c(0, (max(y)+5)),
      #        main=paste0("Trial #", length(no.success)+1," - Lizard #", length(run.steps)))  ## plot the x-y locations for lizard as a line
      #   points(x=homex, y=homey, col='blue', cex=1.5, pch=16)                                ## plot "home" location for lizard (0,0)
      #   points(x=homex, y=homey, cex=1.5)
      #   points(x=startx, y=starty, col='yellow', cex=1.5, pch=16)                            ## plot displacement location
      #   points(x=startx, y=starty, cex=1.5)
      #   points(x=x[length(x)], y=y[length(y)], col='red', cex=1.5, pch=16)                   ## plot final lizard location
      #   points(x=x[length(x)], y=y[length(y)], cex=1.5)
      #   dev.off()
      # }
    
    if (length(run.steps) == n.lizards) {break}              ## exit repeat loop when number of run steps is 20?? AKA more than 20 trials were run??
  }
  
  success <- length(run.steps[run.steps <= n.moves])         ## calculate the number of successful returns within 1000 movement steps
  no.success <- c(no.success, success)                       ## create vector for number of successful returns
  print(paste0('completed trail - ', length(no.success)))    ## track progress for each lizard at the 
  
  lizard.moves <- c(lizard.moves, run.steps)                 ## create vector of number of moves each lizard took to go  steps from last
  if (length(no.success) >= n.trials) {break}                ## exit repeat loop after trails
}
Sys.time() - start_time       ## Calculate the length of time that code was run


## Data Summaries:
  ## Calculate the number of trials where animals successfully homed
    table(no.success)         ## Top row is the number of lizards out of 20 that made it home, second row is the number of trials that had that observation (e.g., 0/20 lizards made it home in 9,000 trials)sleep_for_a_minute()
  
  ## Write output data files: 
    out.df <- data.frame(final.distance = final.dist, total.distance=total.dist, num.moves = lizard.moves, 
                         successful.home = ifelse(lizard.moves > n.moves, yes=0, no=1))                           ## create an output data file for individual lizard movment data 
    write.csv(out.df, paste0(disp.y, "m", home.radius, "mrad", n.trials, "distance_moved_20Klizards.csv"))        ## output total distance moved and final distance from home for each lizard
    write.table(no.success, paste0(disp.y, "m", home.radius, "mrad", n.trials, "success_binary_20Klizards.txt"))  ## output table tracking number of lizards that successfully home in each trial
     
  ## Data Summaries: 
    ddply(out.df, ~as.factor(successful.home), summarise, mean=mean(total.distance, na.rm=T), sd=sd(total.distance, na.rm=T), 
          min=min(total.distance, na.rm=T), max=max(total.distance, na.rm=T))
    
    ddply(out.df, ~as.factor(successful.home), summarise, mean=mean(final.distance, na.rm=T), sd=sd(final.distance, na.rm=T), 
          min=min(final.distance, na.rm=T), max=max(final.distance, na.rm=T))
    

## Histogram Code -- IGNORE FOR NOW
  # tmp <- read.table(paste0(disp.y, "m", home.radius, "mrad", n.trials, "experiments.txt"))     ## read in table of successes???
  # 
  # ## plot histogram of the number of successess
  #   par(mai = c(0.9, 1.0, 0.2, 0.2))
  #   hist(tmp$x, xlim = c(0, 20), ylim = c(0, 50), cex.axis = 1.5, main = "", cex.lab = 1.75, font.lab = 2, xlab = "Number Returned")
  #   lines(rep(16, 51), 0:50, lwd = 5, col = "red")
  # 

  ###Raster package has extract function so can look at slice