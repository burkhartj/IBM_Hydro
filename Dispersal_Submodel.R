##################
# Works well for one individual with fixed step length
# 
# To do:
# Integrate with landscape module
# Add in sensing for deciding when to stop
################################################


library(actuar)

run.steps <- NULL   
angle <- 1:360        ## seq(from = 1, to = 360, by =8) 
startx <- 0           ## x-coordiante of displacement location
starty <- 0          ## y-coordinate of displacement location 

## Parameters to Change
disp.x <- startx             ## distance lizard displaced in the x-dimension
disp.y <- starty             ## distance lizard displaced in the y-dimension
#sensing.dist <- 30          ## radius in which salamander can detect ponds / neighbors (30 m)
n.ind <- 10              ## number of lizards to test per trial
n.trials <- 1            ## number of trials to run the lizard for

#Parameters to save from model runs
max.dispersal.dist <- NULL         ## Stores the max distance can disperse for each individual
total.dist <- NULL           ## stores the total distance traveled for each lizard
ind.moves <- NULL         ## stores the total number of moves for each lizard     
die <- 0

## Simulation Code
start_time <- Sys.time()         ## Start timer for running code 
  run.steps <- NULL              ## create empty vector for storing the number of movment steps made during trial
  plot.lizard <- sample(x = 1:n.ind, 1)   ## select a random lizard to plot
  
  repeat{                ## Loop over number of lizards for each trial
    new.x <- 0                   ## create starting x-coordinate of displacement?
    new.y <- disp.y                  ## create starting y-coordinate of displacement?
    x <- NULL                    ## create empty vector for storing x-coordinates
    y <- NULL                    ## create empty vector for storing y-coordinates 
    dist <- NULL                 ## create empty vector for storing distances from home
    die.disp <- 0
    new.ang <- sample(angle, replace = TRUE, size = 1)     ## randomly select an initial angle
    dist.max <- rllogis(n = 1, shape = 1.5, scale = 30)  ## number of movements allowed per individual
    ang <- NULL                  ## create empty vector for storing movement angles? 
    dist.traveled <- 0             ## create object for storing the individual lizard movement information

    repeat{              ## Loop over number of movment steps (run.steps) for each lizard 
      new.ang <- rnorm(n = 1, mean = new.ang, sd = 10)     ## randomly select an angle
      new.move <- 10                                         ## move 30 m in each step
      rad <- new.ang*(pi/180)                                 ## convert angle from degrees to radians
      new.x <- new.x + new.move*sin(rad)                      ## calculate new x-coordinate 
      new.y <- new.y + new.move*cos(rad)                      ## calculate new y-coordinate
      dist.traveled <- sqrt((abs(new.x))^2 + (abs(new.y))^2)      ## calculate distance from new position to home
      
      x <- c(x, new.x)                                        ## create vector of x-coordinates for each movement step 
      y <- c(y, new.y)                                        ## create vector of y-coordinates for each movement step
      dist <- c(dist, dist.traveled)                              ## create vector of distances from home for each movement step
      ang <- c(ang, new.ang)                                  ## create vector of angles for each movement step
     
      #If find available home, move on to next individual
      if (dist.traveled >= dist.max) {
        die.disp <- die.disp + 1
        break}
      #if (dist.home <= home.radius | length(dist) == n.moves+1) {break}     ## exit repeat loop if distance is less than 20 meters or if greater than 1000 movement steps
      
      dist.traveled <- dist.traveled + new.move    ## calculate the total distance traveled for the lizard
    }
    
    total.dist <- c(total.dist, dist.traveled)                 ## store total distance traveled for each lizard
    max.dispersal.dist   <- c(max.dispersal.dist, dist.max)
    run.steps <- c(run.steps, length(dist))                  ## create vector of movement steps made to reach home
    die <- c(die, die.disp)

    if (length(run.steps) == n.ind) {break}              ## exit repeat loop when number of run steps is 20?? AKA more than 20 trials were run??
  }
  

  ## Write output data files: 
  #   out.df <- data.frame(final.distance = final.dist, total.distance=total.dist, num.moves = lizard.moves, 
  #                        successful.home = ifelse(lizard.moves > n.moves, yes=0, no=1))                           ## create an output data file for individual lizard movment data 
  #   write.csv(out.df, paste0(disp.y, "m", home.radius, "mrad", n.trials, "distance_moved_20Klizards.csv"))        ## output total distance moved and final distance from home for each lizard
  #   write.table(no.success, paste0(disp.y, "m", home.radius, "mrad", n.trials, "success_binary_20Klizards.txt"))  ## output table tracking number of lizards that successfully home in each trial
  #    
  # ## Data Summaries: 
  #   ddply(out.df, ~as.factor(successful.home), summarise, mean=mean(total.distance, na.rm=T), sd=sd(total.distance, na.rm=T), 
  #         min=min(total.distance, na.rm=T), max=max(total.distance, na.rm=T))
  #   
  #   ddply(out.df, ~as.factor(successful.home), summarise, mean=mean(final.distance, na.rm=T), sd=sd(final.distance, na.rm=T), 
  #         min=min(final.distance, na.rm=T), max=max(final.distance, na.rm=T))