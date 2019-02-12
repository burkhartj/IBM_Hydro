##########################
# Dispersal submodel
##########################

#Parameters to save from model runs
run.steps <- NULL              ## create empty vector for storing the number of movment steps
max.dispersal.dist <- NULL         ## Stores the max distance can disperse for each individual
total.dist <- NULL           ## stores the total distance traveled for each animal
ind.moves <- NULL         ## stores the total number of moves for each animal     
die <- NULL
success <- NULL
fell.off <- NULL
new.move <- 30

## Simulation Code
repeat{                ## Loop over all dispersing animals
  new.x <- rasterToPoints(ls[[1]], function(v){v == 1})[,1] #Extract individuals pond x coord
  new.y <- rasterToPoints(ls[[1]], function(v){v == 1})[,2] #Extract individuals pond y coord
  x <- NULL                    ## create empty vector for storing x-coordinates
  y <- NULL                    ## create empty vector for storing y-coordinates 
  dist <- NULL                 ## create empty vector for storing distances from home
  die.disp <- 0
  success.disp <- 0
  fell.off.disp <-0
  disp.to.adj <- 0
  new.ang <- sample(1:360, replace = TRUE, size = 1)     ## randomly select an initial angle
  dist.max <- rtrunc(n = 1, spec = "llogis", a = 30, b = 2000, shape = 1.5, scale = 30) #Max dispersal distance
  ang <- NULL                  ## create empty vector for storing movement angles? 
  dist.traveled <- 0             ## create object for storing the individual animal movement information
  
  repeat{              ## Loop over number of movment steps (run.steps) for each animal 
    new.ang <- rnorm(n = 1, mean = new.ang, sd = 10)     ## randomly select an angle
    new.x <- new.x + new.move*sin(new.ang*(pi/180))      ## calculate new x-coordinate using radians 
    new.y <- new.y + new.move*cos(new.ang*(pi/180))      ## calculate new y-coordinate using radians
    dist.traveled <- sqrt(((new.x-rasterToPoints(ls[[1]], function(v){v == 1})[,1])^2) +  ## calculate distance from new position to pond
                            ((new.y-rasterToPoints(ls[[1]], function(v){v == 1})[,2])^2))
    
    # #If move off of landscape, die and move on to next individual
    # if (new.x > 1000 | new.y > 1000) {
    #   fell.off.disp <- 1
    #   break}
    
    #If move max distance, die and move on to next individual
    if (dist.traveled >= dist.max) {
        die.disp <- 1
        break}
    
    #If find available home (Terrestrial.Resident < Terrestrial.k), stop and move on to next individual
    if (extract( ls[[1]], cbind(new.x,new.y)) == 0 &       #If not on pond cell
        extract( ls[[3]], cbind(new.x,new.y)) < extract( ls[[2]], cbind(new.x,new.y))){ #And k is  greater than number of patch occupants
        success.disp <- 1
        disp.cell <- cellFromXY(ls, cbind(new.x,new.y))
        terrestrial.resident.r <- raster(ls, layer = 3)
        terrestrial.resident.r[disp.cell] <- terrestrial.resident.r[disp.cell]+1
        ls <- stack(pond.r, terrestrial.k.r, terrestrial.resident.r, dist.r)
        break}
    
    #Check neighborhood for available home (Terrestrial.Resident < Terrestrial.k), stop and move on to next individual
    adj.cells <- adjacent(x = terrestrial.resident.r,
                      cells = cellFromXY(terrestrial.resident.r, cbind(new.x,new.y)),
                      directions = sensing.matrix, pairs = FALSE, sorted=TRUE)
    available <- extract(x = terrestrial.k.r, y = adj.cells) - extract(x = terrestrial.resident.r, y = adj.cells)

    #Assign salamander to first available
    i <- 0
    repeat{
      i <- i + 1
      if(i > length(available)){
         break}

      if(available[i] > 0){
         disp.to.adj <- 1
         disp.cell <- adj.cells[i]
         terrestrial.resident.r <- raster(ls, layer = 3)
         terrestrial.resident.r[disp.cell] <- terrestrial.resident.r[disp.cell]+1
         ls <- stack(pond.r, terrestrial.k.r, terrestrial.resident.r, dist.r)
         break}
      }
  
    if(disp.to.adj == 1){
      success.disp <- 1
      break}
    
    dist.traveled <- dist.traveled + new.move    ## calculate the total distance traveled for the animal
  }
  
  total.dist <- c(total.dist, dist.traveled)                 ## store total distance traveled for each animal
  max.dispersal.dist   <- c(max.dispersal.dist, dist.max)
  run.steps <- c(run.steps, length(dist))                  ## create vector of movement steps made to reach home
  die <- c(die, die.disp)
  success <- c(success, success.disp)
  fell.off <- c(fell.off, fell.off.disp)
  
  if (length(run.steps) == n.ind) {break}              ## exit repeat loop when number of run steps is 20?? AKA more than 20 trials were run??
}

count(die)
count(success)
count(fell.off)
count(max.dispersal.dist - total.dist > 0)

plot(terrestrial.resident.r)
