####To do: 
#1) Make pond sizes FLW relevant
#2) Should 

#########################
# Build raster landscape
# Tutorial: https://cran.r-project.org/web/packages/landscapeR/vignettes/landscapeR.html
#########################

library(landscapeR)
library(raster)
library(rgdal)
library(actuar)

# Make an empty landscape
m <- matrix(0, 100, 100)                                  #Set landscape size
r <- raster(m, xmn = 0, xmx = 3000, ymn = 0 , ymx = 3000)   #Set raster min and max y and x coords

# Add in pond patches. Patches are allowed to be contiguous. Can add in pond starting position
# and split background into multiple classes (eg. make ponds first, add in forest / fields).
pond.num <- 5                        ## number of ponds to create
pond.size.mean <- 1                  ## mean size of pond (in 30x30 m grid cells)
pond.size.sd <- 0                    ## st devation of pond
pond.class <- 1                      ## number of pond classes (used for hydroperiod later?)
max.disp.dist <- 5000                ## maximum dispersal distance (in meters) for creating a buffer (quick and dirty fix for world wrapping)

pond.r <- makeClass(r,                      #Raster name
                    pond.num,               #Number of ponds
                    rnorm(pond.num, pond.size.mean, pond.size.sd),      #Pond size
                    val = pond.class)       #Pond class

plot(pond.r, col = c("green", "blue"))


## Make border landscape to contain all dispersers within max dispersal distance of max pond raster extent
bm <- matrix(0, (dim(pond.r)[1]+ceiling(max.disp.dist*2/30)), (dim(pond.r)[2]+ceiling(max.disp.dist*2/30)))
br <- raster(bm, 
             xmn = extent(pond.r)@xmin - (ceiling(max.disp.dist/30) * 30),
             xmx = extent(pond.r)@xmax + (ceiling(max.disp.dist/30) * 30), 
             ymn = extent(pond.r)@ymin - (ceiling(max.disp.dist/30) * 30), 
             ymx = extent(pond.r)@ymax + (ceiling(max.disp.dist/30) * 30))   #Set raster min and max y and x coords

pond.r <- mosaic(br, pond.r, fun=sum)     ## make pond raster with the terrestrial border

#Create the terrestrial carrying capacity raster
terrestrial.k <-180    #Assuming 180 salamanders per 30 x 30 m cell

terrestrial.k.r <- br
terrestrial.k.r[,] <- 1
terrestrial.k.r <- terrestrial.k.r - pond.r
terrestrial.k.r <- terrestrial.k.r * terrestrial.k


#Add empty layer to assign residents to 
terrestrial.resident.r <- br             ##JJB updated to the "br" object instead of "r" since the "br" layer accounts for the buffer distance

## Create and empty distance layer
dist.r <- br

ls <- stack(pond.r, terrestrial.k.r, terrestrial.resident.r, dist.r)  #Form stack of raster layers

names(ls) <- c('Pond.Loc', 'Terrestrial.k', 'Terrestrial.Resident', 'Dist.r')


