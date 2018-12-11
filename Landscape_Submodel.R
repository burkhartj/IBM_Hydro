####To do: 
#1) Make pond sizes FLW relevant

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
pond.num <- 5
pond.size.mean <- 16
pond.size.sd <- 5
pond.class <- 1
pond.r <- makeClass(r,                      #Raster name
                    pond.num,               #Number of ponds
                    rnorm(pond.num, pond.size.mean, pond.size.sd),      #Pond size
                    val = pond.class)       #Pond class

plot(pond.r, col = c("green", "blue"))

#Add layer to raster with terrestrial carrying capacity
terrestrial.k <-180    #Assuming 180 salamanders per 30 x 30 m cell

terrestrial.k.r <- r + 1
terrestrial.k.r <- terrestrial.k.r - pond.r
terrestrial.k.r <- terrestrial.k.r * terrestrial.k


#Add empty layer to assign residents to 
terrestrial.resident.r <- r

## Create and empty distance layer
dist.r

ls <- stack(pond.r, terrestrial.k.r, terrestrial.resident.r, dist.r)  #Form stack of raster layers


