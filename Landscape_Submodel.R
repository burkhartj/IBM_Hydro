#######################
# Build initial landscape
#######################

# Make an empty landscape
r <- raster(matrix(0, landscape.x.rows, landscape.y.columns),
            xmn = landscape.x.coor.min,
            xmx = landscape.x.coor.max,
            ymn = landscape.y.coor.min,
            ymx = landscape.y.coor.max)   #Set raster min and max y and x coords

#Make layer with ponds
pond.r <- makeClass(r,                      #Raster name
                    pond.num,               #Number of ponds
                    rnorm(pond.num, pond.size.mean, pond.size.sd),      #Pond size
                    val = pond.hydro.class)       #Pond class

## Make border landscape to contain all dispersers within max dispersal distance of max pond raster extent
br <- raster(matrix(0, (dim(pond.r)[1]+ceiling(max.disp.dist*2/30)), (dim(pond.r)[2]+ceiling(max.disp.dist*2/30))), 
             xmn = extent(pond.r)@xmin - (ceiling(max.disp.dist/30) * 30),
             xmx = extent(pond.r)@xmax + (ceiling(max.disp.dist/30) * 30), 
             ymn = extent(pond.r)@ymin - (ceiling(max.disp.dist/30) * 30), 
             ymx = extent(pond.r)@ymax + (ceiling(max.disp.dist/30) * 30))   #Set raster min and max y and x coords

pond.r <- mosaic(br, pond.r, fun=sum)     ## make pond raster with the terrestrial border

#Make terrestrial carrying capacity layer
terrestrial.k.r <- ((br + 1) - pond.r) ^2

#Empty layer to assign residents to 
terrestrial.resident.r <- br

## Create an empty distance layer
dist.r <- br

#Stack layers into single object
ls <- stack(pond.r, terrestrial.k.r, terrestrial.resident.r, dist.r)  #Form stack of raster layers