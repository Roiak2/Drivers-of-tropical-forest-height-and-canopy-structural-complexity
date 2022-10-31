###############################################################################
# SCRIPT FOR SAMPLING DESIGN FOR CANOPY STRUCTURE CHAPTER
# 
# Roí Ankori-Karlinsky, April 2021
###############################################################################


############################# THIS SCRIPT #####################################

# 1. Creates 30m-radius (0.28) shapefile from random stratified sampling based 
#    on forest age across Puerto Rico
#
# 2. Checks for spatial autocorrelation of canopy metrics at different distances
#
# 3. Removes overlap within a minimum distance of 250m to avoid autocorrelation
#
# 4. Saves shapefile 

############################## PACKAGES #######################################

#Loading packages 

#Basic data management
library(tidyverse) #data manipulation grammar
library(readxl) #read excel files
library(parallel) # mclapply for multicore processing

#Visualization
library(gridExtra) # allows plotting of multiple ggplots in one
library(grid) #same as above
library(colorRamps) # Colors for graphing
library(RColorBrewer) # Colors for graphing
library(cowplot) #putting plots together
library(ggpubr) #plots together
library(scales) #Helps with visualization
library(kableExtra) #summary tables

#Spatial analyses
library(sp) #provides object oriented classes for spatial data in R 
library(sf) #similar to sp
library(rgdal) # raster analysis
library(raster) #provides classes and methods for raster datasets
library(SDMTools) # Spatial analysis in R
library(spdplyr) #spatial data manipulation
library(spdep) #Moran's I spatial autocorrelation
library(spatialEco) #Spatial autocorrelation


#Making sure R is utilizing all available RAM
memory.limit(size=500000)
rasterOptions(maxmemory = 4e+16)
gc()

########################## WORKING DIRECTORIES ################################

#File path
local_path <- 'C:/Users/roiak/Documents/Damage Project/Rugosity'

#Where local shapefiles are
shp_dir <- paste0(local_path,'/Data/Shapefiles')

#Where island points and buffers will be stored
point_dir <- paste0(shp_dir,"/Island Points")

#Where forest age data is
age_dir <- paste0(shp_dir,"/Forest age")

#Where canopy structure metrics are
csc_dir <- paste0(shp_dir,"/CSC")

#-----------------------------------------------------------------------------#
############################ SAMPLING BY FOREST AGE ###########################

#here I'll upload the forest age raster and perform random stratified sampling
#to sample across the island by age class

#setting working directory to where age data is
setwd(age_dir)

#reading raster
age <- raster("Forest_Age.tif")

#saving coordinate systems
nad<-crs(age)
wgs<-"+proj=longlat +datum=WGS84"

#Transform projection to WGS84 (depreceated for now)
#age <- projectRaster(from=age,crs=wgs)

#checking
str(age)
crs(age)
res(age)

#visualizing
plot(age)

#METADATA FOR FOREST AGES
#0 = not forest
#1 = 5-16 yr
#2 = 17-25 yr
#3 = 26-39 yr
#4 = 40-65 yr
#5 = 66+ yr

#random stratified sampling (7500 points for each age class)
points <- sampleStratified(age,
                           size = 7500,
                           xy = TRUE,
                           sp = TRUE) #return spatial points data frame

#clean space
gc()

#Checking
crs(points) 
str(points)
extent(points)
extent(age)

#Removing non-forest points and those not on mainland Puerto Rico
points <- points %>%
  filter(Forest_Age !=0,
         y > 210000)

#removing unnecessary things
rm(age)
gc()

##################### ADDING CANOPY STRUCTURE METRICS #########################

#Now adding canopy height, top rugosity

####--------------------------- CANOPY HEIGHT -----------------------------####

#Loading canopy height
setwd(csc_dir)
CHM <- raster("CHM_Island.tif")
CHM <- readAll(CHM) #making sure we have values

#checking projections
projection(CHM)

#matching projection of shapefile to raster
Buf30m <- spTransform(points,proj4string(CHM))

#checking it worked
projection(Buf30m)
plot(CHM)
plot(Buf30m,add=T)

#Extracting mean, max, and sd canopy height
Buf30m <- Buf30m %>%
  mutate(
    Mean_Height <- raster::extract(CHM, ., fun = mean, na.rm=T)
  )

#cleaning up space
gc()

#Making columns numeric and better named
Buf30m$Mean_Height <- as.numeric(Buf30m$`Mean_Height <- raster::extract(CHM, ., fun = mean, na.rm=T)`)

#Removing badly named columns
Buf30m <- Buf30m %>%
  dplyr::select(-`Mean_Height <- raster::extract(CHM, ., fun = mean, na.rm=T)`
  )

#clearing space
gc()


####-------------------------- CANOPY RUGOSITY ----------------------------####

#Loading canopy rugosity
setwd(csc_dir)
Rugosity <- raster("Rugosity_Island.tif")
#checking
crs(Rugosity)
plot(Rugosity)

#Extracting mean, max, and sd canopy top rugosity to each dataset

Buf30m <- Buf30m %>%
  mutate(
    Mean_Rugosity <- raster::extract(Rugosity, ., fun = mean, na.rm=T)
  )

#cleaning up space
gc()

#Making columns numeric and better named
Buf30m$Mean_Rugosity <- as.numeric(Buf30m$`Mean_Rugosity <- raster::extract(Rugosity, ., fun = mean, na.rm=T)`)

#Removing badly named columns
Buf30m <- Buf30m %>%
  dplyr::select(-`Mean_Rugosity <- raster::extract(Rugosity, ., fun = mean, na.rm=T)`
  )

#clearing space
gc()


#-----------------------------------------------------------------------------#
######################## CHECKING SPATIAL AUTOCORRELATION #####################

## Spatial Autocorrelation

#Going to run a Moran's I test on all CSC metrics and plot lagged.

#Running this on a random subsample of 3000 points (for computational efficiency).
#Running on 5 increasing distances of radius search for neighbors - 50m, 150m, 250m, 500m, and 1km

#This will give us the optimal minimum distance to have to avoid autocorrelation 
#while maximizing sample size

####----Starting with Rugosity----####


#Running a Moran's I spatial autocorrelation test on rugosity using MC method
#with the spdep package

#loading spatial points data frame
spatial_dat <- Buf30m
rm(Buf30m)
gc()

#saving projections
wgs <- CRS("+proj=longlat +datum=WGS84")
nad <- CRS("+proj=lcc +lat_0=17.8333333333333
+lon_0=-66.4333333333333
+lat_1=18.4333333333333
+lat_2=18.0333333333333 +x_0=200000
+y_0=200000 +datum=NAD83 +units=m
+no_defs")

#subsetting random sample of 3000 points for analysis to start with
set.seed(10)
spatial_dat1 <- spatial_dat[sample(1:length(spatial_dat),3000),]
#convert to wgs
spatial_dat1 <- spTransform(spatial_dat1,wgs)
#removing NAs
spatial_dat1 <- sp.na.omit(spatial_dat1)

#getting coordinates to calculate distances
coo <- coordinates(spatial_dat1)

#defining search radius for neighbors (using 50m now)
S.dist  <-  dnearneigh(coo, 0, 50)

#creating spatial neighbor list (i.e. the weights used in the Moran test)
weight <- nb2listw(S.dist, style='B',zero.policy = T)

#Running Monte-Carlo Moran's I test with 1,000 simulations
MI <- moran.mc(spatial_dat1$Mean_Rugosity,
               weight,
               nsim=999,
               zero.policy=T,
               alternative="less")

#Summary model
MI

#standardizing response variable
x <- as.vector(scale(spatial_dat1$Mean_Rugosity))

#--PLOTTING RESULTS--#

#histogram of residuals, if skewed indicates autocorrelation
hist(MI$res, breaks = 10,
     main="Histogram of Residuals Moran's I",
     xlab="Moran's I test residuals")
abline(v = MI$statistic, col="red", lwd=3, lty=2)

#Mean_Rugosity against rugosity spatially lagged
moran.plot(x,weight,
           main="Mean_Rugosity by Lags Moran's I",
           ylab="Spatially Lagged Mean_Rugosity (standardized)",
           xlab="Canopy Mean_Rugosity (standardized)")


#--NOW DOING FOR LONGER DISTANCE--#
rm(MI,S.dist,weight)
gc()

#defining search radius for neighbors (using 150m now)
S.dist  <-  dnearneigh(coo, 0, 150)

#creating spatial neighbor list (i.e. the weights used in the Moran test)
weight <- nb2listw(S.dist, style='B',zero.policy = T)

#Running Monte-Carlo Moran's I test with 1,000 simulations
MI <- moran.mc(spatial_dat1$Mean_Rugosity,
               weight,
               nsim=999,
               zero.policy=T,
               alternative="less")

#Summary model
MI

#standardizing response variable
x <- as.vector(scale(spatial_dat1$Mean_Rugosity))

#--PLOTTING RESULTS--#

#histogram of residuals, if skewed indicates autocorrelation
hist(MI$res, breaks = 10,
     main="Histogram of Residuals Moran's I",
     xlab="Moran's I test residuals")
abline(v = MI$statistic, col="red", lwd=3, lty=2)

#Mean_Rugosity against rugosity spatially lagged
moran.plot(x,weight,
           main="Mean_Rugosity by Lags Moran's I",
           ylab="Spatially Lagged Mean_Rugosity (standardized)",
           xlab="Canopy Mean_Rugosity (standardized)")


#--NOW DOING FOR LONGER DISTANCE--#
rm(MI,S.dist,weight)
gc()

#defining search radius for neighbors (using 250m now)
S.dist  <-  dnearneigh(coo, 0, 250)

#creating spatial neighbor list (i.e. the weights used in the Moran test)
weight <- nb2listw(S.dist, style='B',zero.policy = T)

#Running Monte-Carlo Moran's I test with 1,000 simulations
MI <- moran.mc(spatial_dat1$Mean_Rugosity,
               weight,
               nsim=999,
               zero.policy=T,
               alternative="less")

#Summary model
MI

#--PLOTTING RESULTS--#

#histogram of residuals, if skewed indicates autocorrelation
hist(MI$res, breaks = 10,
     main="Histogram of Residuals Moran's I",
     xlab="Moran's I test residuals")
abline(v = MI$statistic, col="red", lwd=3, lty=2)

#Mean_Rugosity against rugosity spatially lagged
moran.plot(x,weight,
           main="Mean_Rugosity by Lags Moran's I",
           ylab="Spatially Lagged Mean_Rugosity (standardized)",
           xlab="Canopy Mean_Rugosity (standardized)")


#--NOW DOING FOR LONGER DISTANCE--#
rm(MI,S.dist,weight)
gc()

#defining search radius for neighbors (using 1km now)
S.dist  <-  dnearneigh(coo, 0, 1000)

#creating spatial neighbor list (i.e. the weights used in the Moran test)
weight <- nb2listw(S.dist, style='B',zero.policy = T)

#Running Monte-Carlo Moran's I test with 1,000 simulations
MI <- moran.mc(spatial_dat1$Mean_Rugosity,
               weight,
               nsim=999,
               zero.policy=T,
               alternative="less")

#Summary model
MI

#--PLOTTING RESULTS--#

#histogram of residuals, if skewed indicates autocorrelation
hist(MI$res, breaks = 10,
     main="Histogram of Residuals Moran's I",
     xlab="Moran's I test residuals")
abline(v = MI$statistic, col="red", lwd=3, lty=2)

#Mean_Rugosity against rugosity spatially lagged
moran.plot(x,weight,
           main="Mean_Rugosity by Lags Moran's I",
           ylab="Spatially Lagged Mean_Rugosity (standardized)",
           xlab="Canopy Mean_Rugosity (standardized)")


#removing
rm(MI,S.dist,weight,x)
gc()



####----Canopy Height----####


#defining search radius for neighbors (using 50m to start with)
S.dist  <-  dnearneigh(coo, 0, 50)

#creating spatial neighbor list (i.e. the weights used in the Moran test)
weight <- nb2listw(S.dist, style='B',zero.policy = T)

#Running Monte-Carlo Moran's I test with 1,000 simulations
MI <- moran.mc(spatial_dat1$Mean_Height,
               weight,
               nsim=999,
               zero.policy=T,
               alternative="less")

#Summary model
MI


#--PLOTTING RESULTS--#

#histogram of residuals, if skewed indicates autocorrelation
hist(MI$res, breaks = 50,
     main="Histogram of Residuals Moran's I",
     xlab="Moran's I test residuals")
abline(v = MI$statistic, col="red", lwd=3, lty=2)

#standardizing response variable
x <- as.vector(scale(spatial_dat1$Mean_Height))

#Height against Height spatially lagged
moran.plot(x,weight,
           main="Canopy Height by Lags Moran's I",
           ylab="Spatially Lagged Height (standardized)",
           xlab="Canopy Height (standardized)")

#--NOW DOING FOR LONGER DISTANCE--#
rm(MI,S.dist,weight)


#defining search radius for neighbors (using 150m now)
S.dist  <-  dnearneigh(coo, 0, 150)

#creating spatial neighbor list (i.e. the weights used in the Moran test)
weight <- nb2listw(S.dist, style='B',zero.policy = T)

#Running Monte-Carlo Moran's I test with 1,000 simulations
MI <- moran.mc(spatial_dat1$Mean_Height,
               weight,
               nsim=999,
               zero.policy=T)

#Summary model
MI

#--PLOTTING RESULTS--#

#histogram of residuals, if skewed indicates autocorrelation
hist(MI$res, breaks = 10,
     main="Histogram of Residuals Moran's I",
     xlab="Moran's I test residuals")
abline(v = MI$statistic, col="red", lwd=3, lty=2)

#Height against Height spatially lagged
moran.plot(x,weight,
           main="Height by Lags Moran's I",
           ylab="Spatially Lagged Height (standardized)",
           xlab="Canopy Height (standardized)")


#--NOW DOING FOR LONGER DISTANCE--#
rm(MI,S.dist,weight)


#defining search radius for neighbors (using 250m now)
S.dist  <-  dnearneigh(coo, 0, 250)

#creating spatial neighbor list (i.e. the weights used in the Moran test)
weight <- nb2listw(S.dist, style='B',zero.policy = T)

#Running Monte-Carlo Moran's I test with 1,000 simulations
MI <- moran.mc(spatial_dat1$Mean_Height,
               weight,
               nsim=999,
               zero.policy=T,
               alternative = "less")

#Summary model
MI

#--PLOTTING RESULTS--#

#histogram of residuals, if skewed indicates autocorrelation
hist(MI$res, breaks = 10,
     main="Histogram of Residuals Moran's I",
     xlab="Moran's I test residuals")
abline(v = MI$statistic, col="red", lwd=3, lty=2)

#Height against Height spatially lagged
moran.plot(x,weight,
           main="Height by Lags Moran's I",
           ylab="Spatially Lagged Height (standardized)",
           xlab="Canopy Height (standardized)")


#--NOW DOING FOR LONGER DISTANCE--#
rm(MI,S.dist,weight)
gc()

#defining search radius for neighbors (using 1km now)
S.dist  <-  dnearneigh(coo, 0, 1000)

#creating spatial neighbor list (i.e. the weights used in the Moran test)
weight <- nb2listw(S.dist, style='B',zero.policy = T)

#Running Monte-Carlo Moran's I test with 1,000 simulations
MI <- moran.mc(spatial_dat1$Mean_Height,
               weight,
               nsim=999,
               zero.policy=T)

#Summary model
MI

#--PLOTTING RESULTS--#

#histogram of residuals, if skewed indicates autocorrelation
hist(MI$res, breaks = 10,
     main="Histogram of Residuals Moran's I",
     xlab="Moran's I test residuals")
abline(v = MI$statistic, col="red", lwd=3, lty=2)

#Height against Height spatially lagged
moran.plot(x,weight,
           main="Height by Lags Moran's I",
           ylab="Spatially Lagged Height (standardized)",
           xlab="Canopy Height (standardized)")


#removing and clearing space
rm(list=setdiff(ls(), c("shp_dir","local_path","point_dir","points","nad","wgs")))
gc()

## 250 M OPTIMAL DISTANCE! ##

#-----------------------------------------------------------------------------#
######################## REMOVING NEIGHBORING POINTS ##########################

#Removing points that are less than 250m apart (to anticipate buffers).
#   This means there'll always be at least 50m between buffers

#creating boolean distance matrix of each point combo (TRUE if within 250m)
points_matrix <- gWithinDistance(points, dist = 250, byid = TRUE) 

#Put all true in the upper triangle of matrix and assign as NA
points_matrix[lower.tri(points_matrix, diag=TRUE)] <- NA

#column sums to get whether a point crosses another within 250m 
v <- colSums(points_matrix, na.rm=TRUE) == 0

#subset shapefile by this criterion, keeping only points at least 250m away
points2<- points[v, ]

#Checking how many points were removed
mosaic::tally(~points2@data$Forest_Age) #at least 3000 points per age, nice!
nrow(points2@data) #20,660 points total great!

#-----------------------------------------------------------------------------#
############################ CREATING BUFFERS #################################

#Here I'll create 30m buffers around points (i.e. 30m radius or 0.28ha sites) 

#30m buffer
Buf30 <- buffer(points2, width=30,dissolve=F)

#checking
crs(Buf30)

#making sure still got same distribution and sample size of forest ages 
mosaic::tally(~Buf30@data$Forest_Age) #still the same
nrow(Buf30@data) #still the same


#-----------------------------------------------------------------------------#
################################ SAVING DATAFRAMES ############################

#saving shapefile
writeOGR(Buf30,
         dsn=file.path(point_dir),
         layer = "Buf30m",
         driver = "ESRI Shapefile")

#writing csv of the data from the shapefile
write_csv(Buf30@data,
          paste0(local_path,
                 "/Data/StratifiedSampling.csv"))

#cleaning
rm(list=ls())
gc()
