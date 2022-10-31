###############################################################################

# SCRIPT FOR GETTING TOGETHER ISLAND-WIDE CANOPY STRUCTURE DATASET

# Roí Ankori-Karlinsky, April 2021

###############################################################################


############################# THIS SCRIPT #####################################

# 1. Loads 20,660 30m-radius shapefiles from stratified sampling
# 2. It then adds LiDAR-derived canopy height and top rugosity to each site
# 3. It then adds annual precipitation from PRISM, 
#     chronic wind exposure as calculated by the GWA 3.0 and a 30m-DEM,
#     and exposure to Hurricanes Hugo and Georges
# 4. It then adds elevation, slope, soil type, and soil water availability
# 5. It then saves all of the resultant data to be used for future analyses

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

#Making sure R is utilizing all available RAM
memory.limit(size=500000)
rasterOptions(maxmemory = 4e+16)
gc()

########################## WORKING DIRECTORIES ################################

#File path
local_path <- 'C:/Users/roiak/Documents/Damage Project/Rugosity'

#Where local shapefiles are
shp_dir <- paste0(local_path,'/Data/Shapefiles')

#Where island points are
point_dir <- paste0(shp_dir,"/Island Points/")

#Where precipitation data is
prec_dir <- paste0(shp_dir,"/climate_data")

#Where forest age data is
age_dir <- paste0(shp_dir,"/Forest age")

#Where soil data is
soil_dir <- paste0(shp_dir,"/Soils")

#Where DEM is
dem_dir <- paste0(shp_dir,"/DEM")

#Where canopy structure metrics are
csc_dir <- paste0(shp_dir,"/CSC")

#Where wind exposure is
wind_dir <- paste0(shp_dir,"/Wind Exposure")

#Where past hurricane exposure is
hur_dir <- paste0(wind_dir,"/Hurricanes")


#-----------------------------------------------------------------------------#

############################ LOADING SHAPEFILES ###############################

# Loading shapefile of 20,660 0.28ha forests with 250m minimum distance between centroids

#Load 30m buffer
Buf30m <- readOGR(dsn=paste0(point_dir,"Buf30m.shp"))

#checking
projection(Buf30m)
#saving projections
wgs <- CRS("+proj=longlat +datum=WGS84")
nad <- CRS(projection(Buf30m))


##################### ADDING CANOPY STRUCTURE METRICS #########################

#Now adding canopy height, top rugosity, and aboveground biomass

####--------------------------- CANOPY HEIGHT -----------------------------####

#Loading canopy height
setwd(csc_dir)
CHM <- raster("CHM_Island.tif")
CHM <- readAll(CHM) #making sure we have values

#checking projections
projection(CHM)

#matching projection of shapefile to raster
Buf30m <- spTransform(Buf30m,proj4string(CHM))

#checking it worked
projection(Buf30m)
plot(CHM)
plot(Buf30m,add=T)

#Extracting mean, max, and sd canopy height
Buf30m <- Buf30m %>%
  mutate(
    Mean_Height <- raster::extract(CHM, ., fun = mean, na.rm=T),
    Max_Height <- raster::extract(CHM, ., fun = max, na.rm=T),
    Height_SD <- raster::extract(CHM, ., fun = sd, na.rm=T)
  )

#cleaning up space
gc()

#Making columns numeric and better named
Buf30m$Mean_Height <- as.numeric(Buf30m$`Mean_Height <- raster::extract(CHM, ., fun = mean, na.rm=T)`)
Buf30m$Height_SD <- as.numeric(Buf30m$`Height_SD <- raster::extract(CHM, ., fun = sd, na.rm=T)`)
Buf30m$Max_Height <- as.numeric(Buf30m$`Max_Height <- raster::extract(CHM, ., fun = max, na.rm=T)`)

#Removing badly named columns
Buf30m <- Buf30m %>%
  dplyr::select(-`Mean_Height <- raster::extract(CHM, ., fun = mean, na.rm=T)`,
                -`Height_SD <- raster::extract(CHM, ., fun = sd, na.rm=T)`,
                -`Max_Height <- raster::extract(CHM, ., fun = max, na.rm=T)`
                )

#clearing space
gc()


####----------------------- ABOVEGROUND BIOMASS ---------------------------####

#Running regression on CHM based on plots aboveground biomass (AGB)
#This is based on Jaz's regression from plots comparing logAGB with canopy height
#The R2 is 0.6

#the regression equation is 2.58179 + 0.15241*CHM which gives logAGB, so we exponentiate to get biomass
Biomass <- calc(chm, 
                fun = function(x) (2.58179 + 0.15241*x))

#clearing space
gc()

#matching projection
Buf30m <- spTransform(Buf30m,proj4string(Biomass))
#checking
projection(Buf30m)

#Extracting mean and sd biomass for each sight
Buf30m <- Buf30m %>%
  mutate(
    Mean_Biomass <- raster::extract(Biomass, ., fun = mean, na.rm=T),
    Biomass_SD <- raster::extract(Biomass, ., fun = sd, na.rm=T)
  )

#clearing space
gc()

#Making it numeric
Buf30m$Mean_Biomass <- as.numeric(Buf30m$`Mean_Biomass <- raster::extract(Biomass, ., fun = mean, na.rm = T)`)
Buf30m$Biomass_SD <- as.numeric(Buf30m$`Biomass_SD <- raster::extract(Biomass, ., fun = sd, na.rm = T)`)

#Removing badly named columns
Buf30m <- Buf30m %>%
  dplyr::select(-`Mean_Biomass <- raster::extract(Biomass, ., fun = mean, na.rm = T)`,
                -`Biomass_SD <- raster::extract(Biomass, ., fun = sd, na.rm = T)`)

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
    Mean_Rugosity <- raster::extract(Rugosity, ., fun = mean, na.rm=T),
    Max_Rugosity <- raster::extract(Rugosity, ., fun = max, na.rm=T),
    Rugosity_SD <- raster::extract(Rugosity, ., fun = sd, na.rm=T)
  )

#cleaning up space
gc()

#Making columns numeric and better named
Buf30m$Mean_Rugosity <- as.numeric(Buf30m$`Mean_Rugosity <- raster::extract(Rugosity, ., fun = mean, na.rm=T)`)
Buf30m$Rugosity_SD <- as.numeric(Buf30m$`Rugosity_SD <- raster::extract(Rugosity, ., fun = sd, na.rm=T)`)
Buf30m$Max_Rugosity <- as.numeric(Buf30m$`Max_Rugosity <- raster::extract(Rugosity, ., fun = max, na.rm=T)`)

#Removing badly named columns
Buf30m <- Buf30m %>%
  dplyr::select(-`Mean_Rugosity <- raster::extract(Rugosity, ., fun = mean, na.rm=T)`,
                -`Rugosity_SD <- raster::extract(Rugosity, ., fun = sd, na.rm=T)`,
                -`Max_Rugosity <- raster::extract(Rugosity, ., fun = max, na.rm=T)`
  )

#clearing space
gc()


#-----------------------------------------------------------------------------#
########################## ADDING CLIMATE VARIABLES ###########################

####-------------------LONG TERM MEAN ANNUAL PRECIPITATION ----------------####


#Adding precipitation
setwd(prec_dir)
MAP <- raster("MAP.tif")
#checking
crs(MAP)

#Matching projections
Buf30m <- spTransform(Buf30m, proj4string(MAP))
#checking
crs(Buf30m)

#Extracting precipitation (no need fo functions since resolution is larger than buffer)
Buf30m <- Buf30m %>%
  mutate(
    MAP = raster::extract(MAP,.,fun=mean, na.rm=T)
  )

#making numeric
Buf30m$MAP <- as.numeric(Buf30m$`MAP = raster::extract(MAP,.,fun=mean, na.rm=T)`)

#removing bad names
Buf30m <- Buf30m %>%
  dplyr::select(-`MAP = raster::extract(MAP,.,fun=mean, na.rm=T)`)

#cleaning
rm(MAP)
gc()

####-------------------------- CHRONIC WIND EXPOSURE ----------------------####

setwd(wind_dir)
wind <- raster("WindBinary_PR_30m_15deg_500mdist.tif")
#checking
crs(wind)

#Matching projections
Buf30m <- spTransform(Buf30m, proj4string(wind))
#checking
crs(Buf30m)

#Extracting mean and sd exposure
Buf30m <- Buf30m %>%
  mutate(
    Wind_Binary = raster::extract(wind,., fun=mean, na.rm=T),
    Wind_sd = raster::extract(wind,., fun=sd, na.rm=T)
  )

#cleaning
gc()

#making numeric
Buf30m$Wind_Binary <- as.numeric(Buf30m$`Wind_Binary = raster::extract(wind,., fun=mean, na.rm=T)`)
Buf30m$Wind_sd <- as.numeric(Buf30m$`Wind_sd = raster::extract(wind,., fun=sd, na.rm=T)`)

#removing bad names
Buf30m <- Buf30m %>%
  dplyr::select(-`Wind_Binary = raster::extract(wind,., fun=mean, na.rm=T)`,
                -`Wind_sd = raster::extract(wind,., fun=sd, na.rm=T)`)

#cleaning
rm(wind)
gc()

####---------------------- HURRICANE EXPOSURE ---------------------------####

#loading hurricane hugo
setwd(paste0(hur_dir,"/Hugo"))
hugo <- raster("dblbnd.adf")
#checking
crs(hugo)
plot(hugo)

#Matching projections
Buf30m <- spTransform(Buf30m, proj4string(hugo))
#checking
crs(Buf30m)

#Extracting mean and sd exposure
Buf30m <- Buf30m %>%
  mutate(
    Hugo = raster::extract(hugo,., fun=mean, na.rm=T),
    Hugo_sd = raster::extract(hugo,., fun=sd, na.rm=T),
  )

#cleaning
gc()

#making numeric
Buf30m$Hugo <- as.numeric(Buf30m$`Hugo = raster::extract(hugo,., fun=mean, na.rm=T)`)
Buf30m$Hugo_sd <- as.numeric(Buf30m$`Hugo_sd = raster::extract(hugo,., fun=sd, na.rm=T)`)

#removing bad names
Buf30m <- Buf30m %>%
  dplyr::select(-`Hugo = raster::extract(hugo,., fun=mean, na.rm=T)`,
                -`Hugo_sd = raster::extract(hugo,., fun=sd, na.rm=T)`)

#cleaning
rm(hugo)
gc()

#----------------------------

#loading hurricane georges
setwd(paste0(hur_dir,"/Georges"))
georges <- raster("dblbnd.adf")
#checking
crs(georges)
plot(georges)

#Matching projections
Buf30m <- spTransform(Buf30m, proj4string(georges))
#checking
crs(Buf30m)

#Extracting mean and sd exposure
Buf30m <- Buf30m %>%
  mutate(
    Georges = raster::extract(georges,., fun=mean, na.rm=T),
    Georges_sd = raster::extract(georges,., fun=sd, na.rm=T),
  )

#cleaning
gc()

#making numeric
Buf30m$Georges <- as.numeric(Buf30m$`Georges = raster::extract(georges,., fun=mean, na.rm=T)`)
Buf30m$Georges_sd <- as.numeric(Buf30m$`Georges_sd = raster::extract(georges,., fun=sd, na.rm=T)`)

#removing bad names
Buf30m <- Buf30m %>%
  dplyr::select(-`Georges = raster::extract(georges,., fun=mean, na.rm=T)`,
                -`Georges_sd = raster::extract(georges,., fun=sd, na.rm=T)`)

#cleaning
rm(georges)
gc()


#-----------------------------------------------------------------------------#
######################## ADDING TOPOGRAPHIC VARIABLES #########################

####------------------------------ ELEVATION ------------------------------####

setwd(dem_dir)
#loading elevation map
dem <- raster("DEM_30m_Proj.tif")
#checking
crs(dem)
res(dem)

#Matching projections
Buf30m <- spTransform(Buf30m, proj4string(dem))
#checking
crs(Buf30m)


#Extracting mean and sd elevation
Buf30m <- Buf30m %>%
  mutate(
    Elev_Mean = raster::extract(dem,., fun=mean, na.rm=T),
    Elev_sd = raster::extract(dem,., fun=sd, na.rm=T),
    Elev_Max = raster::extract(dem, ., fun=max, na.rm=T),
    Elev_range = raster::extract(dem, ., fun=range, na.rm=T)
  )

#cleaning up space
gc()

#Making columns numeric and better named
Buf30m$Elev_Mean <- as.numeric(Buf30m$`Elev_Mean <- raster::extract(dem, ., fun = mean, na.rm=T)`)
Buf30m$Elev_sd <- as.numeric(Buf30m$`Elev_sd <- raster::extract(dem, ., fun = sd, na.rm=T)`)
Buf30m$Elev_Max <- as.numeric(Buf30m$`Elev_Max <- raster::extract(dem, ., fun = max, na.rm=T)`)
Buf30m$Elev_range <- as.numeric(Buf30m$`Elev_range = raster::extract(dem, ., fun=range, na.rm=T)`)

#Removing badly named columns
Buf30m <- Buf30m %>%
  dplyr::select(-`Elev_Mean <- raster::extract(dem, ., fun = mean, na.rm=T)`,
                -`Elev_sd <- raster::extract(dem, ., fun = sd, na.rm=T)`,
                -`Elev_Max <- raster::extract(dem, ., fun = max, na.rm=T)`,
                -`Elev_range = raster::extract(dem, ., fun=range, na.rm=T)`
  )


####-------------------------------- SLOPE ------------------------------####

#calculate slope
slope <- raster::terrain(dem,opt="slope",unit="degrees",neighbors=8)

#checking
crs(slope)
res(slope)

#Matching projections
Buf30m <- spTransform(Buf30m, proj4string(slope))
#checking
crs(Buf30m)

#Extracting mean and sd slope
Buf30m <- Buf30m %>%
  mutate(
    Slope = raster::extract(slope,., fun=mean, na.rm=T),
    Slope_sd = raster::extract(slope,., fun=sd, na.rm=T),
  )

#cleaning
gc()

#making numeric
Buf30m$Slope <- as.numeric(Buf30m$`Slope = raster::extract(slope,., fun=mean, na.rm=T)`)
Buf30m$Slope_sd <- as.numeric(Buf30m$`Slope_sd = raster::extract(slope,., fun=sd, na.rm=T)`)

#removing bad names
Buf30m <- Buf30m %>%
  dplyr::select(-`Slope = raster::extract(slope,., fun=mean, na.rm=T)`,
                -`Slope_sd = raster::extract(slope,., fun=sd, na.rm=T)`)

#cleaning
rm(slope)
gc()

####------------------------ TOPOGRAPHIC ROUGHNESS ----------------------####

#Convert to meters not degrees
dem <- raster::projectRaster(from = dem,
                             crs = nad)

#checking
crs(dem)
res(dem)

#calculate topographic ruggedness
TRI <- raster::terrain(dem,
                       opt="TRI")

#checking
crs(TRI)
res(TRI)

#Matching projections
Buf30m <- spTransform(Buf30m, proj4string(TRI))
#checking
crs(Buf30m)

#Extracting mean topographic roughness
Buf30m <- Buf30m %>%
  mutate(
    TRI = raster::extract(TRI,., fun=mean, na.rm=T)
  )

#cleaning
gc()

#making numeric
Buf30m$TRI <- as.numeric(Buf30m$TRI)

#cleaning
rm(TRI,dem)
gc()

####----------------------------- SOIL TYPE -----------------------------####

setwd(soil_dir)
#loading soil shapefile
soil <- readOGR(dsn = paste0(soil_dir),"Soil")

#checking
crs(soil)

#Keeping only soil type
soil <- soil %>%
  dplyr::select(Soil)

#Matching projections
Buf30m <- spTransform(Buf30m, proj4string(soil))
#checking
crs(Buf30m)

#Extracting soil type
Soil <- over(Buf30m, soil)$Soil

Buf30m <- Buf30m %>%
  mutate(
    Soil = Soil
  )

#cleaning
rm(soil,Soil)
gc()

####---------------------- SOIL WATER AVAILABILITY ------------------------####

setwd(soil_dir)
#Loading raster
aws <- raster("aws0_150.tif")

#Matching projections
Buf30m <- spTransform(Buf30m, proj4string(aws))
#checking
projection(Buf30m)
projection(aws)

#Extracting mean and sd water availability
Buf30m <- Buf30m %>%
  mutate(
    AWS_Mean = raster::extract(aws,., fun=mean, na.rm=T),
    AWS_sd = raster::extract(aws,., fun=sd, na.rm=T),
  )

#cleaning
gc()

#making numeric
Buf30m$AWS_Mean <- as.numeric(Buf30m$`AWS_Mean = raster::extract(aws,., fun=mean, na.rm=T)`)
Buf30m$AWS_sd <- as.numeric(Buf30m$`AWS_sd = raster::extract(aws,., fun=sd, na.rm=T)`)

#removing bad names
Buf30m <- Buf30m %>%
  dplyr::select(-`AWS_Mean = raster::extract(aws,., fun=mean, na.rm=T)`,
                -`AWS_sd = raster::extract(aws,., fun=sd, na.rm=T)`)

#cleaning
rm(aws)
gc()

#-----------------------------------------------------------------------------#
################################ SAVING DATAFRAMES ############################

#saving shapefile with all information
writeOGR(Buf30m,
         dsn=file.path(point_dir),
         layer = "Buf30m_Full_Shapefile",
         driver = "ESRI Shapefile")

#writing csv of the data from the shapefile
write_csv(Buf30m@data,paste0(local_path,"/Data/Full_Data_For_Analyses.csv"))

#cleaning
rm(list=ls())
gc()


