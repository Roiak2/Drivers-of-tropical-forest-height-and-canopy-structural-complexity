###############################################
#### PUERTO-RICO PREVALENT WIND DIRECTION #####
###############################################

# Roí Ankori-Karlinsky, June 2020

#---------------------------------------------#

############### THIS SCRIPT #####################

# 1. Plots wind roses from global wind atlas
#    and data from 6 weather stations 2000-2016
# 2. Uses directions to create
#    a wind exposure raster from DEM

############### PACKAGES ########################

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

############### DIRECTORIES #####################

#Directory locations will change depending on one's computer

#File path for most files
local_path <- 'C:/Users/roiak/Documents/Damage Project/Rugosity'

#Directory for GWA data
gwa_dir<- paste0(local_path,'/Data/Wind Models/Global Wind Atlas CFDDA 100m Res Data/')

#Directory for weather station data
station_dir <- paste0(local_path,'/Data/Wind Models/Iowa State Wind Station Data/')

#Where local shapefiles are
shp_dir <- 'C:/Users/roiak/Documents/Damage Project/Rugosity/Data/Shapefiles'

#Where DEM is for exposure calculation
dem_dir <- paste0(shp_dir,"/DEM")

#Where forest age data is (used for extent drawing)
age_dir <- paste0(shp_dir,"/Forest age")

#Where wind exposure raster will be outputted
wind_dir <- paste0(shp_dir,"/Wind Exposure")

############### LOADING DATA ####################

#Loading GWA Data
gwa_data <- read_csv(paste0(gwa_dir,"Merged_Model_Wind.csv")) %>%
  filter(Site=="Island")

#Loading Weather Station Data
station_data <- read_csv(paste0(station_dir,"WindData from Stations With Speed.csv"),
                         col_types = list(.default = col_double(), 
                                          DateTime = col_character(),
                                          Station = col_character())) %>%
  #arranging by station and date
  arrange(Station,DateTime,WindDirection) %>%
  filter(!is.na(WindDirection),
         WindDirection != "null") %>%
  #Adding site direction from each weather station
  mutate(Site_Direction = ifelse(
    Station=="TJBQ","NorthWest",ifelse(
      Station=="TJIG","North",ifelse(
        Station =="TJSJ", "NorthEast",ifelse(
          Station=="TJNR","East",ifelse(
            Station=="TJPS","South",ifelse(
              Station=="TJMZ","West",NA)))))
  ),
  #Dividing by 30-degree sectors
  sector = case_when(
    WindDirection %in% 0:29 ~ 1,
    WindDirection %in% 30:59 ~ 2,
    WindDirection %in% 60:89 ~ 3,
    WindDirection %in% 90:119 ~ 4,
    WindDirection %in% 120:149 ~ 5,
    WindDirection %in% 150:179 ~ 6,
    WindDirection %in% 180:209 ~ 7,
    WindDirection %in% 210:239 ~ 8,
    WindDirection %in% 240:269 ~ 9,
    WindDirection %in% 270:299 ~ 10,
    WindDirection %in% 300:329 ~ 11,
    WindDirection %in% 330:360 ~ 12
  )
  )

#load forest age raster of Puerto Rico to draw extent
setwd(age_dir)
age <- raster("PR_forest_age.img")
#save projection
right_crs <- crs(age)
rm(age)
gc()

#uploading 30m DEM for entirety of Puerto Rico
setwd(dem_dir) #setting working directory
dem <- raster("DEM_30m_Proj.tif")
#save projection
dem_crs <- crs(dem)
dem_crs

#checking
plot(dem)
res(dem)
extent(dem) 

#clearing up memory
gc()

#matching projections to calculate in meters
dem <- projectRaster(dem,crs=right_crs)
#checking
crs(dem)
extent(dem)
res(dem)

#removing unnecessary stuff
rm(age_dir,dem_dir)
#clearing up memory
gc()

 

#------------------------------------------------

############### WIND-ROSES ######################

#Plot model data wind direction for entire island (windrose)
a <- gwa_data %>%
  ggplot(aes(x=center_degree,y=percentage))+
  geom_bar(stat="identity")+
  labs(title="a)",
       x="Compass Direction (Degrees)",
       y="% of Wind")+
  scale_x_continuous(breaks=seq(0,360,30))+
  scale_y_continuous(breaks=seq(0,40,10))+
  geom_hline(yintercept = 10, linetype="dashed")+
  coord_polar(theta="x",start=-0.25)+
  theme_bw()+
  theme(plot.title=element_text(colour="black",face="bold",hjust=-0.1,size=18),
        axis.title = element_blank(),
        axis.text = element_text(size=20,colour="black"),
        panel.border = element_blank(),
        panel.grid.major = element_line(size=1),
        panel.grid.minor = element_blank())

#Plot station data average wind speeds by direction for entire island (windrose)
b <- station_data %>%
  #removing years after study site and hurricane season
  mutate(Date = lubridate::as_datetime(station_data$DateTime)) %>%
  mutate(Month = lubridate::month(Date),
         Year = lubridate::year(Date)) %>%
  filter(Year < 2017,
         Month != 6,
         Month != 7,
         Month != 8,
         Month != 9,
         Month != 10,
         Month != 11) %>%
  group_by(sector) %>%
  summarize(Total_in_sector = sum(sector),
            Wind_Speed = mean(WindSpeed,na.rm=T),
            Max_WindSpeed = max(WindSpeed, na.rm=T)) %>%
  mutate(Total_in_island = sum(Total_in_sector))%>%
  mutate(percent = (Total_in_sector/Total_in_island)*100) %>%
  mutate(Wind_Speed = Wind_Speed * 2.2369,
         Max_WindSpeed = Max_WindSpeed * 2.2369) %>% #converting to m/s
  mutate(sector = sector*30-30) %>%
  arrange(sector)%>% 
  ggplot(aes(x=sector,y=percent,fill=Wind_Speed))+
  geom_bar(stat="identity")+
  scico::scale_fill_scico(breaks=seq(0,25,7.5),
                          palette = "lajolla")+
  labs(title="b)",
       x="Compass Direction (Degrees)",
       y="% of Wind",
       fill="Mean Wind Speed (m/s)")+
  scale_x_continuous(breaks=seq(0,360,30))+
  scale_y_continuous(breaks=seq(0,30,10))+
  geom_hline(yintercept = 10, linetype="dashed")+
  coord_polar(theta="x",start=-0.25)+
  theme_bw()+
  theme(plot.title=element_text(colour="black",face="bold",hjust=-0.1,size=18),
        axis.title = element_blank(),
        axis.text = element_text(size=20,colour="black"),
        panel.border = element_blank(),
        panel.grid.major = element_line(size=1),
        panel.grid.minor = element_blank())

# extract the legend from one of the plots
legend <- get_legend(
  # create some space to the left of the legend
  b+theme(legend.text = element_text(size=24),
          legend.title = element_text(size=18,face="bold"),
          legend.box.margin = margin(-10, 0,-10,0)
  )
)

#combine plots
fig.s1 <- cowplot::plot_grid(a+theme(legend.position = "none"),
                             b+theme(legend.position = "none"),
                             ncol=2,align="hv")+
  theme(plot.margin = margin(t=0,r=0.75,b=0,l=2.75,"cm"))



#output
annotate_figure(fig.s1,
                bottom = text_grob("Compass Direction (degrees)", 
                                   color = "black",
                                   face="bold",
                                   size=24,
                                   vjust=-10.5),
                left = text_grob("% Wind", 
                                 color = "black",
                                 face="bold",
                                 size=24,rot=90,
                                 vjust=0.5),
                right = legend
)


#--------------------------------------------------------#

############### WIND EXPOSURE CALCULATION ################

# Based on EXPOS model by Boose et al., 1994:
# Boose, E. R., D. R. Foster, and M. Fluet. 1994. 
# Hurricane Impacts to Tropical and Temperate Forest Landscapes. 
# Ecological Monographs 64:369-400.
#
# Adapated to R by Dr. Brian Buma: 
# http://www.brianbuma.com/news/2017/12/29/quantifying-wind-exposure-in-r 

####--------------- EXPOS FUNCTION ------------------####


#    Function

#' @param dem # The DEM raster
#' @param deflect # deflection angle of wind after hitting barrier
#' @param angles # the angles of prevalent wind directions
#' @param max.dist # how far up barrier to look for effect

# First, the following code quickly calculates - 
# based on user supplied wind direction, deflection angles, 
# and search distances (how far upwind a barrier should matter) - 
# relative exposure:

windout.iter <- function(dem, deflect, angles, max.dist) {
  #note that the dem raster must be in planar coordinates
  
  #for smaller datasets:
  #dem <- readAll(dem)     #if in memory
  res <- res(dem)
  
  # do not ignore the case where x and y resolution are not equal
  #stopifnot(all.equal(res[1], res[2]))
  xr <- res[1]
  
  #number of distances to check, basically goes every other cell for speed.
  num.dist <- round(max.dist / xr / 2)
  distance <- seq(xr, max.dist, length.out=num.dist) #distance
  result <- list() #empty list to be filled with results
  j <- 1 #j is initially one and then increases and will iterate over 
  
  #iterate over deflection angles
  for (d in deflect) {
    midrow <- cellFromRow(dem,rownr=1)    #note this does the top one
    elev <- extract(dem,midrow) #get elevation 
    coords <- xyFromCell(dem, midrow) #get coordinates
    
    radangle <- (angles+90) * pi/180  #convert wind exposed angles to radians.
    dcosangle <- -cos(radangle) * distance #calculate cosine
    dsinangle <- sin(radangle) * distance #calculate sine
    x <- apply(coords[,1,drop=FALSE], 1, function(j) j + dcosangle) #get x value from cosine
    y <- apply(coords[,2,drop=FALSE], 1, function(j) j + dsinangle) #get y value from sine
    xy <- cbind(as.vector(x), as.vector(y))
    
    comp.elev <- extract(dem, xy) #get elevation for deflection angles
    comp.elev <- matrix(comp.elev, ncol=num.dist, byrow=TRUE) #turn to matrix
    comp.elev <- comp.elev - elev #subtract DEM elevation from deflected angle elevation
    comp.elev <- t(t(comp.elev) / distance) #transpose and divide by distance from location to max distance
    #notAllNA <- rowSums(is.na(comp.elev)) != num.dist
    ang <- atan(comp.elev) * (180 / pi) #converting angles
    
    r <- apply(ang,1,max) #get maximum
    
    r <- r<=d #maximum if it's smaller or equal to deflection angle
    result[[j]] <- r*1 #add to list of results
    j <- j+1 #make sure next j through iteration is one higher than current
  }
  
  output <-simplify2array(result) #simplify result array
  output <- apply(output,1,sum) #get the sum
  output <- output+1 #add one
  outputs <- list() #make list vector
  outputs[[1]] <- output #add the exposure results to it
  outputs[[2]] <- coords #add coordinates
  return(outputs) #return the resultant list
}


#--------------------------------------------------------#


####--------------- SET PARAMETERS -------------------####
#Store matrix where results will come in
storage <- matrix(nrow=nrow(dem),ncol=ncol(dem))

#set parameters
max.dist <-  500 #500m distance grids used to calculate because of heterogeneous topography
deflect <- c(15) #deflection angles when wind hits to determine how slope protects down-wind areas 
angles <- c(60,90,90,120,150) #prevalent wind direction angles
#Prevalent wind directions are mostly coming from NE, E, bit SE
#Therefore the angles are weighted by those

iter <- 1:nrow(dem) #number of rows in dem
r <- res(dem)[1] #resolution

#clearing up memory
gc()
#Making sure R is utilizing all available RAM
memory.limit(size=50000)

#--------------------------------------------------------#

####--------------- RUNNING EXPOS --------------------####


# This next code then loops through a DEM, 
# first subsetting out an area the size of the max distance (that keeps things fast)
# and then calculating all wind directions desired and averaging the results. 

#iterate over each row in the DEM
for (i in iter) {
  temp.extent <- extent(dem) #get extent
  temp.extent@ymax <- temp.extent@ymax-(i*r)    #*res(dem30)[1] to avoid top #avoid top using resolution
  temp.extent@ymin <- temp.extent@ymax-(i*r+max.dist+r) #avoid top using resolution and maximum distance
  
  temp.dem <- crop(dem,temp.extent) #crop the DEM just to this bit of it
  temp1 <- windout.iter(temp.dem,deflect,angles[1],max.dist)[[1]] #run exposure function on first exposure angle
  temp2 <- windout.iter(temp.dem,deflect,angles[2],max.dist)[[1]] #first angle repeated
  temp3 <- windout.iter(temp.dem,deflect,angles[3],max.dist)[[1]] #second angle
  temp4 <- windout.iter(temp.dem,deflect,angles[4],max.dist)[[1]] #third angle
  temp5 <- windout.iter(temp.dem,deflect,angles[5],max.dist)[[1]] #fourth angle
  
  #temp.coords <- SpatialPoints(temp.coords)
  temp <- apply(cbind(temp1,temp2,temp3,temp4,temp5),1,mean,na.rm=T) #get mean of exposure for this bit of the raster
  #t.loc <- cellFromXY(storage,temp.coords)
  
  storage[i,] <- temp #store in results matrix 
  print(i) #print the iteration so user can see code is running properly
  removeTmpFiles(h=0)   #This is needed large processing jobs, which crash otherwise.
  rm(temp.dem) #remove the cropped DEM 
  gc() #clear memory to help run
}

t <- dem   #this creates a place to put the calculated values
t[] <- storage #adding values of exposure

par(mfrow=c(1,2))   #look at some comparisons of DEM and exposure map
plot(t)
plot(dem)


#--------------------------------------------------------#

####--------- SAVING WIND-EXPOSURE RASTER ------------####

#Writing raster to file
writeRaster(t,file.path(wind_dir,"WindExp_PR_30m_15deg_500mdist"),
            format="GTiff",
            overwrite=T)

#Creating binary raster out of the 0-2 map created above
# create classification matrix
reclass_df <- c(0, 1.95, 0,
                1.95, 2, 1)

# reshape the object into a matrix with columns and rows
reclass_m <- matrix(reclass_df,
                    ncol = 3,
                    byrow = TRUE)

# reclassify the raster using the reclass object - reclass_m
wind_classified <- reclassify(t,
                              reclass_m)

#checking
plot(wind_classified)

#saving
setwd(wind_dir)
writeRaster(wind_classified,
            file.path(wind_dir,"WindBinary_PR_30m_15deg_500mdist"),
            format="GTiff",
            overwrite=T)

#cleaning
rm(list=ls())
gc()