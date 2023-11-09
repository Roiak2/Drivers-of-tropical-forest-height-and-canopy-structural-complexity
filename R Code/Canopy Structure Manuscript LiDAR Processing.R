##############################################################################
# Compute forest metrics from LiDAR point clouds.
# Metrics include canooy height model and canopy rugosity
# 
# If you want to view a point cloud with the PointCloudViewer package:
# las <- readLAS(point_cloud_file)
# plot(las, backend = "pcv")
#
# If you want to view coverage of a LASCatalog:
# plot(ctg)
#
# Code by Lora Murphy
##############################################################################

library(lidR)
library(rgdal)
library(raster)

rm(list=ls(all=TRUE))
dtm_dir <- "C:/Users/uriarte/Documents/LoraMurphy/LiDAR/DTM/"
tiles_dir <- "C:/Users/uriarte/Documents/LoraMurphy/LiDAR/gis data"
#tiles_dir <- "C:/users/lora/documents/mu/lidar/2019 forest metrics"
clean_dir <- "U:/USGS LiDAR/Cleaned/"
norm_dir <- "C:/Users/uriarte/Documents/LoraMurphy/LiDAR/Normalized/"
chm_dir <- "C:/Users/uriarte/Documents/LoraMurphy/LiDAR/CHM/"
rug_dir <- "C:/Users/uriarte/Documents/LoraMurphy/LiDAR/Rugosity/"

# This option prevents an error message about exported globals being too big,
# from the "future" package, which is used by lidR. The exported globals 
# appears to be important in the context of parallel processing, and if you
# were using this script in that situation, I believe the next line is a
# terrible idea. But running sequentially on a single machine, this lets you
# proceed.
options(future.globals.maxSize = Inf)

#-----------------------------------------------------------------------------#
# Divide the island up into chunks for processing. 12X7 tiles is what is 
# working OK at the moment. Tiles are approx 0.0135 degrees tall, and 0.01415
# degrees wide
#-----------------------------------------------------------------------------#
tile_shp <- readOGR(dsn = tiles_dir,
                    layer = "lidar_tiles",
                    stringsAsFactors = F)

# Get rid of the islands
minX <- sapply(as.character(tile_shp@data$bndngBx), function(x){
  v <- strsplit(x, ",")[[1]][2]
  v <- strsplit(v, ":")[[1]][2]
  as.numeric(v)
})
tile_shp <- tile_shp[minX > -67.3,]

# Parse out the bounding box
minX <- sapply(as.character(tile_shp@data$bndngBx), function(x){
  v <- strsplit(x, ",")[[1]][2]
  v <- strsplit(v, ":")[[1]][2]
  as.numeric(v)
})
maxX <- sapply(as.character(tile_shp@data$bndngBx), function(x){
  v <- strsplit(x, ",")[[1]][4]
  v <- strsplit(v, ":")[[1]][2]
  v <- gsub("}", "", v)
  as.numeric(v)
})
minY <- sapply(as.character(tile_shp@data$bndngBx), function(x){
  v <- strsplit(x, ",")[[1]][1]
  v <- strsplit(v, ":")[[1]][2]
  as.numeric(v)
})
maxY <- sapply(as.character(tile_shp@data$bndngBx), function(x){
  v <- strsplit(x, ",")[[1]][3]
  v <- strsplit(v, ":")[[1]][2]
  as.numeric(v)
})

# Assign each tile a row and column
xSeq <- c(seq(from = min(minX), to = max(minX), by = 0.01415), max(minX))
ySeq <- c(seq(from = min(minY), to = max(minY), by = 0.0135), max(minY))

tile_shp$row <- NA
tile_shp$col <- NA

for (i in 1:length(xSeq)-1) {
  x <- which(minX >= xSeq[i] & minX <= xSeq[i+1])
  tile_shp$col[x] <- i
}
for (i in 1:length(ySeq)-1) {
  x <- which(minY >= ySeq[i] & minY <= ySeq[i+1])
  tile_shp$row[x] <- i
}

# Group rows and columns into regions
tile_shp$tile_group <- NA
xSeq <- c(0, seq(from=12, to = 121, by = 12))
xSeq[length(xSeq)] <- 121

ySeq <- c(0, seq(from=7, to = 49, by = 7))

count = 1
for (i in 2:length(xSeq)) {
  for (j in 2:length(ySeq)) {
    x <- which(tile_shp$row > ySeq[j-1] & tile_shp$row <= ySeq[j] &
                 tile_shp$col > xSeq[i-1] & tile_shp$col <= xSeq[i])
    tile_shp$tile_group[x] <- count
    count <- count + 1
  }
}

# Hand tuning - combine 57 and 58
tile_shp$tile_group[tile_shp$tile_group == 58] <- 57
x <- which(tile_shp$tile_group > 57)
tile_shp$tile_group[x] <- tile_shp$tile_group[x]-1
rm(i, j, x, xSeq, ySeq, minX, minY, maxX, maxY)

#-----------------------------------------------------------------------------#
# Make a li'l shapefile with a polygon for each tile group
# This method didn't actually work because apparently the shapefiles of the 
# downloaded tiles don't exactly line up. I had to use ArcGIS Integrate to
# make the polygons meet, then Dissolve to create the tile group polygons.
#-----------------------------------------------------------------------------#
#library(maptools)
# Merge polygons by ID
#tg <- tile_shp$tile_group
#tiles <- unionSpatialPolygons(tile_shp, tg, threshold = 0.1)



#-----------------------------------------------------------------------------#
# Work with individual tile groups
#-----------------------------------------------------------------------------#
unique_tile_groups <- unique(tile_shp$tile_group)

for (tg in unique_tile_groups) {
  
  tg_dir <- paste0(norm_dir, "tg_", tg, "/")
  if (!dir.exists(tg_dir)) {
    dir.create(tg_dir)
  }
  
  #---------------------------------------------------------------------------#
  # First step in normalizing the point cloud: merge the DTMs of a tile group
  # (with surrounding tiles for safety) into a single file
  #---------------------------------------------------------------------------#
  
  # Break out just this group of tiles
  gtiles <- tile_shp[tile_shp$tile_group == tg,]
  
  # Figure out row and column range
  minRow <- min(gtiles$row)
  maxRow <- max(gtiles$row)
  minCol <- min(gtiles$col)
  maxCol <- max(gtiles$col)
  
  # Southern border outside of tile - allow that there might not be any
  borderTiles <- NULL
  bt <- tile_shp[tile_shp$row == (minRow - 1) & tile_shp$col %in% (minCol - 1):(maxCol + 1),]
  if (nrow(bt) > 0) borderTiles <- rbind(borderTiles, bt@data)
  
  # Northern border outside of tile
  bt <- tile_shp[tile_shp$row == (maxRow + 1) & tile_shp$col %in% (minCol - 1):(maxCol + 1),]
  if (nrow(bt) > 0) borderTiles <- rbind(borderTiles, bt@data)
  
  # Eastern and western
  bt <- tile_shp[tile_shp$col == (minCol - 1) & tile_shp$row %in% minRow:maxRow,]
  if (nrow(bt) > 0) borderTiles <- rbind(borderTiles, bt@data)
  bt <- tile_shp[tile_shp$col == (maxCol + 1) & tile_shp$row %in% minRow:maxRow,]
  if (nrow(bt) > 0) borderTiles <- rbind(borderTiles, bt@data)
  
  # Get the core filename of each tile
  allTiles <- c(borderTiles$dwnLURL, gtiles$dwnLURL)
  if (any(duplicated(allTiles))) {
    warning(paste0("Dups on tile group", tg))
  } 
  fileCore <- sapply(allTiles, function(x) {
    v <- strsplit(x, "/")[[1]]
    v <- v[length(v)]
    v <- strsplit(v, "_")[[1]][6]
  })
  rm(bt, borderTiles, allTiles, minRow, maxRow, minCol, maxCol)
  
  # Build the list of DTMs to merge. This was a little difficult. One area of
  # the island had missing DTMs - the LiDAR files were there but not the DTMs.
  # Seems unlikely but there it was. So I had to comment out the check below
  # for the tile group affected (52). I elected not to build my own DTM because
  # algorithms do differ and I didn't want it to be inconsistent.
  # Update - decided to build my own because that hole is right over one of
  # our plots.
  args <- list()
  count <- 1
  used_it <- rep(F, length(fileCore))
  for (i in 1:length(fileCore)) {
    dtmf <- paste0(dtm_dir, "USGS_NED_OPR_PR_PuertoRico_2015_",
                   fileCore[i], "_IMG_2018/USGS_NED_OPR_PR_PuertoRico_2015_",
                   fileCore[i], "_IMG_2018.img")
    if (!file.exists(dtmf)) {
      # Check to see if the lidar file exists either - if it doesn't,
      # ignore the fact that this is missing
      if (file.exists(paste0(clean_dir, "USGS_LPC_PR_PuertoRico_2015_", 
                             fileCore[i], "_LAS_2018_1_1.laz"))) {
        # Calculate a DTM
        las <- readLAS(paste0(clean_dir, "USGS_LPC_PR_PuertoRico_2015_", 
                              fileCore[i], "_LAS_2018_1_1.laz"))
        dtmtemp <- grid_terrain(las, res=1, algorithm=tin())
        args[[count]] <- projectRaster(dtmtemp, crs=proj4string(args[[1]]))
        count <- count + 1
        used_it[i] <- T
        
      }
    } else {
      args[[count]] <- raster(dtmf)
      count <- count + 1
      used_it[i] <- T
    }
  }
  alldtm <- do.call(raster::merge, args)
  #writeRaster(alldtm, filename = "tg_52_dtm.tiff", format = "GTiff")
  #plot(gtiles)
  #text(coordinates(gtiles), labels=fileCore)
  
  
  #---------------------------------------------------------------------------#
  # Normalize the point clouds. Normalization is the process of removing 
  # elevation from the point cloud height, leaving just height-above-ground
  #---------------------------------------------------------------------------#
  
  # Create a vector of clipped/cleaned files in the tile group. As you can see
  # I changed the way to assemble the list of files - this makes sure that 
  # every LiDAR tile is matched by a DTM. 
  #fileCore <- fileCore[used_it]
  fileCore <- sapply(gtiles$dwnLURL[gtiles$inForest==1], function(x) {
    v <- strsplit(x, "/")[[1]]
    v <- v[length(v)]
    v <- strsplit(v, "_")[[1]][6]
  })
  in_laz <- paste0(clean_dir, "USGS_LPC_PR_PuertoRico_2015_", fileCore, "_LAS_2018_1_1.laz")
  for (i in 1:length(in_laz)) {
    if (!file.exists(in_laz[i])) stop(paste0("Missing files in tile group", tg))
  }
  
  
  #---------------------------------------------------------------------------#
  # Set processing options and GO
  # 
  # I had to play around heavily with the chunk size and the buffer. Going for
  # chunks being full files didn't work because this ran out of memory. (This
  # is why this is being done in a LASCatalog, which can chunk them up.) But
  # if you don't choose chunk and buffer right, it makes chunks which don't
  # intersect with the point cloud and gets mad when there's nothing in them.
  #---------------------------------------------------------------------------#
  ctg <- readLAScatalog(in_laz)
  opt_output_files(ctg) <- paste0(norm_dir, "tg_", tg, "_{ID}_f")
  opt_chunk_size(ctg) <- 450
  ctg@chunk_options$buffer <- 0.01
  opt_laz_compression(ctg) <- T
  lasnormalize(ctg, alldtm, na.rm = T) 
}


# In case you need to play around or verify results - here's a good place
# to put it
check_on_normalization <- F
if (check_on_normalization) {
  
  # See what tile groups we have represented in the output files
  fnorm <- c(list.files(norm_dir), list.files(paste0(norm_dir, "proc")))
  ftg <- sapply(fnorm, function(x) {strsplit(x, "_")[[1]][2]})
  present_tile_groups <- sort(as.numeric(unique(ftg)))
  todotg <- unique_tile_groups[which(!unique_tile_groups %in% present_tile_groups)]
}
#-----------------------------------------------------------------------------#
# Next step: go run CHM with LASTools script.
# CHM = canopy height model, which is a raster of first returns from LiDAR,
# and thus presumably a map of the top of the canopy. We have removed noise
# from the data, but there will still be artifacts from things like clouds so
# you have to take the data with a grain of salt. I processed this at a 1X1
# m resolution. That means that there are occasional NoData holes where there
# weren't enough points to calculate a height. If the CHM itself was something
# I really needed to work with, I could interpolate these away. I chose not
# to at this point because I want to know when data is actually missing.
# 
# Once the LASTools batch script is done, come back here...
#-----------------------------------------------------------------------------#  





#-----------------------------------------------------------------------------#
# Some CHM post-processing
# We want to assemble all of the little chunks into larger, easier-to-work-with
# files; and we will clean out artifacts by bounding CHM values between 0 and
# 40 meters. Since we didn't make the DTM, the interpolation used for it might
# have caused negative values to be possible in some spots. We also need to set
# a max value for artifacts like clouds. Maria requests 50 m
# 
# We will assemble the CHM pieces into one file per tile group
#-----------------------------------------------------------------------------#
setwd(chm_dir)

# Access the list of directories
tg_groups <- list.files(path =".") 
tg_groups <- tg_groups[!grepl("tif", tg_groups)]

# Process them together into subtiles of the island, one at a time
for (tg in tg_groups) {
  chm_filename <- paste0(tg, "_chm.tif")
  setwd(chm_dir)
  if (!file.exists(chm_filename)) {
    setwd(paste0(chm_dir, tg))
    files <- list.files(".", pattern = "*.tif")
    if (length(files) > 0) {
      args <- list()
      for (i in 1:length(files)) {
        
        # Bound the raster values between 0 and 50
        ras <- raster(files[i])
        values(ras)[values(ras) < 0 ] <- 0
        values(ras)[values(ras) > 50] <- 50
        
        args[[i]] <- ras
      }
      
      # Merge them all together
      chm <- do.call(raster::merge, args)
      setwd(chm_dir)
      
      # Write out the finished larger tile group file
      writeRaster(chm, filename = paste0(tg, "_chm.tif"), format = "GTiff", overwrite = T)
    }
  }
}


#-----------------------------------------------------------------------------#
# Rugosity
# Rugosity is just standard deviation of canopy height.
# As previously mentioned, the CHM wasn't smoothed or interpolated; it will have
# NoData holes. I would rather be aware of the amount of data available rather
# than interpolating it away when calculating rugosity. So I will allow up to 
# 10% NoData values in the rugosity calc.
# 
# Rugosity is calculated in a 15X15 m moving window, and will also be 1X1 m
# in resolution, like CHM.
#-----------------------------------------------------------------------------#

# A function to do a somewhat-NoData-tolerant standard deviation 
rugfunc <- function(x) {
  if (sum(is.na(x)) > 23) {
    return(NA)
  } 
  sd(x, na.rm = T)
}

# Get all the CHM files
setwd(chm_dir)
list.files(".", pattern = "chm.tif")
for (i in 1:length(files)) {
  rug_filename <- gsub("chm", "rugosity", files[i])
  setwd(rug_dir)
  if (!file.exists(rug_filename)) {
    setwd(chm_dir)
    ras <- raster(files[i])
    
    # This is a moving window calculation from the raster package
    rug <- focal(ras, matrix(1, nrow=15, ncol=15), fun=rugfunc)
    
    setwd(rug_dir)
    writeRaster(rug, filename = rug_filename, format = "GTiff", overwrite = T)
  }
}
