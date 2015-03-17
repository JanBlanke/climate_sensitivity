### Name: data_preparation.R
### Author: Jan Blanke
### Description: Preparation data for svd analysis
####################################################

library(sp)
library(raster)
library(rasterVis)
library(plyr)

### Settings
output <- "yield" # lai, cpool, cflux, cmass
mkdir <- "yield_2069_2098"
years <- 2069:2098 # 2050, 2098, 2021:2050, 2069:2098
dir.create(file.path("/media/jan/AddData/Simulations_processed/30y_avg_all/", paste(mkdir)))

#col.extr <- 78 # cmass natural_sum: 83, lai total: 78, cmass total: 78, cpool total: 9, NEE: 10, Veg cflux: 4
#if (length(years) > 1) col.extr <- 3 # ddply creates new dataframe with 3rd column beeing the averaged output
 
### Initialize
df.list <- list()

### Loop through models/rcps
loopstr <- c("CMIP5_CNRM_rcp45", "CMIP5_CNRM_rcp85", "CMIP5_MOHC_rcp45", "CMIP5_MOHC_rcp85", "CMIP5_ICHEC_rcp45", "CMIP5_ICHEC_rcp85", "CORDEX_CNRM_rcp45", "CORDEX_CNRM_rcp85", "CORDEX_MOHC_rcp45", "CORDEX_MOHC_rcp85", "CORDEX_ICHEC_rcp45", "CORDEX_ICHEC_rcp85")

for (i in loopstr) {

  ### Read and prepare files 
  
  ######################
  ### 10 minute data ###
  ######################
  
  setwd("/home/jan/LPJ_GUESS/Paper_1/crop_10min_run_multicore/sims/")                                                                                                                              
  list.files()
  
  ### Read CORDEX 10 arcmin climate data and fracs
  file.to.read <- list.files(pattern = paste(i, "_", output , sep=""))[1]
  
  temp.out <- read.table(file.to.read, header=T, sep="")
    
  ### Convert to spatial objects                                                                                                                              
  # Subset according to years
  years.tmp <- years - 1970 + 500
  temp.out <- subset(temp.out, temp.out[3] >= min(years.tmp) & temp.out[3] <= max(years.tmp))
  
  if (length(years) > 1 ) {
    temp.out <- ddply(temp.out, c(.(Lon), .(Lat)), summarise, out.avg = mean(TeWW, na.rm = TRUE))
    print("30 y average calculated")
    flush.console()
  }
  
  # make it spatial
  col.extr <- which(colnames(temp.out) == "Total") #overwrite col.extr
  if (length(years) > 1) col.extr <- 3
  temp.out <- SpatialPointsDataFrame(temp.out[, 1:2], as.data.frame(temp.out[, c(col.extr)]), proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  
  ######################
  ### 30 minute data ###                                                                                                                              
  ######################
  
  setwd("/media/jan/AddData/Simulations/crop_r3799_netcdf_30min")                                                                                                                              
  list.files()                                                                                                                              
  file.to.read <- list.files(pattern = paste(i, "_", output , sep=""))[1]
  
  ### Read CORDEX 30 arcmin LPJ output files
  if (grepl("csv", file.to.read)) {
    temp.out.2 <- read.table(file.to.read, header=T, sep=",")                                                                                                                              
  } else {
    temp.out.2 <- read.table(file.to.read, header=T, sep="")                                                                                                                              
  }  
  
  ### Convert to spatial objects
  # Subset according to years
  years.tmp <- years - 1951 + 500
  temp.out.2 <- subset(temp.out.2, temp.out.2[3] >= min(years.tmp) & temp.out.2[3] <= max(years.tmp))
  
  if (length(years) > 1 ) {
    temp.out.2 <- ddply(temp.out.2, c(.(Lon), .(Lat)), summarise, out.avg = mean(TeWW, na.rm = TRUE))
  } 
  
  
  col.extr <- which(colnames(temp.out.2) == "Total") #overwrite col.extr
  if (length(years) > 1) col.extr <- 3
  temp.out.2 <- SpatialPixelsDataFrame(temp.out.2[, 1:2], as.data.frame(temp.out.2[, col.extr]), proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  temp.out.2 <- raster(temp.out.2) 
  
  #########################                                                                                                                              
  ### Match resolutions ###                                                                                                                              
  #########################                                                                                                                              
  ### Extract low res pixel values and store in last column - matching both resolutions                                                                                                                              
  temp.out <- extract(temp.out.2, temp.out, sp=T) # take 10 min grid, extract values from 30 min grid                                                                                                                           
    
  ### Create empty grid with resolution of 0.17                                                                                                                              
  #r <- raster(ncols=(bbox(temp.out)[1, 2] - bbox(temp.out)[1, 1]) / 0.171, nrows=(bbox(temp.out)[2, 2] - bbox(temp.out)[2, 1]) / 0.171, xmn=bbox(temp.out)[1, 1], xmx=bbox(temp.out)[1,2], ymn=bbox(temp.out)[2, 1], ymx=bbox(temp.out)[2, 2], crs=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  
  ### Rasterize point data - only 1 layer
  #rast.10min <- rasterize(as.data.frame(temp.out)[, 1:2], r, as.data.frame(temp.out)[, 3], fun=mean)
  #rast.30min <- rasterize(as.data.frame(temp.out)[, 1:2], r, as.data.frame(temp.out)[, 4], fun=mean)
  #rast.brick <- brick(rast.10min, rast.30min)
  #names(rast.brick) <- c("min10", "min30")
  
  ### Save data
  # Save dfs in list
  idx <- which(loopstr == i)
  df.list[[idx]] <- temp.out
  
  # Save raster
  #setwd(paste("/media/jan/AddData/Simulations_processed/", mkdir, sep=""))
  #writeRaster(rast.10min, paste(i,"_", mkdir, "_10min.tiff", sep=""))
  #writeRaster(rast.30min, paste(i,"_", mkdir, "_30min.tiff", sep=""))
  
  # Make plots
  #pdf(paste(i,"_", mkdir, ".pdf", sep=""), width=7, height=5)
  #terrTheme <- rasterTheme(region=rev(terrain.colors(100)))
  #print(levelplot(rast.brick, par.settings=terrTheme, margin=F))
  #print(bwplot(rast.brick))
  #print(histogram(rast.brick))
  #print(densityplot(rast.brick))
  #dev.off()
  
  cat("iteration:", i, " \n")
  flush.console()
  
} # end loop

### Store everything in one data frame
df.out <- lapply(df.list, as.data.frame) # convert to df

SortFun <- function(x) x[ order(x[, 1], x[, 2]), ]
df.out <- lapply(df.out, SortFun) # Sort Lons, Lats

df.out <- do.call(cbind, df.out)
data.to.use <- sort(c(seq(3, 48, 4), seq(4, 48, 4)))

df.coords <- df.out[, - (data.to.use)] # coordinates
coords <- df.coords[, 1:2]
df.out <- df.out[, data.to.use] # delete coordinates

col.names <- loopstr[sort(rep(1:12, 2))]
for (i in seq(1, 24, 2)) col.names[i] <- paste(col.names[i], "_10min", sep="")
for (i in seq(2, 24, 2)) col.names[i] <- paste(col.names[i], "_30min", sep="")
names(df.out) <- col.names  

setwd(paste("/media/jan/AddData/Simulations_processed/30y_avg_all/", mkdir, sep=""))
getwd()
save(df.out, file = "df_out.RData")
save(coords, file = "coords.RData")
