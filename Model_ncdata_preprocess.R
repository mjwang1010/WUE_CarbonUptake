# nc data preprocessing
#
# GIMMS fAPAR data rotation
rm(list = ls())
library(raster)
library(ncdf4)

# File path ---------------------------------------------------------------
outfilepath <- '/Volumes/Land/Data_proj_wuecu/'
infilename1 <- '/Volumes/Land/Data_proj_wuecu/GIMMS_fPAR_monthly_05d.nc' # fAPAR data
infilename2 <- '/Volumes/Land/Data_proj_wuecu/mstmip_driver_global_hd_co2_v1.nc4' # CO2 # Atmospheric CO2

infilepath1 <- '/Volumes/Land/CRUNCEP_post/tair_monthly_mean/'
infilepath2 <- '/Volumes/Land/CRUNCEP_post/press_monthly_mean/'
infilepath3 <- '/Volumes/Land/CRUNCEP_post/swdown_monthly_mean/'
infilepath4 <- '/Volumes/Land/CRUNCEP_post/vpd_monthly/'

# Get data ----------------------------------------------------------------

data1 <- brick(infilename1)
data2 <- brick(infilename2)

# Ta, Press, and Solar radiation data are piled up in a single file

# Data processing ---------------------------------------------------------

### some tests
# data1rot <- rotate(data1[[1]])
# data1t <- t(data1[[1]])
# data1f <- flip(data1t[[1]], 'x')
###
# fPAR need geographical information
data1t <- t(data1)
data1post <- flip(data1t, 'x') # final data
extent(data1post) <- c(-180,180,-90,90)
crs(data1post) <- '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'
index1 <- (2010-1981)*12+6
data1post <- data1post[[1:index1]]

## CO2 concentration data
# from 198107 (the first layer of fPAR data) to 201012
# the first layer represents 170001
totny <- nlayers(data2)/12 # the total number of years
index1 <- (1981-1700)*12+7
index2 <- (2010-1700+1)*12
data2post <- data2[[index1:index2]]

# CRU-NCEP data: temperature, air pressure, and shortwave downward radiation
datapost <- function(path) {
  filename <- list.files(path, full.names = T)
  data <- stack(filename)
  index1 <- (1981-1901)*12 + 7 
  index2 <- (2010-1901+1)*12
  datapost <- data[[index1:index2]]
  #
  return(datapost)
}

data3post <- datapost(infilepath1) # Ta
data4post <- datapost(infilepath2) # Press
data5post <- datapost(infilepath3) # Radiation

data6post <- datapost(infilepath4) # VPD

# infilename3 <- list.files(infilepath1, full.names = T)
# data3 <- stack(infilename3)
# index1 <- (1981-1901)*12 + 7
# index2 <- (2010-1901+1)*12
# data3post <- data3[[index1:index2]]


### Might choose the targeted years later

# Save the data -----------------------------------------------------------


writeRaster(data1post, paste0(outfilepath, 'GIMMS_fPAR_monthly_05d_post.nc'), 
            varname = 'fPAR_monthly', xname = 'longitude', yname = 'latitude', 
            overwrite = T)
writeRaster(data2post, paste0(outfilepath, 'MSTMIP_co2_monthly_05d_post.nc'), 
            varname = 'CO2_monthly', 
            overwrite = T)

writeRaster(data3post, paste0(outfilepath, 'CRUNCEP_Ta_monthly_05d_post.nc'), 
            varname = 'Ta_monthly', 
            overwrite = T)
writeRaster(data4post, paste0(outfilepath, 'CRUNCEP_Press_monthly_05d_post.nc'), 
            varname = 'Press_monthly', 
            overwrite = T)
writeRaster(data5post, paste0(outfilepath, 'CRUNCEP_SWdown_monthly_05d_post.nc'), 
            varname = 'SWdown_monthly', 
            overwrite = T)

writeRaster(data6post, paste0(outfilepath, 'CRUNCEP_VPD_monthly_05d_post.nc'), 
            varname = 'VPD_monthly', 
            overwrite = T)
