# Preprocess data from monthly data to annual data
#
library(raster)
library(ncdf4)

# Use one dataset for a test
# File names --------------------------------------------------------------

infilename1 <- '/Volumes/Land/Data_proj_wuecu/iWUE_monthly_full_05d.nc' # 
infilename2 <- '/Volumes/Land/Data_proj_wuecu/GPP_monthly_full_05d.nc' # GPP with no water stress

outfilepath1 <- '/Volumes/Land/Data_proj_wuecu/'

# Load data ---------------------------------------------------------------

data1 <- brick(infilename1)
data2 <- brick(infilename2)

# Estimate annual values (mean) -------------------------------------------

# start from 1982 to 2010
index1 <- 7
index2 <- (2010-1981)*12+6 # nlayers(data1)
data1select <- data1[[index1:index2]]
data2select <- data2[[index1:index2]]

# estimate annual mean values
# Method 1
# yearindex <- rep(1982:2010, each=12)
# data1ann <- stackApply(data1select, yearindex, fun = mean, na.rm = T)

# Method 2
# consider the number of monthly data in each year
ann_mean_fun <- function(x) {
  num <- length(x)/12 # the number of years
  dim(x) <- c(12, num)
  re.data <- vector('numeric', num)
  for (i in 1:num) {
    # monthly data for each year
    temp1 <- x[,i]
    temp1 <- temp1[complete.cases(temp1)]
    if (length(temp1) < 3) {
      re.data[i] <- NA
    } else {
      re.data[i] <- mean(temp1)
    }
  }
  # return values
  return(re.data)
}
data1ann2 <- calc(data1select, fun = ann_mean_fun)
data2ann2 <- calc(data2select, fun = ann_mean_fun)


# Save the annual mean data -----------------------------------------------

writeRaster(data1ann2, paste0(outfilepath1,'iWUE_annual_full_05d.nc'), 
            varname = 'iWUE_annual', overwrite=T)
writeRaster(data2ann2, paste0(outfilepath1,'GPP_annual_full_05d.nc'), 
            varname = 'GPP_annual_nows', overwrite=T)
