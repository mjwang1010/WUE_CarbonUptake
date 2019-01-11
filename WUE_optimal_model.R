# Estimate IWUE based on an optimal model
# 
# Input - monthly/daily/ data
# Ca - CO2 concentrations, umol mol-1
# VPD (D) - vapor pressure deficit, Pa
# Ta - air temprature, K
# z - elevation, km
# (Optional) P - air pressure, Pa; Can be derived from Ta and z
# Output - 
# modelled iWUE, umol mol-1

rm(list = ls())
library(lubridate)
library(raster)
library(ncdf4)


# Set parameters ----------------------------------------------------------
R <- 8.314 # gas constant, J mol-1 K-1
# beta <- 356.51 # a constant
# beta <- 244.03 # the constant derived from Beni Stocker
mfill <- -9999.0 # the filling value

# Define functions --------------------------------------------------------
# coefficient of Rubisco or photorespiratory compensation
# the unit of K is dependent on 'x25'
K_fun <- function(dH, Ta, x25) {
  # dH, J mol-1 K-1
  # Ta, K
  # x25, Pa or umol mol-1
  R <- 8.314
  K <- x25 * exp((dH*(Ta-298.15))/(298.15*R*Ta))
  return(K)
}

# viscosity of water: eta
# Wang Han et. al., 2017, NPlants
eta_fun <- function(Ta) {
  # Ta, K
  A <- 3.719
  B <- 580
  C <- -138
  eta <- 0.001*exp(A + B/(C+Ta))
  return(eta)
}

# O2 partial pressure
# Wang Han et. al., 2017, NPlants
Po_fun <- function(elev) {
  # elev, elevation, km
  # Po, pressure of O2, Pa
  Po <- 21000 * exp(-0.114*elev)
  return(Po)
}

Po_fun1 <- function(Press) {
  # input - press, surface pressure, Pa
  # Po, partial pressure of O2, Pa
  Po <- Press*0.21
  return(Po)
}

# Estimate the IWUE from the environmental variables
# Model_iWUE_fun <- function(Tk, D, Ca, elev)
Model_iWUE_fun1 <- function(Tk, D, Ca, Press) {
  
  # Comments: the model does not consider the impact of soil moisture.
  # Tk, temperature, (K)
  # D, VPD, (Pa)
  # Ca, CO2 concentration, ppm, umol mol-1
  # elev, elevation, km (to estimate O2 partial pressure)
  # Press, atmospheric pressure, Pa (the second method to estimate O2 partial pressure)
  # iWUE, umol mol-1
  
  Po <- Po_fun1(Press) # unit - Pa; or a function of elevation
  beta <- 244.03 # a constant noted by Beni Stocker
  ###
  Ko <- K_fun(36.38*1000, Tk, 27480) # Pa
  Kc <- K_fun(79.43*1000, Tk, 39.97) # Pa
  gamma_star <- K_fun(37830, Tk, 4.22) # Pa
  K <- Kc*(1 + Po/Ko)
  reta <- eta_fun(Tk)/eta_fun(25+273.15)
  # xi <- (beta*K / (1.6*reta))^0.5 # Pa^0.5
  ### A full model
  xi <- (beta*(K+gamma_star)/(1.6*reta))^0.5
  ###
  iWUE <- Ca / 1.6 / (xi / D^0.5 + 1)
  return(iWUE)
}

# file paths and names --------------------------------------------------
infilename1 <- '/Volumes/Land/Data_proj_wuecu/CRUNCEP_Ta_monthly_05d_post.nc' # air temperature, K
infilename2 <- '/Volumes/Land/Data_proj_wuecu/CRUNCEP_VPD_monthly_05d_post.nc' # vapor pressure deficit, kPa
infilename3 <- '/Volumes/Land/Data_proj_wuecu/mstmip_co2_monthly_05d_post.nc' # atmospheric CO2, umol/mol
infilename4 <- '/Volumes/Land/Data_proj_wuecu/GIMMS_fPAR_monthly_05d_post.nc' # fAPAR, /
infilename5 <- '/Volumes/Land/Data_proj_wuecu/CRUNCEP_Press_monthly_05d_post.nc' # Air press, Pa

outfilepath1 <- '/Volumes/Land/Data_proj_wuecu/'

# Load data ---------------------------------------------------------------
data1 <- brick(infilename1)
data2 <- brick(infilename2) # VPD need a unit transformation
data2 <- data2*1000 # kPa -> Pa
data3 <- brick(infilename3)
data4 <- brick(infilename4)
data5 <- brick(infilename5)

# Present data ------------------------------------------------------------
# In the fPAR data, the missing values were transformed into 'NA'
### test
## comment: this seems not a good way to deal with the problem ... 
## that fPAR data have no geographical information
# data3value <- values(data3)
# data3mat <- as.matrix(data3)
# data4mat <- as.matrix(data4)
# data3matmasked <- data3mat[,1]
# data3matmasked[is.na(data4mat[,1])] <- NA
# dim(data3matmasked) <- c(360, 720)
# data3masked <- raster((data3matmasked))

###
## This works after the geographical information was added
# data3test <- data3[[1:3]]
# data4test <- data4[[1:3]]
# data3test[is.na(data4test)] <- NA

# Clean data first! -------------------------------------------------------
# Criteria:
# fPAR > 0.2 (?)
# Ta > 5 K (?)
# D > 0 Pa

data1masked <- data1
data2masked <- data2
data3masked <- data3
data4masked <- data4
data5masked <- data5
# Mask invalid data
mask1 <- is.na(data1)|is.na(data2)|is.na(data3)|is.na(data4)|
  data4<0.2|data1<(5+273.15)|data2<0
data1masked[mask1] <- NA
data2masked[mask1] <- NA
data3masked[mask1] <- NA
data4masked[mask1] <- NA
data5masked[mask1] <- NA

# Model simulations ------------------------------------------------------------
# a full simulation
iwue_full <- Model_iWUE_fun1(data1masked, data2masked, data3masked, data5masked)
# gpp_globe <- Model_GPP_fun()

# simulation 1: warming or not; Ta
# create a stack
meanstack_fun <- function(n, data1) {
  # n - the number of layers
  # data1 - the raster
  if (n == 2) {
    dataout <- stack(data1, data1)
    return(dataout)
  }
  if (n > 2) {
    dataout <- stack(meanstack_fun((n-1), data1), data1)
    return(dataout)
  }
}
data1clim <- calc(data1, mean, na.rm = T) # Or mean of multiple years for each month?
data1climstack <- meanstack_fun(nlayers(data1), data1clim)
data1climstack[mask1] <- NA 
iwue_sim1_ta <- Model_iWUE_fun1(data1climstack, data2masked, data3masked, data5masked) 

# simulation 2: elevated CO2 or not
data3clim <- calc(data3, mean, na.rm = T)
data3climstack <- meanstack_fun(nlayers(data3), data3clim)
data3climstack[mask1] <- NA
iwue_sim2_co2 <- Model_iWUE_fun1(data1masked, data2masked, data3climstack, data5masked)

# simulation 3: increasing atmospheric water demand or not (D)
data2clim <- calc(data2, mean, na.rm = T)
data2climstack <- meanstack_fun(nlayers(data2), data2clim)
data2climstack[mask1] <- NA
iwue_sim3_vpd <- Model_iWUE_fun1(data1masked, data2climstack, data3masked, data5masked)

# Save the results --------------------------------------------------------

writeRaster(iwue_full, paste0(outfilepath1, 'iWUE_monthly_full_05d.nc'), 
            varname = 'iWUE_monthly', NAflag = mfill, 
            overwrite = T)

# This would be done in a separate program
# writeRaster(iwue_sim1_ta, paste0(outfilepath1, 'iWUE_monthly_ta_constant_05d.nc'), 
#             varname = 'iWUE_monthly_ta_const', NAflag = mfill, 
#             overwrite = T)
# writeRaster(iwue_sim2_co2, paste0(outfilepath1, 'iWUE_monthly_co2_constant_05d.nc'), 
#             varname = 'iWUE_monthly_co2_const', NAflag = mfill, 
#             overwrite = T)
# writeRaster(iwue_sim3_vpd, paste0(outfilepath1, 'iWUE_monthly_vpd_constant_05d.nc'), 
#             varname = 'iWUE_monthly_vpd_const', NAflag = mfill, 
#             overwrite = T)