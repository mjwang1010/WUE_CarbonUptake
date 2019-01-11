# GPP remote sensing model based on the optimal theory
# Input - gridded data

rm(list = ls())

# Load libraries ----------------------------------------------------------
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

Po_fun2 <- function(Press) {
  # input - press, surface pressure, Pa
  # Po, partial pressure of O2, Pa
  Po <- Press*0.21
  return(Po)
}

# Estimate the IWUE from the environmental variables
Model_iWUE_fun1 <- function(Tk, D, Ca, Press) {
  
  # Comments: the model does not consider the impact of soil moisture.
  # Tk, temperature, (K)
  # D, VPD, (Pa)
  # Ca, CO2 concentration, ppm, umol mol-1
  # elev, elevation, km (to estimate O2 partial pressure)
  # Press, atmospheric pressure, Pa (the second method to estimate O2 partial pressure)
  # iWUE, umol mol-1
  
  Po <- Po_fun2(Press) # unit - Pa; or a function of elevation
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

# Estimation GPP from remote sensing data
# Model_GPP_fun <- function(Tk, D, Ca, PAR, fPAR, elev) {
# Model_GPP_fun1 <- function(Tk, D, Ca, PAR, fPAR, Press) {
Model_GPP_fun2 <- function(Tk, D, Ca, PAR, fPAR, iWUE) {
  
  # Tk, temperature, K
  # D, vapor pressure deficit, Pa
  # Ca, CO2, umol mol-1
  # PAR, umol m-2 s-1
  # fPAR, /
  # elev, km (to estimate oxygen partial pressure; instead pressure is used)
  # Press, air pressure, Pa
  ###
  # What about impacts of water availabilities?
  
  ###
  # GPP, umol m-2 s-1
  
  phic <- 0.093 # Quantum efficiency
  # phic <- 0.06 # Or a free parameter (Stocker, sofun model)
  
  # iWUE <- Model_iWUE_fun1(Tk, D, Ca, Press)
  # Here, we used estimated "iWUE" directly.
  
  # the unit of gamma_star is dependent on x25
  # Medlyn et. al. (2002)
  # x25 = 42.75 (umol mol-1), at temperature 25 degree for gamma_star
  gamma_star <- K_fun(37830, Tk, 42.75) # umol mol-1
  APAR <- PAR * fPAR # umol m-2 s-1
  m0 <- (1 - 3*gamma_star / (Ca - 1.6*iWUE + 2*gamma_star))
  # Based on Wang's paper (Nature Plant, 2016)
  cstar <- 4*0.103 # according to the optimal Jmax
  # "m" parameter, a temporal variable
  mtemp <- m0^2 - m0^(4.0/3.0) * cstar^(2.0/3.0) 
  # adjusted prameter m
  mtemp[mtemp<0] <- NA
  # This method caused error: 'object of type 'S4' is not subsettable'
  # ma[mtemp>=0] <- (mtemp[mtemp>=0])^0.5
  ma <- mtemp^0.5
  
  GPP = phic * APAR * ma
  
  return(GPP)
}

Model_GPP_fun3 <- function(Tk, D, Ca, PAR, fPAR, iWUE, WSI) {
  
  # WSI, water stress index, - (PET/P?)
}

# file paths and names --------------------------------------------------
infilename1 <- '/Volumes/Land/Data_proj_wuecu/CRUNCEP_Ta_monthly_05d_post.nc' # air temperature, K
infilename2 <- '/Volumes/Land/Data_proj_wuecu/CRUNCEP_VPD_monthly_05d_post.nc' # vapor pressure deficit, kPa
infilename3 <- '/Volumes/Land/Data_proj_wuecu/mstmip_co2_monthly_05d_post.nc' # atmospheric CO2, umol/mol
infilename4 <- '/Volumes/Land/Data_proj_wuecu/GIMMS_fPAR_monthly_05d_post.nc' # fAPAR, /
infilename5 <- '/Volumes/Land/Data_proj_wuecu/CRUNCEP_Press_monthly_05d_post.nc' # Air press, Pa
infilename6 <- '/Volumes/Land/Data_proj_wuecu/CRUNCEP_SWdown_monthly_05d_post.nc' # Shortwave radiation, W m-2
infilename7 <- '/Volumes/Land/Data_proj_wuecu/iWUE_monthly_full_05d.nc' # iWUE with full inputs, umol mol-1

outfilepath1 <- '/Volumes/Land/Data_proj_wuecu/'

# Load data ---------------------------------------------------------------
data1 <- brick(infilename1)
data2 <- brick(infilename2) # VPD need a unit transformation
data2 <- data2*1000 # kPa -> Pa
data3 <- brick(infilename3)
data4 <- brick(infilename4)
data5 <- brick(infilename5)
data6 <- brick(infilename6)
data7 <- brick(infilename7)

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
# Radiation?

data1masked <- data1
data2masked <- data2
data3masked <- data3
data4masked <- data4
data5masked <- data5
data6masked <- data6
data7masked <- data7

# Mask invalid data
mask1 <- is.na(data1)|is.na(data2)|is.na(data3)|is.na(data4)|
  data4<0.2|data1<(5+273.15)|data2<100

data1masked[mask1] <- NA
data2masked[mask1] <- NA
data3masked[mask1] <- NA
data4masked[mask1] <- NA
data5masked[mask1] <- NA
data6masked[mask1] <- NA
data7masked[mask1] <- NA

fSR <- 4.6 # From W m-2 to umol m-2 s-1
data8masked <- data6masked*fSR*0.45 # from shortwave radiation to PAR, umol m-2 s-1

# Estimate global GPP -----------------------------------------------------

GPP_full_nows <- Model_GPP_fun2(data1masked, data2masked, data3masked, data8masked, data4masked, data7masked)

# Save the data -----------------------------------------------------------
writeRaster(GPP_full_nows, paste0(outfilepath1, 'GPP_monthly_full_nows_05d.nc'), 
            varname = 'GPP_monthly', NAflag = mfill, 
            overwrite = T)

