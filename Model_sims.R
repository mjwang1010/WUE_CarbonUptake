# different simulations


rm(list = ls())

# Load libraries ----------------------------------------------------------
library(raster)
library(ncdf4)

# Set parameters ----------------------------------------------------------
R <- 8.314 # gas constant, J mol-1 K-1
# beta <- 356.51 # a constant
# beta <- 244.03 # the constant derived from Beni Stocker
mfill <- -9999.0 # the filling value
fSR <- 4.6 # From W m-2 to umol m-2 s-1

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
Model_GPP_fun2 <- function(Tk, Ca, PAR, fPAR, iWUE) {
  
  # Tk, temperature, K
  # D, vapor pressure deficit, Pa
  # Ca, CO2, umol mol-1
  # PAR, umol m-2 s-1
  # fPAR, /
  # elev, km (to estimate oxygen partial pressure; instead pressure is used)
  # Press, air pressure, Pa
  ###
  # Here, we used estimated "iWUE" directly
  # What about impacts of water availabilities?
  
  ###
  # GPP, umol m-2 s-1
  
  phic <- 0.093 # Quantum efficiency
  # phic <- 0.06 # Or a free parameter (Stocker, sofun model)
  
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

# input without iWUE
Model_GPP_fun1 <- function(Tk, D, Ca, PAR, fPAR, Press) {
  
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
  
  iWUE <- Model_iWUE_fun1(Tk, D, Ca, Press)
  
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
ta_var_fname <- 'CRUNCEP_Ta_monthly_05d_post.nc'
ta_cons_fname <- 'Ta_clim_monthly_05d.nc'
vpd_var_fname <- 'CRUNCEP_VPD_monthly_05d_post.nc'
vpd_cons_fname <- 'VPD_clim_monthly_05d.nc'
co2_var_fname <- 'mstmip_co2_monthly_05d_post.nc'
co2_cons_fname <- 'CO2_clim_monthly_05d.nc'
fpar_var_fname <- 'GIMMS_fPAR_monthly_05d_post.nc'
fpar_cons_fname <- 'fPAR_clim_monthly_05d.nc'
press_fname <- 'CRUNCEP_Press_monthly_05d_post.nc'
rad_fname <- 'CRUNCEP_SWdown_monthly_05d_post.nc'

infilepath1 <- '/Volumes/Land/Data_proj_wuecu/'
outfilepath1 <- infilepath1

# infilename7 <- '/Volumes/Land/Data_proj_wuecu/iWUE_monthly_full_05d.nc' # iWUE with full inputs, umol mol-1

# Read data ---------------------------------------------------------------

ta_cons_data <- brick(paste0(infilepath1, ta_cons_fname))
vpd_cons_data <- brick(paste0(infilepath1, vpd_cons_fname)) # KPa
vpdpa_cons_data <- vpd_cons_data*1000 # KPa to Pa
co2_cons_data <- brick(paste0(infilepath1, co2_cons_fname))
fpar_cons_data <- brick(paste0(infilepath1, fpar_cons_fname))

press_data <- brick(paste0(infilepath1, press_fname))
rad_data <- brick(paste0(infilepath1, rad_fname))
ta_var_data <- brick(paste0(infilepath1, ta_var_fname))
vpd_var_data <- brick(paste0(infilepath1, vpd_var_fname))
vpdpa_var_data <- vpd_var_data*1000
co2_var_data <- brick(paste0(infilepath1, co2_var_fname))
fpar_var_data <- brick(paste0(infilepath1, fpar_var_fname))

# Simulations -------------------------------------------------------------

# preprocess the data
# PAR
par_data <- rad_data*fSR*0.45 # from shortwave radiation to PAR, umol m-2 s-1
###
# Mask invalid data
# NA data
# FAPAR < 0.2
# Ta < 278.15 K
# VPD < 100 Pa

time1 <- Sys.time()
mask1 <- is.na(ta_var_data)|is.na(vpdpa_var_data)|is.na(fpar_var_data)|
  fpar_var_data<0.2|ta_var_data<(5+273.15)|vpdpa_var_data<100
time2 <- Sys.time()

# ta_cons_data[mask1] <- NA
# ta_var_data[mask1] <- NA
# vpdpa_cons_data[mask1] <- NA
# vpdpa_var_data[mask1] <- NA
# co2_cons_data[mask1] <- NA
# co2_var_data[mask1] <- NA
# fpar_cons_data[mask1] <- NA
# fpar_var_data[mask1] <- NA
# press_data[mask1] <- NA
# par_data[mask1] <- NA

# Baseline: constant Ta, VPD, CO2, and FAPAR

iwue_base <- Model_iWUE_fun1(ta_cons_data, vpdpa_cons_data, co2_cons_data, press_data)
gpp_base <- Model_GPP_fun2(ta_cons_data, co2_cons_data, par_data, fpar_cons_data, iwue_base)

# iWUE Sim1: Ta various, other variables constant
time1 <- Sys.time()
iwue_ta_var <- Model_iWUE_fun1(ta_var_data, vpdpa_cons_data, co2_cons_data, press_data)
# iWUE Sim2: Co2 various, other variables constant
iwue_co2_var <- Model_iWUE_fun1(ta_cons_data, vpdpa_cons_data, co2_var_data, press_data)
# iWUE Sim3: VPD various, other variables constant
iwue_vpd_var <- Model_iWUE_fun1(ta_cons_data, vpdpa_var_data, co2_cons_data, press_data)
time2 <- Sys.time()

# GPP Sim1: iWUE various, other variables constant
iwue_full <- Model_iWUE_fun1(ta_var_data, vpdpa_var_data, co2_var_data, press_data)
gpp_iwue_var <- Model_GPP_fun2(ta_cons_data, co2_cons_data, par_data, fpar_cons_data, iwue_full)

# GPP Sim2: FAPAR various, other variables constant
gpp_fpar_var <- Model_GPP_fun2(ta_cons_data, co2_cons_data, par_data, fpar_var_data, iwue_base)

# GPP Sim3: CO2 various, other variables constant
gpp_co2_var <- Model_GPP_fun2(ta_cons_data, co2_var_data, par_data, fpar_cons_data, iwue_base)

# GPP Sim4: Ta various, other variables constant
gpp_ta_var <- Model_GPP_fun2(ta_var_data, co2_cons_data, par_data, fpar_cons_data, iwue_base)




# Save the data -----------------------------------------------------------
# Mask the data
iwue_base[mask1] <- NA
gpp_base[mask1] <- NA
gpp_iwue_var[mask1] <- NA
gpp_fpar_var[mask1] <- NA
gpp_co2_var[mask1] <- NA
gpp_ta_var[mask1] <- NA
iwue_ta_var[mask1] <- NA
iwue_vpd_var[mask1] <- NA
iwue_co2_var[mask1] <- NA


# Calculate the time
begint <- system.time()
writeRaster(iwue_base, paste0(outfilepath1, 'iWUE_monthly_const_base_05d.nc'), 
            varname = 'iWUE_monthly_base', overwrite = T, NAflag = mfill)
endt <- system.time()
writeRaster(gpp_base, paste0(outfilepath1, 'GPP_monthly_const_base_05d.nc'), 
            varname = 'GPP_monthly_base', overwrite = T, NAflag = mfill)

writeRaster(gpp_co2_var, paste0(outfilepath1, 'GPP_monthly_co2_var_05d.nc'), 
            varname = 'GPP_monthly_co2_var', overwrite = T, NAflag = mfill)
writeRaster(gpp_fpar_var, paste0(outfilepath1, 'GPP_monthly_fpar_var_05d.nc'), 
            varname = 'GPP_monthly_fpar_var', overwrite = T, NAflag = mfill)
writeRaster(gpp_iwue_var, paste0(outfilepath1, 'GPP_monthly_iwue_var_05d.nc'), 
            varname = 'GPP_monthly_iwue_var', overwrite = T, NAflag = mfill)
writeRaster(gpp_ta_var, paste0(outfilepath1, 'GPP_monthly_ta_var_05d.nc'), 
            varname = 'GPP_monthly_ta_var', overwrite = T, NAflag = mfill)

writeRaster(iwue_ta_var, paste0(outfilepath1, 'iWUE_monthly_ta_var_05d.nc'), 
            varname = 'iWUE_monthly_ta_var', overwrite = T, NAflag = mfill)
writeRaster(iwue_co2_var, paste0(outfilepath1, 'iWUE_monthly_co2_var_05d.nc'), 
            varname = 'iWUE_monthly_co2_var', overwrite = T, NAflag = mfill)
writeRaster(iwue_vpd_var, paste0(outfilepath1, 'iWUE_monthly_vpd_var_05d.nc'), 
            varname = 'iWUE_monthly_vpd_var', overwrite = T, NAflag = mfill)