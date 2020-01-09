#############################################################################
# Main script to build a list of profiles, call computations of offset, 
# and create desired plots
# open_profiles.R and file_names.R are taken directly from the work on BGTS
# available at https://github.com/qjutard/time_series_plot
#############################################################################

library(ncdf4)
library(oce)
library(MASS)
library(stringr)
library(parallel)
library(stringi)

source("~/Documents/dark_chla/dark_offset_chla/pathways.R")
source(paste(path_to_source, "file_names.R", sep=""))
source(paste(path_to_source, "open_profiles.R", sep=""))
source(paste(path_to_source, "plot_minima.R", sep=""))
source("~/Documents/cornec_chla_qc/chl_bbp_ttt/process_files.R")
source("~/Documents/cornec_chla_qc/chl_bbp_ttt/error_message.R")

### Set parameters

uf = commandArgs(trailingOnly = TRUE)

WMO = uf[1]
plot_name = uf[2]
median_size = uf[3]
y_zoom_call = uf[4]

# apply default values if necessary

if (plot_name=="NA") {
    plot_name = paste("DARK_", WMO, ".png", sep="")
}
if (median_size=="NA") {
    median_size = 1
}
if (y_zoom_call=="NA") {
	y_zoom = NULL
} else {
	y_zoom = as.numeric(unlist(strsplit(y_zoom_call, ";")))
}


#WMO = "6901524"
#median_size = 5

### Build list of file names from WMO and argo_index

index_ifremer = read.table(path_to_index_ifremer, skip=9, sep = ",")

name_list = file_names(index_ifremer, path_to_netcdf_before_WMO, WMO, path_to_netcdf_after_WMO)
name_meta = paste(path_to_netcdf_before_WMO, WMO, "/", WMO, "_meta.nc", sep="")

### Get a list with information on all profiles
#index_greylist = read.csv(path_to_index_greylist, sep = ",") # if greylist is useful at some point

DEEP_EST = Dark_MLD_table_coriolis(WMO, "/mnt/c/DATA/ftp.ifremer.fr/ifremer/argo/dac/", index_ifremer)

numCores = detectCores()
M = mcmapply(open_profiles, name_list, MoreArgs=list("CHLA", DEEP_EST, index_ifremer), mc.cores=numCores, USE.NAMES=FALSE)

### compute minima

# get the factory dark
metanc = nc_open(name_meta)
calib = ncvar_get(metanc, "PREDEPLOYMENT_CALIB_COEFFICIENT")
id_chla = grep("CHLA", calib) # find chla index
chla_calib = calib[id_chla] # get chla calibration
chla_calib = unlist(strsplit(chla_calib,",")) # separate coefficients
chla_calib_dark = chla_calib[grep("DARK_CHLA",chla_calib)] # get the dark information
chla_calib_dark = unlist(strsplit(chla_calib_dark,"=")) # separate the name from the number
factory_dark = as.numeric(chla_calib_dark[2]) # get the dark coefficient as a number
nc_close(metanc)

n_prof = dim(M)[2]
all_minima = rep(NA, n_prof)
offset_auto = rep(NA, n_prof)
offset_DMMC = rep(NA, n_prof)
for (i in seq(1, n_prof)) {
    chla = M[,i]$PARAM
    chla[which((chla > 50) | (chla < - 0.1))]=NA
    chla = chla[which(!is.na(chla))]
    
    chla_smoothed = runmed(chla, median_size, endrule="constant")
    all_minima[i] = min(chla_smoothed, na.rm = T)
    
    offset_auto[i] = (M[,i]$DARK_CHLA - factory_dark) * M[,i]$SCALE_CHLA
    
    offset_DMMC[i] = M[,i]$DMMC_offset
}
all_minima[which(is.infinite(all_minima))] = NA

offset_1 = all_minima
offset_3 = rep(median(all_minima, na.rm=T), n_prof)

### plot

plot_minima(M, WMO, median_size, offset_1, offset_3, offset_auto, offset_DMMC, plot_name, y_zoom)


