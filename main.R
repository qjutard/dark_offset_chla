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
library(wesanderson)

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
use_DMMC = as.logical(uf[5])
use_kal = as.logical(uf[6])
runmed_size = as.numeric(uf[7])

# apply default values if necessary

if (plot_name=="NA") {
    plot_name = paste("DARK_", WMO, ".png", sep="")
}
if (median_size=="NA") {
    median_size = 5
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
index_greylist = read.csv(path_to_index_greylist, sep = ",") # if greylist is useful at some point

DEEP_EST= NULL
if (use_DMMC) {
    DEEP_EST = Dark_MLD_table_coriolis(WMO, "/DATA/ftp.ifremer.fr/ifremer/argo/dac/", index_ifremer)
}

numCores = detectCores()
M = mcmapply(open_profiles, name_list, MoreArgs=list("CHLA", DEEP_EST, index_ifremer, index_greylist, WMO, use_DMMC=use_DMMC), mc.cores=numCores, USE.NAMES=FALSE)

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
greylist_axis = rep(NA, n_prof)
is_deep = rep(NA, n_prof)
all_minima_var = rep(NA, n_prof)
min_pres = rep(NA, n_prof)
for (i in seq(1, n_prof)) {
    chla = M[,i]$PARAM
    chla_QC = M[,i]$PARAM_QC
    pres = M[,i]$PRES
    chla[which((chla > 50) | (chla < - 0.1))] = NA
    chla[which(chla_QC=="4")] = NA
    chla = chla[which(!is.na(chla))]
    pres = pres[which(!is.na(chla))]
    
    chla_smoothed = runmed(chla, median_size, endrule="constant")
    all_minima[i] = min(chla_smoothed, na.rm = T)
    min_pres[i] = min(pres[which(chla_smoothed==all_minima[i])], na.rm=T)
    all_minima_var[i] = var(chla_smoothed[which(pres>=min_pres[i])], na.rm = T)
    
    is_deep[i] = (max(pres, na.rm=T)>800)
    
    offset_auto[i] = (M[,i]$DARK_CHLA - factory_dark) * M[,i]$SCALE_CHLA
    
    offset_DMMC[i] = M[,i]$DMMC_offset
    
    greylist_axis[i] = M[,i]$is_greylist
}
all_minima[which(is.infinite(all_minima))] = NA

median_axis = which(is.na(greylist_axis) & is_deep & !is.na(all_minima))

offset_3 = rep(median(all_minima[median_axis], na.rm=T), n_prof)
offset_1 = all_minima

offset_runmed = rep(NA, n_prof)
if (!is.na(runmed_size)) {
    offset_runmed[median_axis] = runmed(all_minima[median_axis], runmed_size, endrule="constant")
    offset_runmed[1] = offset_runmed[median_axis[1]]
    for(i in 1:n_prof) {
        if (is.na(offset_runmed[i])) {
            offset_runmed[i] = offset_runmed[i-1]
        }
    }
}


### Kalman filter on min for DM

kal_off = rep(NA, n_prof)
kal_var = rep(NA, n_prof)

if (use_kal) {
    # initialize
    kal_off[1:median_axis[1]] = unique(offset_3) # keep median for the first values
    kal_off[median_axis[1]] = unique(offset_3) # is NA if no valid profile on median_axis
    kal_var[median_axis[1]] = sd(all_minima[median_axis])*10
    
    for (i in 2:length(median_axis)) {
        
        # observation variance
        r = 0.01 + all_minima_var[median_axis[i]]
        # model_variance
        q = (0.001/10) * (M[,median_axis[i]]$JULD - M[,median_axis[i-1]]$JULD) # 0.01 every 10 days
        # model guess
        x = kal_off[median_axis[i-1]] # model by a constant
        # model variance
        P = kal_var[median_axis[i-1]] + q #model by a constant
        
        # innovation
        y = all_minima[median_axis[i]] - x
        
        # innovation covariance
        S = P + r
        
        # Kalman gain
        K = P / S
        
        # updated estimates
        kal_off[median_axis[i]] = x + K*y
        kal_var[median_axis[i]] = (1-K) * P 
    }
    for(i in 1:n_prof) {
        if (is.na(kal_off[i])) {
            kal_off[i] = kal_off[i-1]
        }
    }
}

### plot

plot_minima(M, WMO, median_size, offset_1, offset_3, offset_auto, offset_DMMC, kal_off, offset_runmed,
            plot_name, y_zoom, greylist_axis, runmed_size)


