#############################################################################
# Main script to build a list of profiles, call computations of offset, 
# and create desired plots
# open_profiles.R and file_names.R adapted from the work on BGTS
# available at https://github.com/qjutard/time_series_plot
#############################################################################

library(ncdf4)
library(oce)
library(MASS)
library(stringr)
library(parallel)
library(stringi)
library(RColorBrewer)

source("~/Documents/dark_chla/dark_offset_chla/pathways.R")
source(paste(path_to_source, "file_names.R", sep=""))
source(paste(path_to_source, "open_profiles.R", sep=""))
source(paste(path_to_source, "plot_minima.R", sep=""))
source(paste(path_to_DMMC, "process_files.R", sep="")) # only necessary if -M is going to be used
source(paste(path_to_DMMC, "error_message.R", sep="")) # only necessary if -M is going to be used

### Set parameters

uf = commandArgs(trailingOnly = TRUE)

WMO = uf[1]
plot_name = uf[2]
median_size = uf[3]
y_zoom_call = uf[4]
use_DMMC = as.logical(uf[5])
use_kal = as.logical(uf[6])
runmed_size = as.numeric(uf[7])
date_axis = as.logical(uf[8])
do_write = as.logical(uf[9])

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

num_cores = detectCores()

### Build list of file names from WMO and argo_index

index_ifremer = read.table(path_to_index_ifremer, sep=",", header = T)

names = file_names(index_ifremer, path_to_netcdf, WMO)
name_list = names$name_list
name_meta = names$name_meta

### Get a list with information on all profiles
index_greylist = read.csv(path_to_index_greylist, sep = ",")

DEEP_EST= NULL
if (use_DMMC) {
    DEEP_EST = Dark_MLD_table_coriolis(WMO, path_to_netcdf, index_ifremer, n_cores=num_cores)
}

M = mcmapply(open_profiles, name_list, 
             MoreArgs=list(PARAM_NAME="CHLA", index_ifremer=index_ifremer, index_greylist=index_greylist, WMO=WMO, 
                           use_DMMC=use_DMMC, DEEP_EST=DEEP_EST, path_to_netcdf=path_to_netcdf), 
             mc.cores=num_cores, USE.NAMES=FALSE)

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

offset_med = rep(median(all_minima[median_axis], na.rm=T), n_prof)
offset_min = all_minima

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

offset_kal = rep(NA, n_prof)
kal_var = rep(NA, n_prof)

if (use_kal) {
    # initialize
    offset_kal[1:median_axis[1]] = unique(offset_med) # keep median for the first values
    offset_kal[median_axis[1]] = unique(offset_med) # is NA if no valid profile on median_axis
    kal_var[median_axis[1]] = sd(all_minima[median_axis])*10
    
    for (i in 2:length(median_axis)) {
        
        # observation variance
        r = 0.01 + all_minima_var[median_axis[i]]
        # model_variance
        q = (0.001/10) * (M[,median_axis[i]]$JULD - M[,median_axis[i-1]]$JULD) # 0.01 every 10 days
        # model guess
        x = offset_kal[median_axis[i-1]] # model by a constant
        # model variance
        P = kal_var[median_axis[i-1]] + q #model by a constant
        
        # innovation
        y = all_minima[median_axis[i]] - x
        
        # innovation covariance
        S = P + r
        
        # Kalman gain
        K = P / S
        
        # updated estimates
        offset_kal[median_axis[i]] = x + K*y
        kal_var[median_axis[i]] = (1-K) * P 
    }
    for(i in 1:n_prof) {
        if (is.na(offset_kal[i])) {
            offset_kal[i] = offset_kal[i-1]
        }
    }
}

cut_names = str_split(name_list, "/", simplify = T)
cut_names = cut_names[,length(cut_names[1,])]
cut_names = str_sub(cut_names, 3, 14)

if (do_write) {
    all_offsets = list(offset_min, offset_med, offset_kal, offset_runmed)
    write_names = c("minima", "median", "kalman", "runmed")
    write_filenames = paste("offsets_", WMO, "_", write_names, ".t", sep="")
    
    for (i in 1:length(all_offsets)) {
        if (!all(is.na(all_offsets[[i]]))) {
            write.table(list(cut_names, all_offsets[[i]]), write_filenames[i], col.names = FALSE, row.names = FALSE)
        }
    }
    
}

### plot

plot_minima(M=M, WMO=WMO, median_size=median_size, offset_min=offset_min, offset_med=offset_med, offset_auto=offset_auto, offset_DMMC=offset_DMMC, 
            offset_kal=offset_kal, offset_runmed=offset_runmed, plot_name=plot_name, y_zoom=y_zoom, greylist_axis=greylist_axis, runmed_size=runmed_size, 
            date_axis=date_axis)


