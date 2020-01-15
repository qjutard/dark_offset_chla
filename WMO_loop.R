#############################################################################
# This script performs a loop on the floats to extract offset data relevant
# for future plots, they are then written in a text file for each float
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
source(paste(path_to_source, "zones.R", sep=""))

median_size = 5
num_cores = detectCores()

index_ifremer = read.table(path_to_index_ifremer, skip=9, sep = ",")
index_greylist = read.csv(path_to_index_greylist, sep = ",") # if greylist is useful at some point

treat_WMO <- function(WMO) {
    name_list = file_names(index_ifremer, path_to_netcdf_before_WMO, WMO, path_to_netcdf_after_WMO)
    name_meta = paste(path_to_netcdf_before_WMO, WMO, "/", WMO, "_meta.nc", sep="")
    
    M = mapply(open_profiles, name_list, MoreArgs=list("CHLA", DEEP_EST=NULL, index_ifremer, index_greylist, WMO, use_DMMC=FALSE), USE.NAMES=FALSE)
    
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
    greylist_axis = rep(NA, n_prof)
    zones_axis = rep(NA, n_prof)
    is_deep = rep(NA, n_prof)
    prof_names = rep(NA, n_prof)
    negative_before_auto = rep(NA, n_prof)
    negative_after_auto = rep(NA, n_prof)
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
        if (max(pres, na.rm=T)>800) {
            is_deep[i] = TRUE
        } else { 
            is_deep[i] = FALSE
        }
        
        offset_auto[i] = (M[,i]$DARK_CHLA - factory_dark) * M[,i]$SCALE_CHLA
        
        negative_before_auto[i] = length(which(chla<0))
        negative_after_auto[i] = length(which(chla<offset_auto[i]))
        
        greylist_axis[i] = M[,i]$is_greylist
        
        lat = M[,i]$lat
        lon = M[,i]$lon
        if(!is.na(lat) & !is.na(lon)) {
            zones_axis[i] = zones(lat, lon)
        }
        
        splitted_name = unlist(strsplit(name_list[i],"/"))
        prof_names[i] = splitted_name[length(splitted_name)]
        
    }
    all_minima[which(is.infinite(all_minima))] = NA
    
    offset_med = rep(median(all_minima[which(is.na(greylist_axis) & is_deep)], na.rm=T), n_prof)
    offset_min = all_minima
    
    res=list("profile"=prof_names, "zone"=zones_axis, "off_auto"=offset_auto, "off_min"=offset_min, "off_med"=offset_med, 
             "greylist"=greylist_axis, "is_deep"=is_deep, "negative_before_auto"=negative_before_auto, "negative_after_auto"=negative_after_auto)
    
    output_name = paste(WMO, "_offsets.txt", sep="")
    
    write.table(file=output_name, res, row.names=FALSE, col.names=TRUE)
    
    return(0)
}

#list_WMO = c("6901524","6901527","6901439","6902738")
list_WMO = read.table("WMO_CHLA.list")$V1

M = mcmapply(treat_WMO, list_WMO, mc.cores=num_cores)
