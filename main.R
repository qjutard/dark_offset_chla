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

### Set parameters

WMO = "6901524"

### Build list of file names from WMO and argo_index

index_ifremer = read.table(path_to_index_ifremer, skip=9, sep = ",")

name_list = file_names(index_ifremer, path_to_netcdf_before_WMO, WMO, path_to_netcdf_after_WMO)

### Get a list with information on all profiles
#index_greylist = read.csv(path_to_index_greylist, sep = ",") # if greylist is useful at some point

numCores = detectCores()
M = mcmapply(open_profiles, name_list, MoreArgs=list("CHLA"), mc.cores=numCores, USE.NAMES=FALSE)



### plot


