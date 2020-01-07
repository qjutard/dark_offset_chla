#############################################################################
# This script defines the open_profiles() function which returns a parameter, 
# its pressure axis, its QC, and its date.
#############################################################################

require(ncdf4)
require(MASS)
require(stringr)
require(stringi)
require(oce)

open_profiles <- function(profile_name, PARAM_NAME, DEEP_EST, index_ifremer) {
	# profile_name is a full path to the file, PARAM_NAME is consistent
	# with bgc-argo denomination

    filenc = nc_open(profile_name, readunlim=FALSE, write=FALSE)
	
	### find the profile index	
    parameters = ncvar_get(filenc,"STATION_PARAMETERS")
	param_name_padded = str_pad(PARAM_NAME, 64, "right")
	id_prof = which(parameters==param_name_padded, arr.ind=TRUE)[2]
	if (is.na(id_prof))	{
		return(list("PARAM"=NA, "PRES"=NA, "PARAM_QC"=NA, "JULD"=NA, "param_units"=NA, "profile_id"=NA, "SCALE_CHLA"=NA, "DARK_CHLA"=NA,"DMMC_offset"=NA))
	}

	### get the parameter
	PARAM = ncvar_get(filenc, PARAM_NAME)
	PARAM = PARAM[,id_prof]

	### get the pressure
	PRES = ncvar_get(filenc, "PRES")	
	PRES = PRES[,id_prof]

	### get the QC
	PARAM_NAME_QC = paste(PARAM_NAME, "_QC", sep="")
	PARAM_QC = ncvar_get(filenc, PARAM_NAME_QC)
	PARAM_QC = PARAM_QC[id_prof]
	PARAM_QC = unlist(strsplit(PARAM_QC, ""))

	### get the date (JULD)
	JULD = ncvar_get(filenc, "JULD")
	JULD = JULD[id_prof]
	
	### get the units
	param_units = ncatt_get(filenc, PARAM_NAME, attname = "units")$value
	
	### get the profile index
	len = str_length(profile_name)
	profile_id = str_sub(profile_name, len-6, len-3) # extract 4 characters to '_000' or '000D'
	if (str_sub(profile_id,1,1)=="_") { # ascending profile
	    profile_id = as.numeric(str_sub(profile_id,2))
	} else {
	    profile_id = as.numeric(str_sub(profile_id,1,3)) - 0.5 # descending profiles are pur on half indices preceding the corresponding ascending profile
	}
	
	##### outputs specific to the project on the chla dark offset #####
	
	### get the calib parameters
	
	# get the chla_scale
	calib = ncvar_get(filenc, "SCIENTIFIC_CALIB_COEFFICIENT")
    id_chla = grep("CHLA", calib) # find chla index
    chla_calib = calib[id_chla] # get chla calibration
    chla_calib = unlist(strsplit(chla_calib,",")) # separate coefficients
    chla_calib_scale = chla_calib[grep("SCALE_CHLA",chla_calib)] # get the scale information
    chla_calib_scale = unlist(strsplit(chla_calib_scale,"=")) # separate the name from the number
    chla_scale = as.numeric(chla_calib_scale[2]) # get the scale coefficient as a number
    # get the chla dark
    chla_calib_dark = chla_calib[grep("DARK_CHLA",chla_calib)] # get the dark information
    chla_calib_dark = unlist(strsplit(chla_calib_dark,"=")) # separate the name from the number
    chla_dark = as.numeric(chla_calib_dark[2]) # get the dark coefficient as a number
    
    ### get the DMMC dark offset
    profile_actual = unlist(strsplit(profile_name,'/'))
    profile_actual = profile_actual[length(profile_actual)]
    profile_actual = str_sub(profile_actual,3,14)
    path_to_netcdf = "/DATA/ftp.ifremer.fr/ifremer/argo/dac/"
    L = process_file(profile_actual, index_ifremer, path_to_netcdf, DEEP_EST=DEEP_EST, accept_descent=TRUE, offset_override="dmmc")
    
    if (is.list(L)) {
        DMMC_offset = L$chl_dark_offset
    } else {
        DMMC_offset = NA
    }

    ###################################################################

	nc_close(filenc)
	
	return(list("PARAM"=PARAM, "PRES"=PRES, "PARAM_QC"=PARAM_QC, "JULD"=JULD, "param_units"=param_units, "profile_id"=profile_id, "SCALE_CHLA"=chla_scale, "DARK_CHLA"=chla_dark, "DMMC_offset"=DMMC_offset))
}

