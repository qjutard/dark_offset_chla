#############################################################################
# This script defines the open_profiles() function which returns a parameter, 
# its pressure axis, its QC, and its date.
#############################################################################

require(ncdf4)
require(MASS)
require(stringr)
require(stringi)
require(oce)

open_profiles <- function(profile_name, PARAM_NAME, index_ifremer, index_greylist, WMO, use_DMMC=FALSE, DEEP_EST=NULL, path_to_netcdf=NULL) {
	# profile_name is a full path to the file, PARAM_NAME is consistent
	# with bgc-argo denomination
    #print(profile_name)
    filenc = nc_open(profile_name, readunlim=FALSE, write=FALSE)
	
	### find the profile index	
    parameters = ncvar_get(filenc,"STATION_PARAMETERS")
	param_name_padded = str_pad(PARAM_NAME, 64, "right")
	id_prof_arr = which(parameters==param_name_padded, arr.ind=TRUE)
	if (length(id_prof_arr)==0) {
	    return(list("PARAM"=NA, "PRES"=NA, "PARAM_QC"=NA, "JULD"=NA, "param_units"=NA, "profile_id"=NA, "SCALE_CHLA"=NA, 
	                "DARK_CHLA"=NA,"DMMC_offset"=NA, "is_greylist"=NA, "lat"=NA, "lon"=NA))
	}
	if (length(id_prof_arr)==2) { #if id_prof_arr is a vector, there is only one PARAM_NAME profile
	    id_prof = id_prof_arr[2]
	} else {
	    id_prof = id_prof_arr[1,2] #take the first profile
	    print(paste("Several profiles of", PARAM_NAME,"detected, only using the first one"))
	}
	

	### get the parameter
	PARAM = ncvar_get(filenc, PARAM_NAME, start=c(1,id_prof), count=c(-1,1))

	### get the pressure
	PRES = ncvar_get(filenc, "PRES", start=c(1,id_prof), count=c(-1,1))	

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
	
	lat<- ncvar_get(filenc,"LATITUDE")[1]
	lon<- ncvar_get(filenc,"LONGITUDE")[1]
	
	position_qc<-substr(ncvar_get(filenc,"POSITION_QC"),1,1) # read position QC
	
	# skip the profile if the position QC is bad
	if (position_qc == 3 | position_qc==4) {
	    lat = NA
        lon = NA
	}
	
	### get the calib parameters
	
	# get the chla_scale
	calib = ncvar_get(filenc, "SCIENTIFIC_CALIB_COEFFICIENT")
    id_chla = grep("CHLA", calib) # find chla index
    chla_scale = NA
    chla_dark = NA
    if (length(id_chla)>0) {
        chla_calib = calib[id_chla] # get chla calibration
        chla_calib = unlist(strsplit(chla_calib,",")) # separate coefficients
        chla_calib_scale = chla_calib[grep("SCALE_CHLA",chla_calib)] # get the scale information
        chla_calib_scale = unlist(strsplit(chla_calib_scale,"=")) # separate the name from the number
        chla_scale = as.numeric(chla_calib_scale[2]) # get the scale coefficient as a number
        # get the chla dark
        chla_calib_dark = chla_calib[grep("DARK_CHLA",chla_calib)] # get the dark information
        chla_calib_dark = unlist(strsplit(chla_calib_dark,"=")) # separate the name from the number
        chla_dark = as.numeric(chla_calib_dark[2]) # get the dark coefficient as a number
    }
    
    ### get the DMMC dark offset and greylist
    DMMC_offset = NA
    if (use_DMMC) {
        profile_actual = unlist(strsplit(profile_name,'/'))
        profile_actual = profile_actual[length(profile_actual)]
        profile_actual = str_sub(profile_actual,3,14)
        L = process_file(profile_actual, index_ifremer, path_to_netcdf, DEEP_EST=DEEP_EST, accept_descent=TRUE, 
                         offset_override="dmmc", index_greylist=index_greylist)
        if (is.list(L)) {
            DMMC_offset = L$chl_dark_offset
            if (is.na(DMMC_offset)) {DMMC_offset=0}
        } 
    }
    is_greylist = NA
    if (!is.null(index_greylist) & !is.na(JULD)) {
        
        indices_greylist = which( index_greylist$PLATFORM_CODE==WMO & (index_greylist$PARAMETER_NAME=="CHLA" | index_greylist$PARAMETER_NAME=="BBP700") )
        
        prof_date_trunc = stri_datetime_format(as.Date(JULD, origin='1950-01-01'), format="uuuuMMdd")
        
        for (j in indices_greylist) {
            
            ## is the profile on the greylist ?
            prof_in_greylist = FALSE
            if (is.na(index_greylist$END_DATE[j])) { # all past that date
                if (prof_date_trunc>=index_greylist$START_DATE[j]) {
                    prof_in_greylist = TRUE
                } 
            } else { # date interval
                if (index_greylist$START_DATE[j]<=prof_date_trunc & prof_date_trunc<=index_greylist$END_DATE[j]) {
                    prof_in_greylist = TRUE
                }
            }
            
            ## what is the QC and what to do ?
            if (prof_in_greylist){
                if (index_greylist$QUALITY_CODE[j] == 4) {
                    #print(paste("profile on the greylist with QC 4 at index ", j, " with comment : ", index_greylist$COMMENT[j], sep=""))
                    is_greylist = 4
                } else if (index_greylist$QUALITY_CODE[j] == 3) {
                    is_greylist = 3
                }
            }     
        }
    }
    ###################################################################

	nc_close(filenc)
	
	return(list("PARAM"=PARAM, "PRES"=PRES, "PARAM_QC"=PARAM_QC, "JULD"=JULD, "param_units"=param_units, "profile_id"=profile_id, 
	            "SCALE_CHLA"=chla_scale, "DARK_CHLA"=chla_dark, "DMMC_offset"=DMMC_offset, "is_greylist"=is_greylist, "lat"=lat, "lon"=lon))
}

