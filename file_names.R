#############################################################################
# Function to build a list of profile file names from the argo index, a WMO,
# and the local path to the files
#############################################################################

require(stringr)
require(stringi)

file_names <- function(index_ifremer, path_to_netcdf_before_WMO, WMO, path_to_netcdf_after_WMO) {
	
    files<-as.character(index_ifremer[,1]) #retrieve the path of each netcfd file
    ident<-strsplit(files,"/") #separate the different roots of the files paths
    ident<-matrix(unlist(ident), ncol=4, byrow=TRUE)
    prof_id<-ident[,4] #retrieve all profiles  name as a vector

	# restrict the list to the desired WMO
    prof_id_WMO = substr(prof_id, 3, 9)
	prof_id = prof_id[which(prof_id_WMO==WMO)]
	
	# build file names with '?'
	prof_id = str_sub(prof_id, 3) # remove 'MD'
	prof_id = paste(path_to_netcdf_before_WMO, WMO, path_to_netcdf_after_WMO, "B?", prof_id, sep="")
	
	# identify R or D file
	name_list = rep(NA, length(prof_id))
	for (i in 1:length(prof_id)) {
        ls_match = system2("ls", prof_id[i], stdout=TRUE) 
        if (length(ls_match) == 2) { # if both R and D files exist
            name_list[i] = ls_match[1] # use the D file which is first in alphabetical order
        } else {
            name_list[i] = ls_match
        }
	}
    
	return(name_list)
}
