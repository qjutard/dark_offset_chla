#############################################################################
# Function to build a list of profile file names from the argo index, a WMO,
# and the local path to the files
#############################################################################

require(stringr)
require(stringi)

file_names <- function(index_ifremer, path_to_netcdf_before_WMO, WMO, path_to_netcdf_after_WMO) {
	
	# TODO : implement
	
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
    name_list = system2("ls", prof_id, stdout=TRUE) 
	
	return(name_list)
}
