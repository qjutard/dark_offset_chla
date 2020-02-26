#############################################################################
# Function to build a list of profile file names from the argo index, a WMO,
# and the local path to the files
#############################################################################

require(stringr)
require(stringi)

file_names <- function(index_ifremer, path_to_netcdf, WMO) {
	
    files = as.character(index_ifremer$file) #retrieve the path of each netcfd file
    ident = strsplit(files,"/") #separate the different roots of the files paths
    ident = matrix(unlist(ident), ncol=4, byrow=TRUE)
    dac = ident[,1]
    wod = ident[,2] #retrieve the WMO of all profiles as a vector
	
	name_list = paste(path_to_netcdf, files[which(wod==WMO)], sep="")
	
    all_meta = paste(path_to_netcdf, dac, "/", wod, "/", wod, "_meta.nc", sep="")
    name_meta = all_meta[which(wod==WMO)]
    
	return(list("name_list"=name_list, "name_meta"=name_meta)) 
}
