#!/bin/bash

usage() { 
	echo "Usage: $0 -W <WMO_number> [-n <plot_name>] [-m <median_size>] [-h]
Do '$0 -h' for help" 1>&2
	exit 1 
}
helprint() {
	echo "
#########################################################################################

DARK makes analytics plots to compare methods for the computation of the dark offset
of chla in BGV-ARGO

Usage: $0 -W <WMO_number> [-n <plot_name>] [-m <median_size>] [-h]

### Options

-W <WMO_number> : 7 digits WMO number of the float to consider.
[-n <plot_name>] : Specify a file name for the output (with pathway), if not specified
                   the default is 'DARK_WMO.png' where WMO is replaced by the 7 digit WMO
                   number. Please use a '.png' extension in your file name.
[-m <median_size>] : Specify a size for median running filters, default is 1 which
                     corresponds to no filtering
[-h] : help

#########################################################################################
" 1>&2
	exit 0
}

WMO=NA
plot_name=NA
median_size=NA

while getopts W:n:m:h option
do
case "${option}"
in
W) WMO=${OPTARG};;
n) plot_name=${OPTARG};;
m) median_size=${OPTARG};;
h) helprint;;
*) usage;;
esac
done


Rscript ~/Documents/dark_chla/dark_offset_chla/main.R $WMO $plot_name $median_size