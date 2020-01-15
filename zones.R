#############################################################################
# The script defines the zones() function which associates an ocean or sea
# to a lat/lon couple
#############################################################################


zones <- function(lat, lon) {
    
    lon = ((lon+180)%%360)-180 #get lon in -180;180 format
    
    L=NULL
    
    L[[1]]=list("zone"="SO  ", "lat_min"=-Inf, "lat_max"=-37, "lon_min"=-Inf, "lon_max"=+Inf) #southern ocean
    L[[2]]=list("zone"="ARCT", "lat_min"=66, "lat_max"=+Inf, "lon_min"=-Inf, "lon_max"=+Inf) # arctic ocean
    L[[3]]=list("zone"="PACN", "lat_min"=0, "lat_max"=66, "lon_min"=121, "lon_max"=-90) #pacific north, lon_min>lon_max is acceptable and will be treated adequately
    L[[4]]=list("zone"="PACS", "lat_min"=-37, "lat_max"=0, "lon_min"=121, "lon_max"=-68) #pacific south
    L[[5]]=list("zone"="IND ", "lat_min"=-37, "lat_max"=29, "lon_min"=23, "lon_max"=121) #indian ocean
    L[[6]]=list("zone"="ATLN", "lat_min"=0, "lat_max"=53, "lon_min"=-90, "lon_max"=4) #atlantic north
    L[[7]]=list("zone"="ATLS", "lat_min"=-37, "lat_max"=0, "lon_min"=-68, "lon_max"=23) #atlantic south
    L[[8]]=list("zone"="IRMI", "lat_min"=53, "lat_max"=66, "lon_min"=-45, "lon_max"=10) #irminger sea
    L[[9]]=list("zone"="BALT", "lat_min"=53, "lat_max"=66, "lon_min"=10, "lon_max"=45) #baltic sea
    L[[10]]=list("zone"="MED ", "lat_min"=29, "lat_max"=53, "lon_min"=-4, "lon_max"=45) #mediterranean sea
    ### zones that overwrite others, must be at the end
    L[[11]]=list("zone"="BLAC", "lat_min"=40, "lat_max"=53, "lon_min"=27, "lon_max"=45) #black sea
    L[[12]]=list("zone"="LABR", "lat_min"=53, "lat_max"=80, "lon_min"=-90, "lon_max"=-45) #labrador sea
     
    zone_ret = NA
    for (test in L){

        lat_ok = (test$lat_min <= lat & lat < test$lat_max)
        
        if (test$lon_min < test$lon_max) { #typical case
            lon_ok = (test$lon_min <= lon & lon < test$lon_max)
        } else { #zones going over 180 longitude degrees
            lon_ok = (test$lon_min <= lon | lon < test$lon_max)
        }
        
        if (lat_ok & lon_ok) {
            zone_ret = test$zone
        }
        
    }
    
    return(zone_ret)
}
