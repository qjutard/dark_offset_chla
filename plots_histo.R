library(stringr)
library(stringi)
library(ggplot2)
library(dplyr)
library(maps)


list_WMO = read.table("WMO_CHLA.list")$V1
list_WMO_incois = read.table("incois_WMO_CHLA.list")$V1
list_WMO_bodc = read.table("bodc_WMO_CHLA.list")$V1

list_WMO = c(list_WMO, list_WMO_incois, list_WMO_bodc)

WMO = list_WMO[1]
# file_offsets = paste(WMO, "_offsets.txt", sep="")
table_offsets = read.table(file_offsets, header = TRUE)

ALL_profiles = table_offsets 
ALL_floats = NULL
ALL_floats[[1]] = table_offsets

for (i in 2:length(list_WMO)) {

    WMO = list_WMO[i]
    
    file_offsets = paste(WMO, "_offsets.txt", sep="")
    table_offsets = read.table(file_offsets, header = TRUE)
    
    ALL_profiles = merge(ALL_profiles, table_offsets, all=T)
    ALL_floats[[i]] = table_offsets
}

zone = "IND "

sample = which(is.na(ALL_profiles$greylist) & 
                   ALL_profiles$is_deep & 
                   !(str_sub(ALL_profiles$profile, 3, 9)=="3901066") & 
                   !(str_sub(ALL_profiles$profile, 3, 9)=="3901067") &
                   !is.na(ALL_profiles$off_auto) &
                   !(ALL_profiles$zone=="BLAC")
               )

created_negatives = ALL_profiles$negative_after_auto - ALL_profiles$negative_before_auto

auto_to_med = ALL_profiles$off_auto - ALL_profiles$off_med


hist(created_negatives[sample], plot = T, breaks = 50, col="slateblue1")
hist(auto_to_med[sample], plot = T, breaks = 50, col="slateblue1")
hist(ALL_profiles$off_auto[sample], plot=T, breaks = 50, col="slateblue1")
hist(ALL_profiles$off_min[sample], plot=T, breaks = 50, col="slateblue1")
hist(ALL_profiles$off_med[sample], plot=T, breaks = 50, col="slateblue1")

mirrored_hist(ALL_profiles$off_auto[sample], ALL_profiles$off_min[sample], 50)
mirrored_hist(ALL_profiles$off_auto[sample], ALL_profiles$off_med[sample], 50)

all_zones = c("ATLN", "ATLS", "IND ", "SO ", "LABR", "MED ", "BLAC", "ARCT", "BALT", "IRMI", "PACN", "PACS")

for (zone in all_zones) {
    plot_zone(ALL_profiles, zone)
}

# tests on negative auto_to_med

auto_to_med_sample = auto_to_med[sample]
off_auto_sample = ALL_profiles$off_auto[sample]
off_med_sample = ALL_profiles$off_med[sample]
prof_name_sample = ALL_profiles$profile[sample]
zone_sample = ALL_profiles$zone[sample]

cutoff = -0.02

mirrored_hist(off_auto_sample[which(auto_to_med_sample<cutoff)], off_med_sample[which(auto_to_med_sample<cutoff)], 70)

negative_diff_floats = unique(str_sub(prof_name_sample[which(auto_to_med_sample<cutoff & off_auto_sample==0)],3,9))

# map plot
long=ALL_profiles$lon[sample]
lat=ALL_profiles$lat[sample]
plot_profiles = data.frame(diff=auto_to_med_sample, long=ALL_profiles$lon[sample], lat=ALL_profiles$lat[sample])
cutoff=-0.02
plot_profiles2 = data.frame(diff=auto_to_med_sample[which(auto_to_med_sample>=cutoff)], long=long[which(auto_to_med_sample>=cutoff)], lat=lat[which(auto_to_med_sample>=cutoff)])
plot_profiles2$diff[which(plot_profiles2$diff<0)] = 0
plot_profiles2$diff[which(plot_profiles2$diff>0.07)] = 0.07


map <- ggplot(world_map, aes(x = long, y = lat, group=group)) +
    geom_polygon(fill="darkgray") + theme_dark() +
    scale_x_continuous(limits = c(-180,-70)) + 
    scale_y_continuous(limits = c(-60,60))
    
map + geom_point(data = plot_profiles2, aes(x = long, y = lat, group=NULL, color=diff), size=0.1) +
    #scale_color_gradient2(midpoint=0, low="blue", mid="white", high="red", space ="Lab" )
    scale_color_gradientn(colors=viridisLite::plasma(10))

mirrored_hist <- function(x1, x2, br) {
    par(mfrow=c(2,1))
  
    hist1 = hist(x1, plot=F, breaks = br)
    hist2 = hist(x2, plot=F, breaks = br)
    
    Xrange = range(c(x1, x2))
    Yrange = range(c(hist1$counts, hist2$counts))
    
    #Make the plot
    par(mar=c(0,5,3,3))
    hist(x1 , main="" , xlim=Xrange, xlab="", ylim=Yrange , xaxt="n", las=1 , col="slateblue1", breaks=br)
    par(mar=c(5,5,0,3))
    hist(x2 , main="" , xlim=Xrange, xlab="Value of my variable", ylim=rev(Yrange) , las=1 , col="tomato3"  , breaks=br)
    
}

plot_zone <- function(ALL_profiles, zone) {
    sample = which(is.na(ALL_profiles$greylist) & 
                       ALL_profiles$is_deep & 
                       !(str_sub(ALL_profiles$profile, 3, 9)=="3901066") & 
                       !(str_sub(ALL_profiles$profile, 3, 9)=="3901067") &
                       !is.na(ALL_profiles$off_auto) &
                       (ALL_profiles$zone==zone)
    )
    
    auto_to_med = ALL_profiles$off_auto - ALL_profiles$off_med
    
    br=70
    
    x1 = ALL_profiles$off_auto[sample]
    x2 = ALL_profiles$off_med[sample]
    
    hist1 = hist(x1, plot=F, breaks = br)
    hist2 = hist(x2, plot=F, breaks = br)
    
    Xrange = range(c(x1, x2))
    Yrange = range(c(hist1$counts, hist2$counts))
    
    plot_name = paste("histograms_",zone,".png", sep="")
    
    png(plot_name, width = 600, height = 600)
    
    par(mfrow=c(3,1))
    
    #Make the plot
    par(mar=c(0,5,3,3))
    hist(ALL_profiles$off_auto[sample] , main=paste("Histograms in zone :",zone) , xlim=Xrange, xlab="", ylim=Yrange , xaxt="n", las=1 , col="slateblue1", breaks=br, ylab="Automatic")
    par(mar=c(5,5,0,3))
    hist(ALL_profiles$off_med[sample] , main="" , xlim=Xrange, xlab="mg/m³", ylim=rev(Yrange) , las=1 , col="tomato3"  , breaks=br, ylab="Median")
    par()
    hist(auto_to_med[sample], breaks = br, col="slateblue1", main="", xlab="mg/m³", ylab="Auto - Med")
    
    dev.off()
    
}

dev.off()



