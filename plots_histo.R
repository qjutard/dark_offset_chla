list_WMO = read.table("WMO_CHLA.list")$V1

WMO = list_WMO[1]
file_offsets = paste(WMO, "_offsets.txt", sep="")
table_offsets = read.table(file_offsets, header = TRUE)

ALL_profiles = table_offsets 
ALL_floats = NULL
ALL_floats[[1]] = table_offsets

for (i in 2:length(list_WMO)) {
    #print(WMO)
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
                   !(ALL_profiles$zone=="BLAC") &
                   (ALL_profiles$zone==zone)
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



