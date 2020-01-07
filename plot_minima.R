#############################################################################
# Defines a function that plots minima of profiles from a float to
# compare methods 1 and 3 from the ADMT20 presentation by Xing
#############################################################################

plot_minima <- function(M, WMO, median_size, offset_1, offset_3, offset_auto, offset_DMMC, plot_name) {
    n_prof = dim(M)[2]
    
    juld = rep(NA, n_prof)
    percent_qc4 = rep(NA, n_prof)
    for (i in 1:n_prof) {
        juld[i] = M[,i]$JULD
        profile_qc = M[,i]$PARAM_QC
        profile_qc = profile_qc[which(!is.na(profile_qc)&profile_qc!=" ")]
        if (length(profile_qc)>0) {
            percent_qc4[i] = length(which(profile_qc=="4")) / length(profile_qc)
        }
    }
    dates = as.Date(juld, origin='1950-01-01')
    QC_colors = rep(rgb(1,1,1),n_prof)
    QC_colors[which(!is.na(percent_qc4))] = rgb(percent_qc4[which(!is.na(percent_qc4))], 1-percent_qc4[which(!is.na(percent_qc4))], rep(0,length(which(!is.na(percent_qc4)))))
    
    
    png(plot_name, width = 800, height = 400)
    
    Xrange = range(as.numeric(dates), na.rm=T)
    Yrange = range(c(offset_1, offset_auto, offset_DMMC), na.rm=T)
    Yrange[2] = Yrange[2]+(Yrange[2]-Yrange[1])*0.4
    
    plot(dates, offset_1, xlab = "time", ylab="chla offset",xlim=Xrange, ylim=Yrange)
    title(main=paste("Visualisation of the median of minima for the \ncomputation of the dark offset of",WMO), sub=paste("median size =",median_size))
    points(dates, offset_auto, pch="x", col="green")
    points(dates, offset_DMMC, pch="+", col="blue")
    points(dates, rep(Yrange[2]+(Yrange[2]-Yrange[1])*0.02, n_prof), col=QC_colors, pch=15)
    lines(dates, offset_3, col="red")
    legend(x=Xrange[2]*0.8+Xrange[1]*0.2, y=Yrange[2], legend=c("profile minium","median of minima","automatic offset","DMMC offset"), pch=c('o','_','x','+'), col=c('black','red','green','blue'))
    
    dev.off()
    
    return(0)
}