#############################################################################
# Defines a function that plots minima of profiles from a float to
# compare methods 1 and 3 from the ADMT20 presentation by Xing
#############################################################################

plot_minima <- function(M, WMO, median_size, offset_1, offset_3, offset_auto, offset_DMMC, plot_name) {
    n_prof = dim(M)[2]
    
    juld = rep(NA, n_prof)
    for (i in 1:n_prof) {
        juld[i] = M[,i]$JULD
    }
    dates = as.Date(juld, origin='1950-01-01')
    
    png(plot_name, width = 800, height = 400)
    
    
    
    plot(dates, offset_1, xlab = "time", ylab="chla offset", ylim=range(c(offset_1, offset_auto, offset_DMMC),na.rm=T))
    title(main=paste("Visualisation of the median of minima for the \ncomputation of the dark offset of",WMO), sub=paste("median size =",median_size))
    points(dates, offset_auto, pch="x", col="green")
    points(dates, offset_DMMC, pch="+", col="blue")
    lines(dates, offset_3, col="red")
    
    dev.off()
    
    return(0)
}