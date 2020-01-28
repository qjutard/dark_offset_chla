#############################################################################
# Defines a function that plots minima of profiles from a float to
# compare methods 1 and 3 from the ADMT20 presentation by Xing
#############################################################################
require(wesanderson)
library(RColorBrewer)

plot_minima <- function(M, WMO, median_size, offset_1, offset_3, offset_auto, offset_DMMC, offset_kal, plot_name, y_zoom, greylist_axis) {
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
    QC_colors[which(greylist_axis==3)] = rgb(0.5,0.5,0.5)
    QC_colors[which(greylist_axis==4)] = rgb(0,0,0)
    
    png(plot_name, width = 800, height = 400)
    
    Xrange = range(as.numeric(dates), na.rm=T)
	if (is.null(y_zoom)) {
    	Yrange = range(c(offset_1, offset_auto, offset_DMMC), na.rm=T)
    	Yrange[2] = Yrange[2]+(Yrange[2]-Yrange[1])*0.4
    } else {
		Yrange = y_zoom
    }
    
    #my_colors = wes_palette("Darjeeling1", 3)
    #col_min = my_colors[3]
    #col_med = "#000000"
    #col_auto = my_colors[2]
    #col_DMMC = my_colors[1]
    my_colors = brewer.pal(5, "Set1")
    col_min = my_colors[5]
    col_med = my_colors[1]
    col_auto = my_colors[3]
    col_DMMC = my_colors[2]
    col_kal = my_colors[4]
    
    
    plot(dates, offset_1, xlab = "time", ylab="chla offset",xlim=Xrange, ylim=Yrange, col=col_min)
    title(main=paste("Visualisation of the different methods for the computation of the dark offset of",WMO), 
          sub=paste("median size =",median_size))
    points(dates, offset_auto, pch="x", col=col_auto)
    points(dates, offset_DMMC, pch="+", col=col_DMMC)
    points(dates, offset_kal, pch="+", col=col_kal)
    points(dates, rep(Yrange[2]+(Yrange[2]-Yrange[1])*0.02, n_prof), col=QC_colors, pch=15)
    lines(dates, offset_3, col=col_med)
    
    leg_text = c("profile minimum","median of minima","automatic offset")
    leg_symb = c('o','_','x')
    leg_col = c(col_min, col_med, col_auto)
    if (!all(is.na(offset_DMMC))) {
        leg_text = c(leg_text, "DMMC offset")
        leg_symb = c(leg_symb, '+')
        leg_col = c(leg_col, col_DMMC)
    }
    if (!all(is.na(offset_kal))) {
        leg_text = c(leg_text, "Kalman filtered minima")
        leg_symb = c(leg_symb, '+')
        leg_col = c(leg_col, col_kal)
    }
    legend(x=Xrange[2]*0.8+Xrange[1]*0.2, y=Yrange[2], legend=leg_text, pch=leg_symb, col=leg_col)
    
    dev.off()
    
    return(0)
}
