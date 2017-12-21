#### plotLocalizations.Dense.R
#### Wu Lab, Johns Hopkins University
#### Author: Xiaona Tang
#### Date: Dec 20, 2017

## plotLocalizations.Dense-methods
##
###############################################################################
##' @name plotLocalizations.Dense
##' @aliases plotLocalizations.Dense
##' @title plotLocalizations.Dense
##' @rdname plotLocalizations.Dense-methods
##' @docType methods
##' @description Plot localization map of molecules from a list of files in a folder with color coded 
##'              by local density of each molecule. 
##'              
##' @usage
##'
##' plotLocalizations.Dense(scale=256,r=125,file.No=0,point.scale=0.15)
##'
##' @param scale The pixel scale of image data.
##' @param r Radius within each molecule to calculate density.
##' @param file.No Select file(s) in the folder to plot. Default 0 for plotting all files in the folder.
##' @param point.scale Size of the dots representing the molecules.
##' @return
##' \itemize{
##' \item{PDF:} One PDF file with one plot on each page.
##' }
##' @details Plot localization map of molecules from a list of files in a folder with color coded 
##'              by local density of each molecule. The localization of molecule is considered as the first position of its track. 
##'              
##'          Upon running of the function, users will be prompted to input the name of the track list (trackll).
##'          Input un-merged trackll and the plotting will start.
##'          The local density of each molecule is calculated by counting the number of molecules within a given radius
##'          around the position of the molecule. The higher the number, the higher the local density.

##'
##' @examples
##'
##' # Generate trackll, and process, 
##' # e.g. mask region of interest, tracks from multiple files should not be merged.
##' folder=system.file("extdata","SWR1",package="smt")
##' trackll=createTrackll(interact=F,folder,input=0)
##' trackll=maskTracks(folder,trackll)
##' 
##' # Plot localization map,
##' plotLocalizations.Dense(scale=256,r=125,file.No=0,point.scale=0.15)


#####################################################################################
#####################################################################################


## Function for plotting molecule localizations with color coded by local molecule density.

.plotLocalizations.Dense<-function(scale=256,r=125,file.No=0,point.scale=0.15){
  
  library(sp)
  library(sampSurf)

  ## Import trackll (un-merged) information  
  if (length(file.No)>1&min(file.No)==0){
    stop("Wrong file.No input. Ceased.",call. = FALSE, domain = NULL)
  }
  
  trackll.label<-readline(cat("Enter the un-merged trackll you want to plot:    "))
  trackll<-get(paste(trackll.label))
  
  if (length(file.No)>=1&file.No[1]>0){
    trackll<-trackll[file.No]
  }
  
  ## Set plot area as black background and white frontground.
  oldpar <- par
  par(mar=c(3, 4, 4, 3),xpd=FALSE)
  
  par(mfrow=c(1,1),bg="black",fg="white")
  cat("Plotting...A PDF file will be output in the working directory.\n")
  
  ## Get trackl info and plot localization map for each file in the trackll.
  for (i in c(1:length(trackll))){
    plot.new()
    plot.window(xlim=c(0,scale*0.107),ylim=c(0,scale*0.107),xaxs = "i", yaxs = "i")
    
    axis(1,cex.axis=1,col.axis="white")
    axis(2,cex.axis=1,col.axis="white")
    mtext(gsub(".mat","",names(trackll[i])),side=3,line=0.5,cex=2,col.main="white")             
    mtext(expression(paste("X (",mu,"m)")),side=1,line=2,cex.lab=1,col="white")
    mtext(expression(paste("Y (",mu,"m)")),side=2,line=2,cex.lab=1,col="white")
    box()
    
    ## Get the localization of each track (molecule) as the first position of the track.
    localizations<-data.frame(matrix(ncol = 3, nrow = 0))
    for (j in c(1:length(trackll[[i]]))){
      localizations<-rbind(localizations,setNames(data.frame(paste0(names(trackll[[i]][j])),trackll[[i]][[j]]$x[[1]]*0.107,trackll[[i]][[j]]$y[[1]]*0.107), c("Trajectory","x", "y")))
    }
    ## Calculate local molecule density by counting the molecule number within a given radius.
    coordinates(localizations) <- c("x", "y")
    for (k in c(1:length(localizations$x))){
      sp.n = spCircle(r/1000, centerPoint=c(x=localizations[k,]$x,y=localizations[k,]$y), spID='tree.1') 
      count.n <- over(localizations, sp.n$spCircle)
      localizations$density[k]=sum(count.n,na.rm = T)
    }
    
    ## Generate color gradient/ramp based on the highest molecule density in each file.
    #cl=grDevices::heat.colors(max(localizations$density))
    #cl=colorspace::diverge_hsv(max(localizations$density))
    #cl=colorspace::diverge_hcl(max(localizations$density))
    cl <- colorRampPalette(c("blue4", "white", "red"))(n = max(localizations$density))
    
    ## plot the molecules as dots color coded by it's local density.
    for(i in c(1:length(localizations$x))){
      points(localizations[i,]$x,localizations[i,]$y,pch=16,col=cl[localizations[i,]$density],cex=point.scale)
    }
    ## Add color gradient legend to the right edge of each plot.
    .legend.col(col = cl, lev = localizations$density)
    ## Add radius and molecule number (n) as text legend to the topright corner of each plot.
    legend("topright",paste(rep(c("r = ","n = ")), rep(c(r,nrow(localizations))),rep(c(" nm",""))),col="white",bty = "n")
    
  }
  
  ## Reset plotting area parameters.
  par(oldpar)
  par(mfrow=c(1,1),bg="white",fg="black")
}




## Function for outputting the plots into one multipage PDF file in the working directory.

plotLocalizations.Dense<-function(scale=256, r=125, file.No=0, point.scale=0.15){
  ## Output the plots into one PDF file in the working directory.
  pdf(paste("plotLocalization.Dense.Heatmap--",format(Sys.time(),"%Y%m%d.%H%M%S"),".pdf",sep=""),width=11.7,height=11.7)

  .plotLocalizations.Dense(scale=scale, r=r, file.No=file.No, point.scale=point.scale)

  dev.off()
  
}

