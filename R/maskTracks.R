## maskTracks-methods
##
##
###############################################################################
##' @name maskTracks
##' @aliases maskTracks
##' @title maskTracks
##' @rdname maskTracks-methods
##' @docType methods

##' @description apply binary image masks to lists of track lists
#
##' @usage
##' maskTracks(folder, trackll)

##' @param folder Full path to the output files.
##' @param trackll A list of track lists.

##' @examples
##' #Basic masking with folder path with image masks
##' folder = "/FOLDERPATH"
##' trackll.masked <- maskTracks(folder = folder, trackll = trackll)
##'
##' #Compare the masking effect
##' plotTrackOverlay(trackll)
##' plotTrackOverlay(trackll.masked)
##'
##' #Plot mask
##' mask.list=list.files(path=folder,pattern="_MASK.tif",full.names=T)
##' plotMask(mask.file=mask.list[[1]])
##' 
##' #If Nuclear image is available
##' plotNucTrackOverlay(folder=folder,trackll=trackll)
##' plotNucTrackOverlay(folder=folder,trackll=trackll.masked)
##'
##' #Plot mask
##' plotMask(folder=folder)

##' @details
##' IMPORTANT: It will take an extremely long time to mask unfiltered data. Filter first using filterTrack(trackll,filter=c(min=7,max=Inf)), then mask using maskTracks(folder, trackll)!

##' Note the mask file should have the same name as the output files with a "_MASK.tif" ending. 
##' Users can use plotMask() and plotTrackOverlay() to see the mask and its effect on screening tracks.

##' @export maskTracks

##------------------------------------------------------------------------------
##
# read in mask and derive positive pix form a mask
maskPoint=function(mask.file,plot=F){
    
    mask.file.name=basename(mask.file)
    # read in tiff mask
    # library(rtiff)
    cat("Reading mask file    ",mask.file.name,"\n")
    mask=rtiff::readTiff(fn=mask.file)
    # plot(mask)
    
    pospt=which(mask@red!=0,arr.ind=T)
    pos.point=with(data.frame(pospt),data.frame(x=col,y=row))
    
    # horizontal is the same vertical is fliped as the pixels is counted from
    # upper left in the image, but is counted from lower left in the plot.
    if (plot) {
        plot(mask)
        plot(x=pos.point$x,y=pos.point$y)
        #ggplot(pos.point,aes(x=x,y=y))+geom_point()
    }
    
    
    return(pos.point)
}

##------------------------------------------------------------------------------
##
# Use each trajectory's geometric center as unit for clusterization.
# Each data point is an averaged trajectory.

trackCenter=function(trackll){
    
    # arithmetic mean for geometric center (centroid)
    track.center=list()
    length(track.center)=length(trackll)
    names(track.center)=names(trackll)
    
    for (i in 1:length(trackll)){
        track.center[[i]]=lapply(trackll[[i]],function(x){
            # round coords
            apply(x,2,function(coord){round(mean(coord))})})
    }
    
    return(track.center)
}


##------------------------------------------------------------------------------
##
# track.center and pos.point should be one to one cooresponding
posTracks=function(track.center,pos.point){
    
    # convert list to data.frame
    track.center.df=do.call(rbind.data.frame,track.center)
    names(track.center.df)=c("x","y","z", "Frame")
    pos.tracks=plyr::match_df(track.center.df,pos.point,on=c("x","y"))
    return(pos.tracks)
    
}

##------------------------------------------------------------------------------
##

maskTracks=function(folder, trackll){
    
    # read in mask
    maskl=list.files(path=folder,pattern="_MASK.tif",full.names=T)
    
    if (length(maskl)==0){
        cat("No image mask file ending '_MASK.tif' found.\n")
        
    }
    
    mask.track.index=list()
    length(mask.track.index)=length(trackll)
    names(mask.track.index)=names(trackll)
    
    masked.tracks=list()
    length(masked.tracks)=length(trackll)
    names(masked.tracks)=names(trackll)
    
    
    for (i in 1:length(trackll)){
        track.center=trackCenter(trackll)[[i]]
        pos.point=maskPoint(maskl[[i]],plot=F)
        mask.track.index[[i]]=posTracks(track.center,pos.point)
        
        index=rownames(mask.track.index[[i]])
        
        masked.tracks[[i]]=lapply(trackll[i],function(x){x[as.numeric(index)]})[[1]]
        
    }
    cat("\nAll files masked.\n")
    return(masked.tracks)
}