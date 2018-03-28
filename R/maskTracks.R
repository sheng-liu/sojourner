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
##' 
##' indexCell(folder, trackll, areaFilter = c(0, Inf), intensityFilter = c(0, Inf), export = F, max.pixel = 128)
##' 
##' filterOnCell(trackll, numTracks = 0)
##' 
##' sampleTracks(trackll, num = 0)


##' @param folder Full path to the output files.
##' @param trackll A list of track lists.
##' @param max.pixel Pixel dimension of image
##' @param areaFilter Range of cell areas (pixel sq) to keep in filtering
##' @param intensityFilter Range of avg cell intensities (grayscale) to keep in filtering
##' @param numTracks Minimum number of required tracks in the trackll
##' @param num Number of tracks to randomly sample per trackl in trackll

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
##' If there are more mask files than trackll, masking will fail. If there are less mask files, trackls without masks will be deleted.
##' Users can use plotMask() and plotTrackOverlay() to see the mask and its effect on screening tracks.
##' 
##' indexCell() will mask a trackll, separate each cell into a trackl, display all cell areas and mean intensities, and then apply any area and intensity filters.
##' There is also the capability to export the final areas/intensities as "indexCell.csv" to the home directory and the pixel dimensions can be changed.
##' 
##' filterOnCell() eliminates all trackl in trackll that has less than numTracks tracks.
##' 
##' sampleTracks() randomly samples num number of tracks for each trackl in trackll.

##' @export maskTracks
##' @export indexCell
##' @export filterOnCell
##' @export sampleTracks

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
    
    # change rownames to number to avoid PC mac difference in rowname when do.call rbind
    rownames(track.center.df)=seq(1,dim(track.center.df)[1])
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
    
    if (length(maskl) > length(trackll)){
        stop("More masks than trackl.\n")
    }
    # make mask list and trackll one-to-one, delete extra trackl
    maskl.check = list()
    maskl.names = gsub("_MASK.tif","",basename(maskl))
    trackll.names = gsub("[.].*","",names(trackll))
    for (i in 1:length(trackll.names)) {
        found = FALSE
        for (j in 1:length(maskl.names)) {
            if (trackll.names[[i]] == maskl.names[[j]]) {
                maskl.check[[length(maskl.check)+1]] = maskl[j]
                found = TRUE
                break
            }
        }
        if (!found) {
            cat(paste(names(trackll)[[i]], "mask not found. Trackl deleted from masked output.\n", sep = " "))
            trackll[[i]] <- NULL
        }
    }
    maskl = maskl.check
    
    mask.track.index=list()
    length(mask.track.index)=length(trackll)
    names(mask.track.index)=names(trackll)
    
    masked.tracks=list()
    length(masked.tracks)=length(trackll)
    names(masked.tracks)=names(trackll)
    
    for (i in 1:length(trackll)){
        
        ## Finds all track centers
        track.center=trackCenter(trackll)[[i]]
        
        ## Returns all positive mask pixel locations
        pos.point=maskPoint(maskl[[i]],plot=F)
        
        ## Filters all positive track centers
        mask.track.index[[i]]=posTracks(track.center,pos.point)
        
        ## Collects track indexes to keep 
        index=rownames(mask.track.index[[i]])
        
        ## Filters only such indexes in the raw trackll[i]
        masked.tracks[[i]]=lapply(trackll[i],function(x){x[as.numeric(index)]})[[1]]
        
    }
    cat("\nAll files masked.\n")
    return(masked.tracks)
}

# Function masks, separates by cell, displays cell areas and mean intensities, and filters.
indexCell=function(folder, trackll, areaFilter = c(0, Inf), intensityFilter = c(0, Inf), export = F, max.pixel = 128){
    
    # Read in mask
    maskl=list.files(path=folder,pattern="_MASK.tif",full.names=T)

    # Read in nuclear image
    glowl=list.files(path=folder,pattern="_Nuclei.tif",full.names=T)
    
    if (length(maskl)==0){
        cat("No image mask file ending '_MASK.tif' found.\n")
        
    }
    
    if (length(maskl) > length(trackll)){
        stop("More masks than trackl.\n")
    }

    # Make mask list and trackll one-to-one and onto, delete extra trackl
    maskl.check = list()
    maskl.names = gsub("_MASK.tif","",basename(maskl))
    trackll.names = gsub("[.].*","",names(trackll))
    for (i in 1:length(trackll.names)) {
        found = FALSE
        for (j in 1:length(maskl.names)) {
            if (trackll.names[[i]] == maskl.names[[j]]) {
                maskl.check[[length(maskl.check)+1]] = maskl[j]
                found = TRUE
                break
            }
        }
        if (!found) {
            cat(paste(names(trackll)[[i]], "mask not found. Trackl deleted from masked output.\n", sep = " "))
            trackll[[i]] <- NULL
        }
    }
    maskl = maskl.check
    
    # Instantiate new variables
    masked.trackll=list()
    raw.names=names(trackll)
    cell.count = list()
    length(cell.count)=length(trackll)
    
    areas=list()
    intensities=list()
    
    ## Find all track centers
    track.center=trackCenter(trackll)
    
    # Loop through each trackl
    for (i in 1:length(trackll)){
        
        ## Returns all positive mask pixel locations
        invisible(capture.output(pos.point <- maskPoint(maskl[[i]],plot=F)))
        
        #Instantiate empty
        binary.mat = matrix( rep( 0, len=max.pixel*max.pixel), nrow = max.pixel)
        
        #Fill with binary pospoints
        for (m in 1:nrow(pos.point)){
            binary.mat[pos.point[[1]][[m]], pos.point[[2]][[m]]] = 1
        }
        
        # Calculate connected components and label accordingly
        labeled.mat <- SDMTools::ConnCompLabel(binary.mat)
        
        # FOR VISUALIZING THE LABELED MASK
        # image(t(labeled.mat[128:1,]),col=c('grey',rainbow(length(unique(labeled.mat))-1)))
        
        pos.points.indexed = list()
        length(pos.points.indexed) = max(labeled.mat)

        # Convert labeled.mat to lists of pos.points per cell
        for (j in 1:nrow(pos.point)) {
            pos.points.indexed[[labeled.mat[pos.point[j,]$x, pos.point[j,]$y]]] <- rbind( pos.points.indexed[[labeled.mat[pos.point[j,]$x, pos.point[j,]$y]]], pos.point[j,])
        }

        ## Display Area ##
        cat(maskl.names[[i]], "\n", "Cell Areas (pixels squared): ", "\n", sep = "")
        areal <- sapply(pos.points.indexed, nrow)
        cat(cat(areal, sep = "\n"))
        areas <- append(areas, areal) 

        
        ## Display Intensity ##
        glow = t(matrix(pixmap::getChannels(rtiff::readTiff(glowl[[i]])), nrow = max.pixel, ncol = max.pixel))
        raw.intensities = rep(list(list()), max(labeled.mat))
        for(row in 1:nrow(labeled.mat)) {
            for(col in 1:ncol(labeled.mat)) {
                if (labeled.mat[row, col] != 0){
                    raw.intensities[[labeled.mat[row, col]]][[length(raw.intensities[[labeled.mat[row, col]]])+1]] = glow[row, col]
                }
            }
        }
        cat("Cell Intensities (grayscale): ", "\n", sep = "")
        intensityl <- list()
        for (cell in 1:length(raw.intensities)){
            intensityl[[cell]] <- mean(unlist(raw.intensities[[cell]]))
            cat(paste(intensityl[[cell]], "\n"))
        }
        cat("\n")
        intensities <- append(intensities, intensityl) 


        # Update cell count per trackl
        cell.count[[i]] <- length(pos.points.indexed)
        
        # Loop through each cell
        for (k in 1:length(pos.points.indexed)){
            
            mask.track.index=list()
            length(mask.track.index)=length(trackll)
            names(mask.track.index)=names(trackll)
            
            ## Filters all positive track centers
            mask.track.index[[k]]=posTracks(track.center[[i]], pos.points.indexed[[k]])
            
            ## Collects track indexes to keep 
            index=rownames(mask.track.index[[k]])
            
            ## Filters only such indexes in the raw trackll[i]
            masked.trackll[[length(masked.trackll)+1]]=lapply(trackll[i],function(x){x[as.numeric(index)]})[[1]]
        }
        
    }
    
    #print(areas, row.names = FALSE)
    #print(intensities, row.names = FALSE)

    # Edit the names with concatenated cell number to trackl names
    masked.names = list()
    for (c in 1:length(cell.count)){
        for (v in 1:cell.count[[c]]){
            masked.names[[length(masked.names)+1]] <- paste(raw.names[c], toString(v), sep = "_")
        }
    }
    names(masked.trackll) <- masked.names
    
    if (export){
        df <- data.frame(matrix(ncol = 3, nrow = 0))
        colnames(df) <- c("Cell", "Area (pixel sq)", "Intensity (grayscale)")
    }
    
    # Filter
    j = 1
    for(i in 1:length(areas)){
        if (areas[[i]] > areaFilter[[2]] || areas[[i]] < areaFilter[[1]] || intensities[[i]] > intensityFilter[[2]] || intensities[[i]] < intensityFilter[[1]]){
            masked.trackll[[j]] <- NULL;
        } else {
            if (export){
                df[nrow(df) + 1,] = list(names(masked.trackll)[[j]], areas[[i]], intensities[[i]])
            }
            j = j + 1
        }
    }
    cat("\nOnly cell with areas between ", areaFilter[[1]], " and ", areaFilter[[2]], " pixels sq kept.", sep = "")
    cat("\nOnly cell with avg intensities between ", intensityFilter[[1]], " and ", intensityFilter[[2]], " grayscale kept.", sep = "")
    if (export){
        write.csv(df, file = "indexCell.csv")
        cat("\nindexCell.csv exported to working directory.")
    }
    
    cat("\n\nAll files masked, separated by cell indexes, and filters applied.\n")
    return(masked.trackll)
}

filterOnCell=function(trackll, numTracks = 0){
    if (numTracks == 0) {
        cat("\nEnter a lower limit for number of tracks.\n")
    } else {
        i = 1
        while(i != length(trackll)){
            if (length(trackll[[i]]) < numTracks){
                trackll[[i]] <- NULL
            } else {
                i = i + 1
            }
        }
    }
    return(trackll)
}

sampleTracks = function(trackll, num = 0){
    if (num == 0) {
        cat("\nEnter sample size.\n")
    } else {
        for (i in 1:length(trackll)){
            trackll[[i]] <- sample(trackll[[i]], num)
        }
        return(trackll)
    }
}