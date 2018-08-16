## segmentTrackll-methods
##
###############################################################################
## segmentTrackll
###############################################################################
##' @name segmentTrackll
##' @aliases segmentTrackll
##' @title segmentTrackll
##' @rdname segmentTrackll-methods
##' @docType methods
##' @description 
##' generate track segments of desired length with given trackll variable.
##' segment length of 6 or greater is needed. Tracks that are shorter than the
##' indicated seg value will be sorted out.

#helper function for single track segmentation. Used in segmentTrackll
segmentTrack=function(track,seg){
    segs=list()
    track.len = dim(track)[1]
    window=1:seg
    maxStart = track.len-seg
    for (i in 0:maxStart){
        segs[[i+1]] <- track[window+i,]
    }
    names(segs) = paste("",1:(maxStart+1),sep = "")
    return(segs)
}

#segments all tracks in trackll in to segments of indicated length
segmentTrackll=function(trackll,seg=7){
    min.len.trackll=filterTrack(trackll,filter=c(min=seg,max=Inf))
    t.length = length(trackll)
    segmentll=list()
    for (i in 1: t.length){
        tl.length=min.len.trackll[[i]]
        segmentl = lapply(min.len.trackll[[i]], segmentTrack, seg=seg)
        mapply(function(trackname,segl){names(segl) = paste(trackname,names(segl),seg="")}, names(min.len.trackll[[i]]), segmentl)
        segmentl = unlist(segmentl, recursive = FALSE)
        segmentll[[i]] = segmentl
    }
    names(segmentll)=names(trackll)
    return(segmentll)
}

#dcoef calculation for each segment in segmentll
find.dcoef=function(segmentll, dcoef){
    slopes = lapply(dcoef, function(x){x[,"slope"]})
    for (i in 1:2) {
        nn = names(segmentll[[i]])
        for(j in 1:length(nn)){
            print(nn[j])
            print(slopes[[i]][nn[j]])
        }
    }
}

#cutoff value calculation Currently not in use.
getCutoff=function(minVal, maxVal){
    #assuming no ouliers, we use max and min
    low = (minVal*minVal*maxVal) ^ (1/3)
    high = (minVal*maxVal*maxVal) ^ (1/3)
    return (c(low,high))
}

#calculates average displacement, displacement variance, and dcoef for all tracks in dcoef.
getSpecs=function(segmentll){
    segSize = dim(segmentll[[1]][[1]])[1]
    MSD = msd(segmentll, dt=segSize-1,filter = c(min=segSize,max=Inf))
    sll.dcoef = Dcoef(MSD=MSD, dt=segSize-1,filter = c(min=segSize-1,max=Inf), rsquare = 0)
    sll.dispvar = dispVariance(segmentll, min=segSize)
    sll.avgdisp = lapply(displacement.trackll(segmentll), function(x){reshape2::melt(x$SummarizedDisplacement)})
    specList = list()
    for (i in 1:length(sll.avgdisp)) {
        newframe = sll.avgdisp[[i]]
        dispvar.frame = reshape2::melt(sll.dispvar[[i]])
        rownames(dispvar.frame) = dispvar.frame$L1
        newframe = cbind(newframe, dispvar = dispvar.frame$value)
        rownames(newframe) = dispvar.frame$L1
        newframe = cbind(newframe, dcoef = NA) #now this is unnecessary since rsquare=0 and it wont generate NAs
        slots = match(rownames(sll.dcoef[[i]]), rownames(dispvar.frame))
        newframe[,"dcoef"][slots]=sll.dcoef[[i]][,"slope"]
        newframe$L1 = NULL
        newframe$L2 = NULL
        names(newframe) = c("avgdisp","dispvar","dcoef")
        specList[[i]] = newframe
    }
    return(specList)
}

#returns vector with one of "B","F","I" indicating the state of the track.
seg.distinguish=function(specl){
    min.avgdisp=min(specl$avgdisp)
    max.avgdisp=max(specl$avgdisp)
    min.dispvar=min(specl$dispvar)
    max.dispvar=max(specl$dispvar)
    cut.avgdisp=getCutoff(min.avgdisp, max.avgdisp)
    cut.dispvar=getCutoff(min.dispvar, max.dispvar)
    
    avgdisp.bound = specl$avgdisp < sqrt(min.avgdisp*max.avgdisp)
    dispvar.bound = specl$dispvar < sqrt(min.dispvar*max.dispvar)
    dcoef.bound = specl$dcoef <= 0
    boundvec = avgdisp.bound & dispvar.bound
    
    avgdisp.free = specl$avgdisp > sqrt(min.avgdisp*max.avgdisp)
    dispvar.free = specl$dispvar > sqrt(min.dispvar*max.dispvar)
    dcoef.free = specl$dcoef > 0
    freevec = avgdisp.free & dispvar.free & dcoef.free
    
    boundvec[boundvec] = "B"
    boundvec[freevec] = "F"
    boundvec[boundvec ==FALSE] = "I"
    
    return(boundvec)
}

#returns an extended form of segll with all the characteristics, trackname
#and a rough-guess classification.
seg.classify=function(segmentll){
    donelist = list()
    specs = getSpecs(segmentll)
    for (i in 1:length(segmentll)){
        statevec = seg.distinguish(specs[[i]])
        segl = segmentll[[i]]
        for (j in 1:length(segl)){
            segl[[j]] = cbind(segl[[j]],avgdisp=specs[[i]][,"avgdisp"][[j]])
            segl[[j]] = cbind(segl[[j]],dispvar=specs[[i]][,"dispvar"][[j]])
            segl[[j]] = cbind(segl[[j]],dcoef=specs[[i]][,"dcoef"][[j]])
            segl[[j]] = cbind(segl[[j]],state=statevec[j])
            segl[[j]] = cbind(segl[[j]],name=names(segmentll[[i]])[j])
            segl[[j]] = cbind(segl[[j]],tracknum=strsplit(names(segmentll[[i]])[j],".",fixed=T)[[1]][5])
        }
        donelist[[i]]=segl
    }
    return(donelist)
}

#export classified segements: modified version of exportTrackll
.exportRowWise2=function(track.list){
    
    #Confirmation text of function call
    cat("\nWriting .csv row-wise output in current directory for", getTrackFileName(track.list), "...\n");
    
    #Collect track file name
    track.file.name <- getTrackFileName(track.list);
    
    #Check for frame record column
    if (ncol(track.list[[1]]) == 11){
        
        #Rename track list as trajectory numbers
        names(track.list) <- c(1:length(track.list));
        
        #Combine track data frames by trajectory and reorder column
        df <- bind_rows(track.list, .id = "Trajectory")[, c("Trajectory", "Frame", "x", "y", "z", "avgdisp", "dispvar", "dcoef", "state", "name", "tracknum")]
        
    } else{ ###THIS PART WAS AVOIDED
        #Empty data frame df to be written into the .csv
        df <- NULL;
        
        #Loop through every trajectory in input track.list
        for (i in 1:length(track.list)){
            
            #Create a data frame temp with trajectory, frame, and track coordinate data 
            temp <- data.frame("Trajectory" = i, "Frame" = getStartFrame(track.list, i), track.list[[i]][1:3], track.list[[i]][,"state"]);
            
            #Append data frame df with data frame temp
            df <- rbind(df, temp);
        }
    }
    #Write the data frame df into the .csv and display confirmation text
    file.name = paste(track.file.name, format(Sys.time(), format = "_%y-%m-%d_%H-%M-%S"), ".csv", sep = "")
    write.csv(df, file=file.name);
    cat(paste("\n", file.name, " placed in current directory.\n", sep =""))
}

#export classified segements: modified version of exportTrackll
exportSegll = function(trackll, cores = 1){
    
    # detect number of cores
    max.cores=parallel::detectCores(logical=T)
    
    if (cores == 1){
        export = lapply(trackll,function(x){
            .exportRowWise2(track.list = x)
        })
    } else {
        # parallel excecute above block of code
        if (cores>max.cores)
            stop("Number of cores specified is greater than maxium: ",
                 max.cores)
        
        cat("Initiated parallel execution on", cores, "cores\n")
        
        # use outfile="" to display result on screen
        cl <- parallel::makeCluster(spec=cores,type="PSOCK",outfile="")
        # register cluster
        parallel::setDefaultCluster(cl)
        
        # pass environment variables to workers
        parallel::clusterExport(cl,
                                varlist=c(".exportRowWise2"),
                                envir=environment())
        
        export = parallel::parLapply(cl,trackll,function(x){
            .exportRowWise2(track.list = x)
        })
        
        # stop cluster
        cat("\nStopping clusters...\n")
        parallel::stopCluster(cl)
    }
}

#Used to generate pdf output file
plotSeg=function(segll){
    ab.segll=convert.abtrackll(segll)
    plotTrack(ab.segll, frame.min = 1)
}