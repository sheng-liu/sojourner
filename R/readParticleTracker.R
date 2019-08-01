
## readParticleTracker-methods
##' @name readParticleTracker
##' @aliases readParticleTracker
##' @title readParticleTracker
##' @rdname readParticleTracker-methods
##' @docType methods
##'
##' @description Read output file (tracks/trajectories in csv format) from 
##' ParticleTracker (a program of ImageJ plugin MosaicSuit).

##' @usage
##' readParticleTracker(folder, ab.track=FALSE, cores=1, frameRecord=TRUE)
##'
##'
## @method # this roxygen directive does not working
##' @param folder Full path to ImageJ .csv files output folder.
##' @param ab.track Use absolute coordinates for tracks.
##' @param cores Number of cores used for parallel computation. This can be 
##' the cores on a workstation, or on a cluster. Tip: each core will be 
##' assigned to read in a file when paralleled.
##' @param frameRecord Add a fourth column to the track list after the 
##' xyz-coordinates for the frame that coordinate point was found (especially 
##' helpful when linking frames).
##' @return trackll

##' @examples
##' # reading in tracks
##' #folder=system.file('extdata','ImageJ',package='sojourner')
##' #trackll=readParticleTracker(folder)
##' #str(trackll,max.level=2)

##' @details
##' The usage of readParticleTracker() is equivalent to readDiatrack().
##'
##' Note: the folder name should not contain '.', as it is a key character for 
##' subsequent indexing of file names.
##'
##' trackID=fileID.frameID.duration.indexPerFile.indexPerTrackll
##'
##' This 'indexPerFile' is the index within a track file.
##'
##' This 'indexPerTrackll' is the index within a trackll, which is unique.
##'
##' The macro used for generating the csv file is also included in ImageJ 
##' folder of the package: 
##' folder=system.file('extdata','ImageJ',package='sojourner')
##'

##' 
# @export readParticleTracker

############################################################################### 

##-----------------------------------------------------------------------------
## .readParticleTracker a function to read ParticleTracker (a program of
## ImageJ plugin MosaicSuite) output .csv file and returns a list of
## tracks


.readParticleTracker = function(file, interact = FALSE, ab.track = FALSE, 
    frameRecord = FALSE) {
    
    
    # interactively open window
    if (interact == TRUE) {
        file = file.choose()
    }
    
    file.name = basename(file)
    cat("\nReading ParticleTracker file: ", file.name, "...\n")
    
    # reading in data
    data = read.csv(file = file, header = TRUE)
    
    vars = c("Trajectory", "Frame", "x", "y", "z")
    track.data = dplyr::select(data, one_of(vars))
    
    track.list = split(track.data, f = track.data$Trajectory)
    # trackl.list=noquote(track.list) reset row.names of data.frame
    track.list = lapply(track.list, function(df) {
        row.names(df) = NULL
        return(df)
    })
    
    ## name the track.list the track list already has name, which is named
    ## in the particle tracker, to be compatible with downstream operation,
    ## need to first remove the name, and rename it.  names(track.list)=NULL
    
    # somehow flow thorugh the same code the readParticleTracker output
    # trackl has quotation the readDiatrack output trackl doesn't have
    # quotation the quotation can't be removed, and all works
    
    # init.frame.id
    init.fm = function(df) {
        df$Frame[1]
    }
    duration.fm = function(df) {
        length(df$Frame)
    }
    
    frame.id = sapply(track.list, init.fm)
    duration = sapply(track.list, duration.fm)
    
    
    ## TODO, modify the macro to remove the last affix .tif in xxx.tif.csv
    ## for now
    
    # file.subname=substr(file.name, start=nchar(file.name)-12,
    # stop=nchar(file.name)-8) # file.subname=substr(file.name, #
    # start=nchar(file.name)-8, # stop=nchar(file.name)-4)
    
    
    if (substr(file.name, start = nchar(file.name) - 7, 
        stop = nchar(file.name) - 4) == ".tif") {
        file.subname = substr(file.name, start = nchar(file.name) - 12, 
            stop = nchar(file.name) - 8)
    } else {
        file.subname = substr(file.name, start = nchar(file.name) - 8, 
            stop = nchar(file.name) - 4)
    }
    
    # file.id
    file.id = rep(file.subname, length(duration))
    
    # indexPerFile
    indexPerFile = seq(from = 1, to = length(duration))
    
    ## trackID=fileID.frameID.duration.indexPerFile
    track.name = paste(file.id, frame.id, duration, indexPerFile, sep = ".")
    # names(track.list)=noquote(track.name)
    names(track.list) = track.name
    
    # remove columns
    
    if (frameRecord) {
        subset.df = function(df) {
            df[, c("x", "y", "z", "Frame")]
        }
        track.list = lapply(track.list, subset.df)
    } else {
        subset.df = function(df) {
            df[, c("x", "y", "z")]
        }
        track.list = lapply(track.list, subset.df)
    }
    
    
    # convert normal trackll to ab.trackll for plotting can be paralleled
    # this has to be a seperate function here, as reading in it is not
    # going through loops.
    
    if (ab.track == TRUE) {
        
        cat("\nConverting to ab.trackl for plotting\n")
        abTrack = function(track) {
            data.frame(x = track$x - min(track$x), y = track$y - min(track$y))
        }
        ab.track.list = lapply(track.list, abTrack)
        
    }
    
    cat("\n", file.subname, "read and processed.\n")
    
    if (ab.track == TRUE) 
        return(ab.track.list) else return(track.list)
    
}


readParticleTracker = function(folder, ab.track = FALSE, cores = 1, 
    frameRecord = TRUE) {
    
    trackll = list()
    track.holder = c()
    
    # getting a file list of ParticleTracker files in a directory
    file.list = list.files(path = folder, pattern = ".csv", full.names = TRUE)
    file.name = list.files(path = folder, pattern = ".csv", full.names = FALSE)
    folder.name = basename(folder)
    
    # read in tracks list of list of data.frames, first level list of file
    # names and second level list of data.frames
    
    max.cores = parallel::detectCores(logical = TRUE)
    
    if (cores == 1) {
        
        # TODO: if cores are not assigned and detected cores more than one
        # automatic assign 1/4 of max.cores if (cores == 1 & max.cores>1)
        # switch cores=c(1,n,auto)
        
        for (i in seq_along(file.list)) {
            
            
            track = .readParticleTracker(file = file.list[i], 
                ab.track = ab.track, frameRecord = frameRecord)
            
            # add indexPerTrackll to track name
            indexPerTrackll = seq_along(track)
            names(track) = mapply(paste, names(track), indexPerTrackll, 
                sep = ".")
            
            trackll[[i]] = track
            names(trackll)[i] = file.name[i]
        }
        
    } else {
        
        # parallel this block of code assign each file to a CPU for reading in
        # using .readParticleTracker
        
        if (cores > max.cores) 
            stop(paste("Number of cores specified is", 
                "greater than recomended maxium: "), max.cores)
        
        cat("Initiated parallel execution on", cores, "cores\n")
        # use outfile=' to display result on screen
        cl <- parallel::makeCluster(spec = cores, type = "PSOCK", outfile = "")
        # register cluster
        parallel::setDefaultCluster(cl)
        
        # pass environment variables to workers
        
        # parallel::clusterExport(cl,varlist=c('.readParticleTracker',
        # 'ab.track'),envir=environment())
        # trackll=parallel::parLapply(cl,file.list,function(fname){
        # track=.readParticleTracker(file=fname,ab.track=ab.track)
        
        parallel::clusterExport(cl, varlist = c(".readParticleTracker", 
            "ab.track", "frameRecord"), envir = environment())
        
        trackll = parallel::parLapply(cl, file.list, function(fname) {
            track = .readParticleTracker(file = fname, ab.track = ab.track, 
                frameRecord = frameRecord)
            
            # add indexPerTrackll to track name
            indexPerTrackll = seq_along(track)
            names(track) = mapply(paste, names(track), indexPerTrackll, 
                sep = ".")
            return(track)
        })
        
        # stop cluster
        
        cat("\nStopping clusters...\n")
        
        parallel::stopCluster(cl)
        
        names(trackll) = file.name
        
    }
    
    cat("\nProcess complete.\n")
    
    return(trackll)
}




