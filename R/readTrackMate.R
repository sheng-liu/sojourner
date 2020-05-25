
## readTrackMate-methods
##
##
###############################################################################
##' @name readTrackMate
##' @aliases readTrackMate
##' @title readTrackMate
##' @rdname readTrackMate-methods
##' @docType methods
##'
##' @description Read output file (tracks/trajectories in csv format) from TrackMate (a program of ImageJ plugin).

##' @usage
##' readTrackMate(folder, ab.track = F, cores = 1, frameRecord = T)
##'
##' .readTrackMate(file, interact = F, ab.track = F, frameRecord = F)
##'
## @method # this roxygen directive does not working
##' @param folder Full path to ImageJ .csv files output folder.
##' @param ab.track Use absolute coordinates for tracks.
##' @param cores Number of cores used for parallel computation. This can be the cores on a workstation, or on a cluster. Tip: each core will be assigned to read in a file when paralleled.
##' @param frameRecord Add a fourth column to the track list after the xyz-coordinates for the frame that coordinate point was found (especially helpful when linking frames).
##' @param file Full path to Diatrack .mat session file.
##' @param interact Open menu to interactively choose file.
##' @examples
##' ## Not run:
##' trackMate.file.url="https://drive.google.com/file/d/1abYagKpu1yea_7FM8tJ7nkQdOktFyS6j/view?usp=sharing"
##' # with the upgrade of cloud protocol, this command does not work anymore, 
##' # download and save file to a location on local computer 
##' # download.file(trackMate.file.url, "~/test/test.csv")
##' # dir.create("~/test/")
##' # folder=~/test/"
##' # trackll=readTrackMate(folder)
##' # str(trackll,max.level=2)
##' 
##' @details
##' Note: the folder name should not contain ".", as it is a key character for subsequent indexing of file names.
##'
##' trackID=fileID.frameID.duration.indexPerFile.indexPerTrackll
##'
##' This "indexPerFile" is the index within a diatrackFile.
##'
##' This "indexPerTrackll" is the index within a trackll, which is unique.
##' 
##' The usage of readTrackMate() is similar to ReadDiatrack(), please reference for more details.

##' @export readTrackMate
##' @export .readTrackMate

###############################################################################

##------------------------------------------------------------------------------
## .readTrackMate
## a function to read TrackMate (a program of ImageJ plugin) output .csv file and returns a list of tracks

.readTrackMate = function(file, interact = FALSE, ab.track = FALSE, 
                                frameRecord = FALSE) {
    
    
    # interactively open window
    if (interact == TRUE) {
        file = file.choose()
    }
    
    file.name = basename(file)
    cat("\nReading trajectory file: ", file.name, "...\n")
    
    # reading in data
    data = read.csv(file = file, header = TRUE)
    
    names(data)[names(data) == "TRACK_ID"] <- "Trajectory"
    names(data)[names(data) == "POSITION_X"] <- "x"
    names(data)[names(data) == "POSITION_Y"] <- "y"
    names(data)[names(data) == "POSITION_Z"] <- "z"
    names(data)[names(data) == "FRAME"] <- "Frame"
    
    # if (frameRecord) {
    #     names(data)[names(data) == "FRAME"] <- "Frame"}

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
    
    # TODO: old comment, still valid?
    # somehow flow thorugh the same code the read output
    # trackl has quotation the readDiatrack output trackl doesn't have
    # quotation the quotation can't be removed, and all works
    
    # init.frame.id
    init.fm = function(df) {
        df$Frame[1]
    }
    duration.fm = function(df) {
        length(df$Frame)
    }
    
    frame.id = vapply(track.list, init.fm, integer(1))
    duration = vapply(track.list, duration.fm, integer(1))
    
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
    
    # remove Frame columns when needed
    # dplyr::select already selected three or four columns
    # it is either no change, or remove Frame column
    # 
    # # TODO: remove above code trunk for other reading in functions
    # # replace it with adding new ones.
    # # in the description, this is not add new by calculating, it is select
    # existing frame to be included in calculaton. 'Add' infers for data that
    # doesnt have frame number, add new ones accordingly. this
    # manipulation is unecessary, can by defualt add frame, when there is no
    # frame in the data, add new ones.

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


readTrackMate=function(folder,ab.track=F,cores=1, frameRecord=T){
    
    trackll=list()
    track.holder=c()
    
    # getting a file list of Diatrack files in a directory
    file.list=list.files(path=folder,pattern=".csv",full.names=T)
    file.name=list.files(path=folder,pattern=".csv",full.names=F)
    folder.name=basename(folder)
    
    # read in tracks
    # list of list of data.frames,
    # first level list of file names and
    # second level list of data.frames
    
    max.cores=parallel::detectCores(logical=T)
    
    if (cores == 1){
        
        # TODO: if cores are not assigned and detected cores more than one
        # automatic assign 1/4 of max.cores
        # if (cores == 1 & max.cores>1)
        # switch cores=c(1,n,auto)
        
        for (i in 1:length(file.list)){
            
            
            track=.readTrackMate(file=file.list[i],ab.track=ab.track, frameRecord = frameRecord)
            
            # add indexPerTrackll to track name
            indexPerTrackll=1:length(track)
            names(track)=mapply(paste,names(track),indexPerTrackll,sep=".")
            
            trackll[[i]]=track
            names(trackll)[i]=file.name[i]
        }
        
    }else{
        
        # parallel this block of code
        # assign each file to a CPU for reading in using .readTrackMate
        
        if (cores>max.cores)
            stop("Number of cores specified is greater than recomended maxium: ",max.cores)
        
        cat("Initiated parallel execution on", cores, "cores\n")
        # use outfile=" to display result on screen
        cl <- parallel::makeCluster(spec=cores,type="PSOCK",outfile="")
        # register cluster
        parallel::setDefaultCluster(cl)
        
        # pass environment variables to workers
        
        parallel::clusterExport(cl,varlist=c(".readTrackMate","ab.track", "frameRecord"),envir=environment())
        
        trackll=parallel::parLapply(cl,file.list,function(fname){
            track=.readTrackMate(file=fname,ab.track=ab.track, frameRecord = frameRecord)
            
            # add indexPerTrackll to track name
            indexPerTrackll=1:length(track)
            names(track)=mapply(paste,names(track),indexPerTrackll,sep=".")
            return(track)
        })
        
        # stop cluster
        
        cat("\nStopping clusters...\n")
        
        parallel::stopCluster(cl)
        
        names(trackll)=file.name
        
    }
    
    cat("\nProcess complete.\n")
    
    return(trackll)
}

