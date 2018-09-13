## readUtrack-methods
##
##
###############################################################################
##' @name readUtrack
##' @aliases readUtrack
##' @title readUtrack
##' @rdname readUtrack-methods
##' @docType methods

##' @description take in a single channel Utrack file as input, along with several other user-configurable parameters and output options, to return a track list of all the trajectories found in the file

##' @usage 
##' readUtrack(folder, ab.track = F, cores = 1, frameRecord = T)
##' 
##' .readUtrack(file, interact = F, ab.track = F, frameRecord = F)

##' @param folder Full path to Utrack files output folder.
##' @param ab.track Use absolute coordinates for tracks.
##' @param cores Number of cores used for parallel computation. This can be the cores on a workstation, or on a cluster. Tip: each core will be assigned to read in a file when paralleled.
##' @param frameRecord Add a fourth column to the track list after the xyz-coordinates for the frame that coordinate point was found (almost mandatory for Utrack).
##' @param file Full path to Utrack file.
##' @param interact Open menu to interactively choose file.

##' @details
##' The naming scheme for each track is as follows:
##' 
##' [Last five characters of the file name].[Start frame #].[Length].[Track #]
##' 
##' (Note: The last five characters of the file name, excluding the extension, cannot contain ".")

##' @examples
##' #Basic function call of .readUtrack
##' #trackll <- readUtrack(folder = /FILEPATH/, cores = 2)
##' 
##' #Basic function call of .readUtrack
##' #trackl <- .readUtrack(interact = TRUE)

##' @export .readUtrack
##' @export readUtrack

##' @importFrom R.matlab readMat

###############################################################################

#### readUtrack ####

.readUtrack = function(file, interact = F, ab.track = F, frameRecord = F){
    
    #Interactively open window
    if (interact == TRUE) {
        file=file.choose();
    }
    
    #Collect file name information
    file.name = basename(file);
    file.subname = substr(file.name, start=nchar(file.name)-8, stop=nchar(file.name)-4);
    
    #Display starter text
    cat("\nReading Utrack file: ",file.name,"...\n");
    
    #Port .mat data
    data <- readMat(file)$tracksFinal;
    
    #Stop, since data has to be in sections of 3
    if (length(data) %% 3 != 0){
        stop("Error in data file.\n")
    }
    
    num.tracks = length(data)/3
    track.list = list();
    frame.list = list();
    length.list = list();
    
    #Loop through each track
    for (i in 1:num.tracks){
        #Track data exists every 3 lists, starting from list 2
        counter = i*3 - 1;
        
        #Track data exists every 3 lists, starting from list 3
        frame.data = data[[counter+1]]
        
        #Frame data has to be a 2x4 data frame (length = 8)
        #Second column has to be 1 (for start frame) and then 2 (for end frame)
        if (length(frame.data) != 8 || frame.data[[3]] != 1 || frame.data[[4]] != 2){
            stop("Only 1 channel readings accepted.\n")
        }
        start.frame = frame.data[[1]]
        end.frame = frame.data[[2]]
        
        track <- data.frame("x" = numeric(), "y" = numeric(), "z" = integer());
        
        #Add data to track
        frame = start.frame

        track.length = length(data[[counter]])/8
        for(j in 1:track.length){
            #xyz coordinates start every 8 frames
            point = 8*j-7
            x = data[[counter]][[point]]
            y = data[[counter]][[point+1]]
            z = data[[counter]][[point+2]]
            
            #Skip if NaN (denoting skipped frame track)
            if (!is.nan(x)){
                if (frameRecord){
                    track <- rbind(track, data.frame("x" = x, "y" = y, "z" = z, "Frame" = frame));
                } else {
                    track <- rbind(track, data.frame("x" = x, "y" = y, "z" = z));
                }
            }
            #Count frames
            frame = frame + 1;
        }
        track.list[[i]] <- track;
        frame.list[[length(frame.list) + 1]] <- start.frame;
        length.list[[length(length.list) + 1]] <- end.frame - start.frame + 1
        if (ab.track){
            track <- abTrack(track);
        }
    }
    #Name track list:
    #[Last five characters of the file name without extension (cannot contain ".")].[Start frame #].[Length].[Track #]
    names(track.list) = paste(file.subname, frame.list, length.list, c(1:length(track.list)), sep=".");
    
    #File read and processed confirmation text
    cat("\n", file.subname, "read and processed.\n")
    
    #Return track list
    return(track.list);
} 

#### readUtrack ####

readUtrack = function(folder, ab.track = F, cores = 1, frameRecord = T){
    
    trackll = list()
    track.holder = c()
    
    # getting a file list of Diatrack files in a directory
    file.list = list.files(path = folder, pattern = ".mat", full.names = T)
    file.name = list.files(path = folder, pattern = ".mat", full.names = F)
    folder.name=basename(folder)
    
    
    # read in tracks
    # list of list of data.frames,
    # first level list of file names and
    # second level list of data.frames
    
    max.cores = parallel::detectCores(logical=T)
    
    if (cores == 1){
        
        for (i in 1:length(file.list)){
            
            track.list = .readUtrack(file = file.list[i], ab.track = ab.track, frameRecord = frameRecord)
            
            # add indexPerTrackll to track name
            indexPerTrackll = 1:length(track.list)
            names(track.list) = mapply(paste, names(track.list), indexPerTrackll,sep = ".")
            
            trackll[[i]] = track.list
            names(trackll)[i] = file.name[i]
        }
        
    } else {
        
        # parallel this block of code
        # assign reading in using .readDiatrack to each CPUs
        
        # detect number of cores
        # FUTURE: if more than one, automatic using multicore
        
        if (cores>max.cores)
            stop("Number of cores specified is greater than recomended maximum: ", max.cores)
        
        cat("Initiated parallel execution on", cores, "cores\n")
        # use outfile="" to display result on screen
        cl <- parallel::makeCluster(spec = cores,type = "PSOCK", outfile = "")
        # register cluster
        parallel::setDefaultCluster(cl)
        
        # pass environment variables to workers
        parallel::clusterExport(cl,varlist=c(".readUtrack","ab.track", "frameRecord"),envir=environment())
        
        # trackll=parallel::parLapply(cl,file.list,function(fname){
        trackll=parallel::parLapply(cl,file.list,function(fname){
            track=.readUtrack(file=fname,ab.track=ab.track, frameRecord = frameRecord)
            # add indexPerTrackll to track name
            indexPerTrackll=1:length(track)
            names(track)=mapply(paste,names(track),indexPerTrackll,sep=".")
            return(track)
        })
        
        # stop cluster
        cat("\nStopping clusters...\n")
        parallel::stopCluster(cl)
        
        names(trackll)=file.name
        # names(track)=file.name
        
    }
    
    cat("\nProcess complete.\n")
    
    return(trackll)
}

