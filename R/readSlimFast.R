## readSlimFast-methods
##
##
###############################################################################
##' @name readSlimFast
##' @aliases readSlimFast
##' @title readSlimFast
##' @rdname readSlimFast-methods
##' @docType methods
##'
##' @description take in a SlimFast .txt session file as input, along with several other user-configurable parameters and output options, to return a track list of all the trajectories

##' @usage 
##' readSlimFast(folder, ab.track = FALSE, cores = 1, frameRecord = TRUE)
##' .readSlimFast(file, interact = FALSE,  ab.track = FALSE, frameRecord = FALSE)

##' @param folder Full path to SlimFast .txt files output folder.
##' @param ab.track Use absolute coordinates for tracks.
##' @param cores Number of cores used for parallel computation. This can be the cores on a workstation, or on a cluster. Tip: each core will be assigned to read in a file when paralleled.
##' @param frameRecord Add a fourth column to the track list after the xyz-coordinates for the frame that coordinate point was found (especially helpful when linking frames).
##' @param file Full path to track file.
##' @param interact Open menu to interactively choose file.
##' @return trackll
##' @details
##' The naming scheme for each track is as follows:
##' 
##' [Last five characters of the file name].[Start frame #].[Length].[Track #]
##' 
##' (Note: The last five characters of the file name, excluding the extension, cannot contain ".")

##' @examples
##' # Basic function call of .readSlimFast
##' # trackll <- readSlimFast(folder = folder, cores = 1)
##'
##' # Basic function call of .readSlimFast
##' # trackl <- .readSlimFast(interact = TRUE)
##' 

##' @export .readSlimFast
##' @export readSlimFast
###############################################################################

## TO DO: More efficient way of reading by using dplyr split

#### .readSlimFast ####

.readSlimFast = function(file, interact = FALSE,  ab.track = FALSE, frameRecord = FALSE){
    
    #Interactively open window
    if (interact == TRUE) {
        file=file.choose();
    }
    
    #Collect file name information
    file.name = basename(file);
    file.subname = substr(file.name, start=nchar(file.name)-8, stop=nchar(file.name)-4);
    
    #Display starter text
    cat("\nReading SlimFast file: ",file.name,"...\n");
    
    #Read first four columns of input SLIMFAST file
    data <- as.data.frame(subset(read.table(file,skip=1), select=c(1:4)));
    
    #Name columns and add z column of 1s in the appropriate location
    colnames(data) <- c("x","y","Frame", "track");
    data <- cbind(data, data.frame("z" = 1));
    data <- data[,c("x","y","z","Frame", "track")]
    
    #Instantiate track, start frame, and length lists
    track.list = list();
    frame.list = list();
    length.list = list();
    
    #Instantiate counter indexing variable and extract last trajectory
    counter = 1
    last.trajectory = data[nrow(data), ][[5]];
    
    #Loop through every given trajectory number
    for (i in 1:last.trajectory){
        
        #Create track data frame
        track <- data.frame();
        
        #Loop through every line of the input file with a globally updating counter index
        repeat{
            
            #If trajectory number is equal to the track number, add line to track and update counter
            if (i == as.integer(data[counter, ][[5]]) && counter <= nrow(data)){
                track <- rbind(track, data[counter, ]);
                counter = counter + 1;
            } else { #Break out of loop if they are unequal and the track ends
                break;
            }
        }
            
        #Remove track number column from track
        track <- track[-c(5)];
            
        #Option to add/remove frame record
        if (!frameRecord){
            track <- track[-c(4)];
        }
        
        #Calcualte absolute track coordinates if desired
        if (ab.track){
            track <- abTrack(track);
        }
        
        #Rename row names of track to appropriate index values
        print(track)
        print(nrow(track))
        rownames(track) <- sapply(1:nrow(track),toString);
        
        #Add start frame of track to frame list
        frame.list[[length(frame.list) + 1]] <- track[[4]][[1]];
        
        #Add track length to length list
        length.list[[length(length.list) + 1]] <- nrow(track);
        
        #Append temporary track data frame into track list
        track.list[[i]] <- track;
    }
    
    #Name track list:
    #[Last five characters of the file name without extension (cannot contain ".")].[Start frame #].[Length].[Track #]
    names(track.list) = paste(file.subname, frame.list, length.list, c(1:length(track.list)), sep=".");
    
    cat("\n", file.subname, "read and processed.\n")
    #Return track list
    return (track.list);
}

#### readDiaSessions ####

readSlimFast = function(folder, ab.track = FALSE, cores = 1, frameRecord = TRUE){
    
    trackll = list()
    track.holder = c()
    
    # getting a file list of Diatrack files in a directory
    file.list = list.files(path = folder, pattern = ".txt", full.names = TRUE)
    file.name = list.files(path = folder, pattern = ".txt", full.names = FALSE)
    folder.name=basename(folder)
    
    
    # read in tracks
    # list of list of data.frames,
    # first level list of file names and
    # second level list of data.frames
    
    max.cores = parallel::detectCores(logical=TRUE)
    
    if (cores == 1){
        
        for (i in 1:length(file.list)){
            
            track.list = .readSlimFast(file = file.list[i], ab.track = ab.track, frameRecord = frameRecord)
            
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
        parallel::clusterExport(cl,varlist=c(".readSlimFast","ab.track", "frameRecord"),envir=environment())
        
        # trackll=parallel::parLapply(cl,file.list,function(fname){
        trackll=parallel::parLapply(cl,file.list,function(fname){
            track=.readSlimFast(file=fname,ab.track=ab.track, frameRecord = frameRecord)
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
