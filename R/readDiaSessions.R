#### readDiaSessions.R
#### Wu Lab, Johns Hopkins University
#### Author: Sun Jay Yoo
#### Date: July 10, 2017

## readDiaSessions-methods
##
##
###############################################################################
##' @name readDiaSessions
##' @aliases readDiaSessions
##' @title readDiaSessions
##' @rdname readDiaSessions-methods
##' @docType methods

##' @description take in a Diatrack .mat session file as input, along with several other user-configurable parameters and output options, to return a track list of all the trajectories found in the session file
                                         
##' @usage 
##' readDiaSessions(folder, ab.track = F, cores = 1, frameRecord = T)
##' 
##' .readDiaSessions(file, interact = F, ab.track = F, frameRecord = F)

##' @param folder Full path to Diatrack .mat session files output folder.
##' @param ab.track Use absolute coordinates for tracks.
##' @param cores Number of cores used for parallel computation. This can be the cores on a workstation, or on a cluster. Tip: each core will be assigned to read in a file when paralleled.
##' @param frameRecord Add a fourth column to the track list after the xyz-coordinates for the frame that coordinate point was found (especially helpful when linking frames).
##' @param file Full path to Diatrack .mat session file.
##' @param interact Open menu to interactively choose file.

##' @details
##' The naming scheme for each track is as follows:
##' 
##' [Last five characters of the file name].[Start frame #].[Length].[Track #]
##' 
##' (Note: The last five characters of the file name, excluding the extension, cannot contain “.”)

##' @examples
##' #Basic function call of .readDiaSessions
##' trackll <- readDiaSessions(folder = /FILEPATH/, cores = 2)
##' 
##' #Basic function call of .readDiaSessions
##' trackl <- .readDiaSessions(interact = T)

##' @export .readDiaSessions
##' @export readDiaSessions

##' @importFrom R.matlab readMat

###############################################################################

#### Note ####

#This script takes Diatrack .mat files as input, and returns a list of data frames (a track list) of all the particle trajectories.
#The aim is to optimize and un-censor this process, instead of having to use MATLAB to extract a large .txt file which is then fed into R.

#Additional features:
#Adding frame records, removing frame records, outputing column-wise and row-wise to .csv files, linking skipped frames

#### Testing ####

#A .mat session file with 10117 frames was used to test both scripts.

#Using the MATLAB script, a 272.6MB .txt file was first created and was then fed into the readDiatrack() script to output track lists. 
#Automating this process using "matlabr" resulted in 4488 censored tracks (should be 4487 tracks since the script does not censor first frame) in 3:48 mins.

#Using readDiaSessions, the intermediate .txt file was no longer needed to be created and the session file directly results in track lists.
#This script resulted in 34689 uncensored tracks in 2:01 mins. 

#### readDiaSessions ####

#Install packages and dependencies
#install.packages("R.matlab")
#library(R.matlab)

.readDiaSessions = function(file, interact = F, ab.track = F, frameRecord = F){
    
    #Interactively open window
    if (interact == TRUE) {
        file=file.choose();
    }
    
    #Start timer
    #if (timer == TRUE) {
    #    start.time = Sys.time()
    #}
    
    #Collect file name information
    file.name = basename(file);
    file.subname = substr(file.name, start=nchar(file.name)-8, stop=nchar(file.name)-4);
    
    #Display starter text
    cat("\nReading Diatrack session file: ",file.name,"...\n");
    
    #Pre-process data (for both newer and older session file versions)
    #Successor and predecessor rows of first frame switched for consistency
    #(Unsure why Diatrack reverses the ordering of these two rows for the first frame)
    data <- readMat(file)$tracks;
    if (length(data[1][[1]][[1]]) == 7){
        temp <- data[1][[1]][[1]][[7]]; 
        data[1][[1]][[1]][[7]] <- data[1][[1]][[1]][[6]];
        succ = 7
    } else if (length(data[1][[1]][[1]]) == 8){
        temp <- data[1][[1]][[1]][[8]]; 
        data[1][[1]][[1]][[8]] <- data[1][[1]][[1]][[7]];
        succ = 8;
    } else {
        cat("ERROR: Use a different Diatrack version.")
    }
    data[1][[1]][[1]][[6]] <- temp;
        pred = 6;

    #Data structure of data for future reference:
    #data[FRAME][[1]][[1]][[ROW]][[COL]] 
    
    #Instantiate indexing variables, the track list, and the start frame list
    startIndex = 0;
    startFrame = 1;
    trajectoryIndex = 1;
    track.list = list();
    frame.list = list();
    length.list = list();
    
    #Loop for each trajectory track to be saved into track list
    repeat{
        
        #Loop until the next particle without a predecessor is found
        repeat{
            
            #Basic iteration through indexes and then through frames
            if (startIndex == 0 || startIndex < ncol(data[[startFrame]][[1]][[1]])){
                startIndex = startIndex + 1;
            } else {
                
                
                startFrame = startFrame + 1;
                startIndex = 1;
            }
            
            #Check at each iteration
            if (startFrame > length(data)){ #Break at end frame
                break;
            } else if (length(data[startFrame][[1]][[1]][[1]]) == 0){ #Iterate to next frame at empty frames
                next;
            } else if (data[startFrame][[1]][[1]][[pred]][[startIndex]] == 0) { #Break if particle is found
                break;
            }
            #Do nothing and iterate to next indexed particle if no particle is found
        }
        
        #Break track loop at end frame
        if (startFrame > length(data)){
            break;
        }
        
        #Instantiate initial frame and index coordinates into looping frame and index coordinates
        frame = startFrame;
        index = startIndex;
        
        #Create temporary track to insert into track list
        track <- data.frame("x" = numeric(), "y" = numeric(), "z" = integer());
        
        #Loop through every instance the particle exists and add its data to track
        #Break once it no longer has successors
        repeat{
            RefinedCooX = round(data[frame][[1]][[1]][[2]][[index]], 2);
            RefinedCooY = round(data[frame][[1]][[1]][[1]][[index]], 2);
            RefinedCooZ = round(data[frame][[1]][[1]][[3]][[index]], digits = 1);
            if (frameRecord){
                track <- rbind(track, data.frame("x" = RefinedCooX, "y" = RefinedCooY, "z" = RefinedCooZ, "Frame" = frame));
            } else {
                track <- rbind(track, data.frame("x" = RefinedCooX, "y" = RefinedCooY, "z" = RefinedCooZ));
            }
            if (data[frame][[1]][[1]][[succ]][[index]] != 0) {
                index = data[frame][[1]][[1]][[succ]][[index]];
                frame = frame + 1;
            } else
                break;
        }
        
        #Add start frame to frame list
        frame.list[[length(frame.list) + 1]] <- startFrame;
            
        #Add track length to length list
        length.list[[length(length.list) + 1]] <- nrow(track);
            
        #Calcualte absolute track coordinates if desired
        if (ab.track){
            track <- abTrack(track);
        }
        
        #Append temporary track for particle into track list and iterate to the next trajectory
        track.list[[trajectoryIndex]] <- track;
        trajectoryIndex = trajectoryIndex + 1;
    }
    #Name track list:
    #[Last five characters of the file name without extension (cannot contain ".")].[Start frame #].[Length].[Track #]
    names(track.list) = paste(file.subname, frame.list, length.list, c(1:length(track.list)), sep=".");
    
    #File read and processed confirmation text
    cat("\n", file.subname, "read and processed.\n")
    
    #Display stop timer
    #if (timer == TRUE) {
    #    end.time = Sys.time();
    #    time.taken = end.time - start.time;
    #   cat("Duration: ");
    #    cat(time.taken);
    #    cat(" mins\n");
    #}
    
    #Return track list
    return(track.list);
} 

#### readDiaSessions ####

readDiaSessions = function(folder, ab.track = F, cores = 1, frameRecord = T){
    
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
    
    max.cores = parallel::detectCores(logical = F)
    
    if (cores == 1){
        
        for (i in 1:length(file.list)){
            
            track.list = .readDiaSessions(file = file.list[i], ab.track = ab.track, frameRecord = frameRecord)
            
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
        parallel::clusterExport(cl,varlist=c(".readDiaSessions","ab.track", "frameRecord"),envir=environment())
        
        # trackll=parallel::parLapply(cl,file.list,function(fname){
        trackll=parallel::parLapply(cl,file.list,function(fname){
            track=.readDiaSessions(file=fname,ab.track=ab.track, frameRecord = frameRecord)
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
