## linkSkippedFrames-methods
##
##
###############################################################################
##' @name linkSkippedFrames
##' @aliases linkSkippedFrames
##' @title linkSkippedFrames
##' @rdname linkSkippedFrames-methods
##' @docType methods
##'
##' @description link trajectories skipped (or do not appear for) a number of frames

##' @usage 
##' linkSkippedFrames(trackll, tolerance, maxSkip, cores = 1)
##' 
##' @param trackll A list of track lists.
##' @param tolerance Distance tolerance level measured in pixels after the frame skip.
##' @param maxSkip Maximum number of frames a trajectory can skip.
##' @param cores Number of cores used for parallel computation. This can be the cores on a workstation, or on a cluster. Tip: each core will be assigned to read in a file when paralleled.
##' @return linked trackll
##' @details
##' Given user input for a tolerance level to limit how far the next point after the skip can deviate from the last point in pixel distance 
##' and a maximum number of frame skips possible, all trajectories falling within these parameters are automatically linked, renamed, and ordered accordingly. 
##' For a maxSkip example, if the maxSkip for a trajectory ending in frame 7 was 3, the next linked trajectory can start up to a maximum frame of 11.
##' 
##' Although not required, in order for the output to have a frame record column (recommended), the input must have one as well.
##' 
##' The naming scheme for each linked track is as follows:
##' 
##' [Last five characters of the file name].[Start frame #].[Length].[Track #].[# of links]
##' 
##' Track List: [full name of input file]
##' 
##' Track: [Last five characters of the file name].[Start frame].[Length].[Track].[# of links].[Index in overall list (will differ from Track # when merging)]
##' 
##' (Note: The last five characters of the file name, excluding the extension, cannot contain ".")

##' @examples
##' folder=system.file("extdata","SWR1",package="sojourner")
##' trackll=createTrackll(folder=folder, input=3)
##' #Basic function call of linkSkippedFrames
##' trackll.linked <- linkSkippedFrames(trackll, tolerance = 5, maxSkip = 10)

##' @export linkSkippedFrames

###############################################################################

#### .linkSkippedFrames ####

.linkSkippedFrames = function(track.list, tolerance, maxSkip){
    
    #Get abbreviated file name
    file.name = getTrackFileName(track.list);
    
    #Confirmation text of function call
    cat("\nLinking", file.name, "...\n");
    
    #Instantiate empty linked track list
    track.list.linked = list();
    
    #Instantiate frame/length/linknum list for track naming
    frame.list = list();
    length.list = list();
    linknum.list = list();
    
    #Instantiate starting frame of the last trajectory, max frame after skips, and index
    lastTrajFrame = getStartFrame(track.list, length(track.list));
    maxFrame = maxSkip + 1;
    trajectoryIndex = 1;
    
    #Extract file sub-name from given track list
    file.subname = substr(names(track.list[1]), 1, gregexpr(pattern = '\\.', names(track.list[1]))[[1]][1] - 1);
    
    #Loop through each new linked trajectory
    repeat{
        
        #Extract first trajectory in track list into data frame temp for linking trajectories
        temp <- track.list[[1]];
        
        #Record the start amd last frame of the starting trajectory from track name
        startFrame = getStartFrame(track.list, 1);
        lastFrame = getStartFrame(track.list, 1) + nrow(temp) - 1;
        
        #Delete the extracted first trajectory from track list and shift list to replace
        #The first trajectory in the track list is the next trajectory after the data frame temp
        track.list[[1]] <- NULL
        
        #Record the last X and Y coordinate values of data frame temp
        lastX = temp[[1]][[nrow(temp)]];
        lastY = temp[[2]][[nrow(temp)]];
        
        #Set link counter to 0;
        linkCounter = 0;
        
        #Instantiate index i and loop through trajectories in track list to look for frame skips 
        i = 1;
        repeat{
            
            #If the starting frame of the trajectory is beyond the maximum skips possible or the last frame, break
            nextFrame =  getStartFrame(track.list, i);
            if (i > length(track.list) || nextFrame > lastFrame + maxFrame){
                break;
            }
            
            #Record the first X and Y coordinate values of trajectory
            nextX = track.list[[i]][[1]][[1]];
            nextY = track.list[[i]][[2]][[1]];
            
            #Check if the distance difference between the first coordinate of trajectory and the last coordinate of temp are within the set tolerance
            if (sqrt((lastX-nextX)^2 + (lastY-nextY)^2) <= tolerance){
                
                #Update lastFrame
                lastFrame = getStartFrame(track.list, i) + nrow(track.list[[i]]) - 1;
                
                #Append trajectory to temp
                temp <- rbind(temp, track.list[[i]]);
                
                #Delete trajectory from track list and shift
                track.list[[i]] <- NULL
                
                #Recalculate last X and Y coordinate values of appended data frame temp
                lastX = temp[[1]][[nrow(temp)]];
                lastY = temp[[2]][[nrow(temp)]];
                
                #Increment link counter
                linkCounter = linkCounter + 1;
                
            } else {
                
                #Increment index to next trajectory if nothing is found
                i = i + 1;
            }
        }
        
        #Add start frame to frame list
        frame.list[[length(frame.list) + 1]] <- startFrame;
        
        #Add track length to length list
        length.list[[length(length.list) + 1]] <- nrow(temp);
        
        linknum.list[[length(linknum.list) + 1]] <- linkCounter;
        
        #Append data frame temp of the fully linked trajectory into track.list.linked and increment index
        track.list.linked[[trajectoryIndex]] <- temp;
        trajectoryIndex = trajectoryIndex + 1;
        
        #Break if the track list of decreasing length is empty 
        if (length(track.list) == 0){
            break;
        }
    }
    
    #Name track list:
    #[File name from input].[Start frame #].[Length].[Track #].[# of links]
    names(track.list.linked) = paste(file.subname, frame.list, length.list, c(1:length(track.list.linked)), linknum.list, sep=".");
    
    #Return linked track list and confirmation text
    cat("\n", Reduce("+", linknum.list), "links found in", file.name, "\n");
    
    return (track.list.linked);
}

#### linkSkippedFrames ####

linkSkippedFrames = function(trackll, tolerance, maxSkip, cores = 1){
    
    # detect number of cores
    max.cores=parallel::detectCores(logical=TRUE)
    
    if (cores == 1){
        link.trackll = lapply(trackll,function(x){
            .linkSkippedFrames(track.list = x, tolerance = tolerance, maxSkip = maxSkip)
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
                                varlist=c(".linkSkippedFrames", "tolerance", "maxSkip"),
                                envir=environment())
        
        link.trackll = parallel::parLapply(cl,trackll,function(x){
            .linkSkippedFrames(track.list = x, tolerance = tolerance, maxSkip = maxSkip)
        })
        
        # stop cluster
        cat("\nStopping clusters...\n")
        parallel::stopCluster(cl)
    }
    cat("\nProcess complete.\n")
    return (link.trackll);
    
}

#### Extra and useless functions ####

findMaxLinks = function(track.list, maxTolerance = 10, maxMaxSkip = 500){
    minLength = length(track.list)
    minTolerance = maxTolerance
    minSkip = maxMaxSkip
    cat("Searching...");
    for (i in 1:maxTolerance){
        for (j in 1:maxMaxSkip){
            length = length(linkSkippedFrames(track.list, tolerance = i, maxSkip = j));
            if (length < minLength){
                minLength = length
                minTolerance = i
                minSkip = j
            }
            
        }
    }
    cat("For the maximum number of links, with minimum tolerance and then minimum skips:")
    cat("\nTolerance =",  minTolerance)
    cat("\nSkips =",  minSkip)
    cat("\nLength =",  minLength)
}


