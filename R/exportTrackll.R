## exportTrackll-methods
##
##
###############################################################################
##' @name exportTrackll
##' @aliases exportTrackll
##' @title exportTrackll
##' @rdname exportTrackll-methods
##' @docType methods
##'

##' @description take in a list of track lists (trackll) and export it into row-wise (ImageJ/MOSAIC) .csv files in the working directory

##' @usage 
##' exportTrackll(trackll, cores = 1)
##' 
##' .exportRowWise(track.list)

##' @param trackll A list of track lists.
##' @param cores Number of cores used for parallel computation. This can be the cores on a workstation, or on a cluster. Tip: each core will be assigned to read in a file when paralleled.
##' @param track.list A track list (a list of trajectory data frames).

##' @details
##' The reason why ImageJ/MOSAIC style .csv export was chosen is because it fully preserves track frame data, while maintaining short computation time and easy readability in Excel/etc.
##' 
##' In order to import this .csv export back into a trackll at any point (while preserving all information), select input = 3 in createTrackll.
##' 
##' If the track list does not have a fourth frame record column (not recommended), it will just output the start frame of each track instead and will take noticeably longer.
##' 
##' It is not recommended that exportTrackll be run on merged list of track lists (trackll).
##' Also, ensure that the input trackll is a list of track lists and not just a track list.
##' 
##' The naming scheme for each export is as follows:
##' 
##' [Last five characters of the file name]_[yy-MM-dd]_[HH-mm-ss].csv

##' @examples
##' #Basic function call to exportTrackll with 2 cores into current directory
##' exportTrackll(trackll, cores = 2)
##' 
##' #Export one track list
##' .exportRowWise(trackl)
##' 
##' #Get current working directory
##' getwd()
##' 
##' #Import export save back into a trackll
##' trackll.2 <- createTrackll(folder = getwd(), input = 3, cores = 2)
##' @importFrom rowr cbind.fill
##' @export .exportRowWise
##' @export exportTrackll

###############################################################################

#### .exportRowWise ####

.exportRowWise = function(track.list){
    
    #Confirmation text of function call
    cat("\nWriting .csv row-wise output in current directory for", getTrackFileName(track.list), "...\n");
    
    #Collect track file name
    track.file.name <- getTrackFileName(track.list);
    
    #Check for frame record column
    if (ncol(track.list[[1]]) != 3){
        
        #Rename track list as trajectory numbers
        names(track.list) <- c(1:length(track.list));
        
        #Combine track data frames by trajectory and reorder column
        df <- bind_rows(track.list, .id = "Trajectory")[, c("Trajectory", "Frame", "x", "y", "z")]
    
    } else{
        #Empty data frame df to be written into the .csv
        df <- NULL;
        
        #Loop through every trajectory in input track.list
        for (i in 1:length(track.list)){
            
            #Create a data frame temp with trajectory, frame, and track coordinate data 
            temp <- data.frame("Trajectory" = i, "Frame" = getStartFrame(track.list, i), track.list[[i]][1:3]);

            #Append data frame df with data frame temp
            df <- rbind(df, temp);
        }
    }
    #Write the data frame df into the .csv and display confirmation text
    file.name = paste(track.file.name, format(Sys.time(), format = "_%y-%m-%d_%H-%M-%S"), ".csv", sep = "")
    write.csv(df, file=file.name);
    cat(paste("\n", file.name, " placed in current directory.\n", sep =""))
}

#### .exportColWise ####

#Function unused for now- may be helpful when lossy Diatrack .txt export is neccessary

#Install packages and dependencies
#library(rowr)

.exportColWise = function(track.list){
    
    #Confirmation text of function call
    cat("\nWriting .csv column-wise output in current directory for", getTrackFileName(track.list), "...\n");

    frame.list <- list()
    
    #Loop through every trajectory in input track.list
    for (i in 1:length(track.list)){
        
        start.frame = getStartFrame(track.list[i])
        
        frame.list <- c(frame.list, start.frame, 0, 0)
        
        temp <- track.list[[i]][1:3]
        
        if (i != 1){
            df <- cbind.fill(df, temp, fill = 0) 
        } else {
            df <- temp
        }
    }
    
    colnames(df) <- frame.list
    
    #header = "format (columnwise): Frame1 row n+1: (y(tn) x(tn) z(tn)), row n+1: (y(t(n+1)) x(t(n+1)) z(t(n+1))), row n+2: (y(t(n+2)) x(t(n+2) z(t(n+2)) y(t(n+3)).... where Frame1 is the frame number where the target is seen for the first time, and the columns define trajectories. Beware! the number of tracks is limited by the width of the widest text file on your machine. Rowwise export preferred"

    #Write the data frame df into the .csv and display confirmation text
    file.name = paste("COL", getTrackFileName(track.list), ".csv", sep = "")
    #write(header, file = file.name, append = T)
    write.table(df, file = file.name, row.names = FALSE, sep = ","); #,append = T
    
    cat(paste("\n", file.name, "placed in current directory.\n\n", sep =""))
}

#### exportTrackll ####

exportTrackll = function(trackll, cores = 1){
    
    # detect number of cores
    max.cores=parallel::detectCores(logical=T)
    
    if (cores == 1){
        export = lapply(trackll,function(x){
            .exportRowWise(track.list = x)
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
                                varlist=c(".exportRowWise"),
                                envir=environment())
        
        export = parallel::parLapply(cl,trackll,function(x){
            .exportRowWise(track.list = x)
        })
        
        # stop cluster
        cat("\nStopping clusters...\n")
        parallel::stopCluster(cl)
    }
    
}
