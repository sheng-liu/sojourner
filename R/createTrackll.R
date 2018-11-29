## createTrackll-methods
##
##
###############################################################################
##' @name createTrackll
##' @aliases createTrackll
##' @title createTrackll
##' @rdname createTrackll-methods
##' @docType methods
##'
##' @description take in Diatrack (.txt or .mat), ImageJ (.csv), SlimFast (.txt), or Utrack (.mat) input from a folder to output a list of track lists with the option to record frames and use multiple cores.

##' @usage 
##' createTrackll(folder, interact = FALSE, input = 1, ab.track = FALSE, 
##' cores = 1, frameRecord = TRUE)

##' @param interact Open interactive menu to choose the desired folder by selecting any file in it and select input type (script will process all files of that type in this folder).
##' @param folder Full path output file folder (if they are .txt, ensure that they are either all Diatrack or all SlimFast).
##' @param input Input file type (Diatrack .txt file = 1; Diatrack .mat session file = 2; ImageJ .csv file = 3; SlimFast .txt file = 4; Utrack .mat file = 5).
##' @param ab.track Use absolute coordinates for tracks.
##' @param cores Number of cores used for parallel computation. This can be the cores on a workstation, or on a cluster. Tip: each core will be assigned to read in a file when paralleled.
##' @param frameRecord Add a fourth column to the track list after the xyz-coordinates for the frame that coordinate point was found (especially helpful when linking frames). Highly recommended to leave on.
##' @return trackll
##' @details
##' 
##' (Note: When reading only Diatrack .mat sessipn files (input = 2), intensities will be saved after the frame column)
##' 
##' It is highly advised that the frame record option be left on to preserve the most information, especially when linking frames and when using Utrack.
##' If the frame record option is turned on for reading Diatrack .txt files (input = 1), take note that the frame record is artificially created as consecutive frames after the given start frame.
##' Otherwise, all other data types naturally record the frames of every coordinate point.
##'
##' The pre-censoring of single-frame tracks is dependent on the tracking software. For complete lossless track data, use Diatrack (.mat) session files.
##' If the initial creation of the trackll does not have a frame record, future exports and imports of the trackll will only preserve the start frames.
##'
##' If the cores are set to the maximum number of cores available on the system, the script may return a error after processing all the files. 
##' This error is due to the requirement of some systems to have one core open for system functions. 
##' This error will not affect the trackll output, but to avoid it, one can input one less than the maximum number of cores available.
##'
##' The naming scheme for the list of track list is as follows:
##' 
##' Track List: [full name of input file]
##' 
##' Track: [Last five characters of the file name].[Start frame].[Length].[Track].[Index in overall list (will differ from Track # when merging)]
##' 
##' (Note: The last five characters of the file name, excluding the extension, cannot contain ".")

##' @examples
##' # select track folder interactively and specify using 2 cores
##' #trackll <- createTrackll(interact = TRUE, cores = 2)
##' 
##' # Specify trackll folder path programably, specify the file format to be "2"
##' #(i.e. Diatrack txt files ending with .txt), cores to be 2.
##' folder=system.file("extdata","SWR1",package="sojourner")
##' trackll <- createTrackll(folder=folder, input = 1, cores = 2)

##' @export createTrackll

###############################################################################

### createTrackll ###

createTrackll=function(folder, interact = FALSE, input = 1, ab.track = FALSE, cores = 1, frameRecord = TRUE){
    
    #Interactive menu to select file in desired folder and input type
    if (interact == TRUE){
        cat("Choose one file in the folder for processing... \n")
        folder = dirname(file.choose());
        input = 0;
        if (input == 0){
            cat("Folder selection:", folder, "\n");
            cat("Enter input file type and press ENTER: \n")
            cat("1. Diatrack .txt file \n")
            cat("2. Diatrack .mat session file: \n")
            cat("3. ImageJ/MOSAIC or exported save .csv file\n")
            cat("4. SlimFast .txt file \n")
            cat("5. Utrack .mat file \n")
            input <- readline();
        }
    }

    #Error if no input
    if (input > 5 || input < 1){
        cat("Restart script with correct input.")
    }

    #Designate file types
    if (input == 1){
        return(readDiatrack(folder, ab.track = ab.track, cores = cores, frameRecord = frameRecord));
    } else if (input == 2){
        return(readDiaSessions(folder, ab.track = ab.track, cores = cores, frameRecord = frameRecord));
    } else if (input == 3){
        return(readParticleTracker(folder, ab.track = ab.track, cores = cores, frameRecord = frameRecord));
    } else if (input == 4){
        return(readSlimFast(folder, ab.track = ab.track, cores = cores, frameRecord = frameRecord));
    } else if (input == 5){
        return(readUtrack(folder, ab.track = ab.track, cores = cores, frameRecord = frameRecord));
    }
}
