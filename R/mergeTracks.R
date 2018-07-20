## mergeTracks-methods
##
##
###############################################################################
##' @name mergeTracks
##' @aliases mergeTracks
##' @title mergeTracks
##' @rdname mergeTracks-methods
##' @docType methods

##' @description merge track lists in a list into one

##' @usage
##' mergeTracks(folder, trackll)

##' @param folder Full path to the output files.
##' @param trackll A list of track lists.

##' @examples
##' #Basic masking with folder path with image masks
##' folder = "/FOLDERPATH"
##' trackll.merged <- mergeTracks(folder = folder, trackll = trackll)


##' @details
##' IMPORTANT: Once a trackll has been merged, it cannot be masked using maskTracks().
##' 
##' Merging creates a list list of data.frames. 
##' First level is the folder name, second level is a list of data.frames/tracks from all output files merged into one.
##' 
##' If not merged, track lists takes the name of individual files in the folder.
##' If merged, the single merged track list takes the folder name.

##' @export mergeTracks

##------------------------------------------------------------------------------
##

mergeTracks=function(folder, trackll){
    
    first.name = names(trackll)[[1]]
    pattern = substr(first.name, nchar(first.name)-3, nchar(first.name))
    
    track.holder=c()
    
    # getting a file list of Diatrack files in a directory
    file.list=list.files(path=folder,pattern=pattern,full.names=T)
    folder.name=basename(folder)
    
    # concatenate track list into one list of data.frames
    for (i in 1:length(file.list)){
        track.holder=c(track.holder,trackll[[i]])
    }
    
    # rename indexPerTrackll of index
    # extrac index
    Index=strsplit(names(track.holder),split="[.]")  # split="\\."
    
    # remove the last old indexPerTrackll
    Index=lapply(Index,function(x){
        x=x[1:(length(x)-1)]
        x=paste(x,collapse=".")})
    
    # add indexPerTrackll to track name
    indexPerTrackll=1:length(track.holder)
    names(track.holder)=mapply(paste,Index,
                               indexPerTrackll,sep=".")
    
    # make the result a list of list with length 1
    trackll=list()
    trackll[[1]]=track.holder
    names(trackll)[[1]]=folder.name
    
    cat(paste("\nMerging of folder", folder.name, "complete.\n", sep = " "))
    
    return(trackll)
}
