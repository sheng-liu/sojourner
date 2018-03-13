#### combineTrackll.R
#### Wu Lab, Johns Hopkins University
#### Author: Xiaona Tang
#### Date: Sep 20, 2017

## combineTrackll-methods
##
###############################################################################
##' @name combineTrackll
##' @aliases combineTrackll
##' @title combineTrackll
##' @rdname combineTrackll-methods
##' @docType methods
##' @description Combine multiple tracklls into one trackll.
##'
##' @usage
##'
##' combineTrackll(trackll=c(trackll1,trackll2),name="combined trackll",merged=T)
##'
##' @param trackll The tracklls to be combined together.
##' @param name a character string given to set the "names" attribute for the combined trackll
##' @param merged An Logical indicate if the tracklls to combine are merged or not.
##' @return
##' \itemize{
##' \item{trackll:} combined trackll.
##' }
##' @details Combine multiple track lists (tracklls) from multiple folders into one trackll, i.e. combining track infomation
##'          from files in multiple folders (replicates) together as if they are in one folder. The tracklls can be
##'          either merged or un-merged.
##'          
##'          The name argument sets the "names" attribute for the combined trackll, which will be used in the same way
##'          as the folder names for the original tracklls, e.g., displayed as legend when plotting Dcoef or MSD for 
##'          the combined trackll. 
##'
##'
##' @examples
##'
##' # Generate trackll, and process, 
##' # e.g. mask region of interest, merge tracks from multiple files.
##' folder1=system.file("extdata","HSF",package="sojourner")
##' trackll1=createTrackll(interact=F,folder1,input=2, cores = 2)
##' trackll1=maskTracks(folder1,trackll1)
##' trackll1=mergeTracks(folder1,trackll1)
##' 
##' folder2=system.file("extdata","HSF_2",package="sojourner")
##' trackll2=createTrackll(interact=F,folder2,input=2, cores = 2)
##' trackll2=maskTracks(folder2,trackll2)
##' trackll2=mergeTracks(folder2,trackll2)
##'
##' # Combine the tracklls together, input trackll names when prompted,
##' trackll=combineTrackll(trackll=c(trackll1,trackll2),merged=T)

##' @export combineTrackll

#####################################################################################
#####################################################################################



combineTrackll<-function(trackll=c(trackll1,trackll2),name="combined trackll",merged=T){
    totalTracklls <- as.numeric(length(trackll))
    temp<-c()
    if(merged==T){
        for (i in 1:totalTracklls){
            temp<-append(temp,trackll[i][[1]])
        }
        temp<-list(temp)
        names(temp)<-name
    }
    else if(merged==F){
        for (i in 1:totalTracklls){
            temp<-append(temp,trackll[i])
        }
    }
    else{
        cat("Wrong input. Start again.\n")
        stop()
    }
    return (temp)
    
}
