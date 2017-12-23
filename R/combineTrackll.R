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
##' combineTrackll(trackll=c(trackll1,trackll2),merged=T)
##'
##' @param total The number of tracklls to combine.
##' @param merged An Logical indicate if the tracklls to combine are merged or not.
##' @return
##' \itemize{
##' \item{trackll:} combined trackll.
##' }
##' @details Combine multiple track listw (tracklls) from multiple folders into one trackll, i.e. combining track infomation
##'          from files in multiple folders (replicates) together as if they are in one folder. The tracklls can be
##'          either merged or un-merged.
##'          
##'          If the combined tracklls are merged, users will be prompted to input a combined attribute name for the combined trackll. 
##'
##'
##' @examples
##'
##' # Generate trackll, and process, 
##' # e.g. mask region of interest, merge tracks from multiple files.
##' folder1=system.file("extdata","HSF",package="smt")
##' trackll1=createTrackll(interact=F,folder1,input=2)
##' trackll1=maskTracks(folder1,trackll1)
##' trackll1=mergeTracks(folder1,trackll1)
##' 
##' folder2=system.file("extdata","HSF_2",package="smt")
##' trackll2=createTrackll(interact=F,folder2,input=2)
##' trackll2=maskTracks(folder2,trackll2)
##' trackll2=mergeTracks(folder2,trackll2)
##'
##' # Combine the tracklls together, input trackll names when prompted,
##' trackll=combineTrackll(trackll=c(trackll1,trackll2),merged=T)

##' @export combineTrackll

#####################################################################################
#####################################################################################



combineTrackll<-function(trackll=c(trackll1,trackll2),merged=T){
  totalTracklls <- as.numeric(length(trackll))
  temp<-c()
  if(merged==T){
    for (i in 1:totalTracklls){
      temp<-append(temp,trackll[i][[1]])
    }
    temp<-list(temp)
    names(temp)<-readline(cat("Set the attribute name of the combined trackll, \ne.g. HSF_WT_combined trackll:   "))
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
