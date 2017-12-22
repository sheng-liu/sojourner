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
##' combineTrackll(total=2,merged=T)
##'
##' @param total The number of tracklls to combine.
##' @param merged An Logical indicate if the tracklls to combine are merged or not.
##' @return
##' \itemize{
##' \item{trackll:} combined trackll.
##' }
##' @details Combine multiple track listw (tracklls) from multiple folders into one trackll, i.e. combining track infomation
##'          from files in multiple folders (replicates) together as if they are in one folder. The tracklls
##'          either merged or un-merged.
##'
##'          Users will be prompted to input the name of each trackll to combine.
##'          If the combined tracklls are merged, users will also be prompted to input a combined attribute name for the combined trackll.
##'
##'
##' @examples
##'
##' # Generate trackll, and process,
##' # e.g. mask region of interest, merge tracks from multiple files.
##' folder1=system.file("extdata","SWR1",package="smt")
##' trackll1=createTrackll(interact=F,folder,input=1)
##' trackll1=maskTracks(folder,trackll)
##' trackll1=mergeTracks(folder,trackll)
##'
##' folder2=system.file("extdata","SWR1_2",package="smt")
##' trackll2=createTrackll(interact=F,folder,input=1)
##' trackll2=maskTracks(folder,trackll)
##' trackll2=mergeTracks(folder,trackll)
##'
##' # Combine the tracklls together, input trackll names when prompted,
##' trackll=combineTrackll(total=2,merged=T)

##' @export combineTrackll

#####################################################################################
#####################################################################################



combineTrackll<-function(total=2,merged=T){
  totalTracklls <- total
  totalTracklls <- as.numeric(totalTracklls)
  temp<-c()
  if(merged==T){
    for (i in 1:totalTracklls){
      trackll.n<-readline("Enter the name of the trackll needed to combine:   ")
      cat("Trackll",i,":  ",trackll.n,"\n")
      temp<-append(temp,get(paste(trackll.n))[[1]])
    }
    temp<-list(temp)
    names(temp)<-readline(cat("Set the attribute name of the combined trackll, \ne.g. H2Av_WT_combined trackll:   "))
  }
  else if(merged==F){
    for (i in 1:totalTracklls){
      trackll.n<-readline("Enter the name of the trackll needed to combine:   ")
      cat("Trackll",i,":  ",trackll.n,"\n")
      temp<-append(temp,get(paste(trackll.n)))
    }
  }
  else{
    cat("Wrong input. Start again.\n")
    stop()
  }
  return (temp)

}
