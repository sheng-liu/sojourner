## filterTrack-methods
##
##
###############################################################################
##' @name filterTrack
##' @aliases filterTrack trimTrack trackLength
##' @title filterTrack
##' @rdname filterTrack-methods
##' @docType methods
##'
##' @description methods for filter and trim tracks based on track length.

##' @usage
##' filterTrack(trackll,filter=c(min=7,max=Inf))
##' trimTrack(trackll,trimmer=c(min=1,max=32))
##' trackLength(trackll)

##' @param trackll a list of track lists.
##' @param filter range of possible track lengths to keep
##' @param trimmer range of track lengths allowed in output, otherwise trimmed.


##' @return
##' \itemize{
##' \item{trackll} filtered or trimmed tracks.
##' \item{len} list of track lengths.
##' }

##' @details filterTrack() is used to filter out tracks that has length within 
##' a specified range (default 7~Inf). On the other hand, despite the lengths 
##' of tracks, trimTrack() is used to trim /cutoff all tracks to a specified 
##' range (default 1~32).
##'
##'
##' @examples
##' folder=system.file("extdata","SWR1",package="sojourner")
##' trackll=createTrackll(folder=folder, input=3)
##'
##' trackll.filter=filterTrack(trackll,filter=c(7,Inf))
##' trackll.trim=trimTrack(trackll,trimmer=c(1,20))
##'
##' # see the min and max length of the trackll
##' # trackLength() is a helper function output track length of trackll
##' lapply(trackLength(trackll),min)
##' lapply(trackLength(trackll.filter),min)
##'
##' lapply(trackLength(trackll),max)
##' lapply(trackLength(trackll.trim),max)
##'



##-----------------------------------------------------------------------------
## filterTrack

## a function to filter trackll based on specified fitler value 
## (filterTrack on track length), default 6 frames/steps to Inf

##' @export filterTrack
filterTrack=function(trackll,filter=c(min=7,max=Inf)){

    #filter=match.arg(filter)

    # reinforce name
    names(filter)=c("min","max")

    cat("applying filter, min",filter["min"],"  max",filter["max"],"\n")


    track.len=list()
    for (i in seq_along(trackll)){
        track.len[[i]]=sapply(trackll[[i]],function(track){dim(track)[1]})
        trackll[[i]]=trackll[[i]][ track.len[[i]]>=filter["min"] & 
                                       track.len[[i]]<filter["max"]]
    }

    return(trackll)
}


# no need for the focus swtich, as one can simply filter on a number that is 
# bigger than the dt he wanted to draw on

##-----------------------------------------------------------------------------
## trim long tracks into shorter ones

.trimTrack=function(track,min=1,max=32){
    t=track[min:max,]
    cc.t=complete.cases(t)
    t=t[cc.t,]
    return(t)
}

##' @export trimTrack
trimTrack=function(trackll,trimmer=c(min=1,max=32)){

    # reinforce name
    names(trimmer)=c("min","max")

    cat("applying trimmer, min",trimmer["min"],"  max",trimmer["max"],"\n")


    trackll.trim=list()
    length(trackll.trim)=length(trackll)
    names(trackll.trim)=names(trackll)
    for (i in seq_along(trackll)){
        trackll.trim[[i]]=lapply(trackll[[i]],.trimTrack,
                            min=trimmer["min"],max=trimmer["max"])
    }

    return(trackll.trim)

}

# max(sapply(x[[1]],dim))

##' @export trackLength
trackLength=function(trackll){
    len=list()

    length(len)=length(trackll)
    names(len)=names(trackll)

    for (i in seq_along(trackll)){
        len[[i]]=sapply(trackll[[i]],function(track){
            dim(track)[1]}
        )
    }
    return(len)
}

