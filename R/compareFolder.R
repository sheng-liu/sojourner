
## compareFolder-methods
##
##
###############################################################################
##' @name compareFolder
##' @aliases compareFolder
##' @title compareFolder
##' @rdname compareFolder-methods
##' @docType methods
##' 
##' @description compare folders with Diatrack output files. merge track files
##' in each folder into one item of a track list. This list can then be fed
##' into other functions for comparison. It keeps folder information as names
##' of the resulting list.

##' @usage compareFolder(folders,input=1,ab.track=FALSE,cores=1)
##' @param folders a vector storing paths to the folders location."..."
##'   indicates multiple (unlimited) folders can be added into the function.
##' @param input Input file type (Diatrack .txt file = 1; Diatrack .mat session
##'   file = 2; ImageJ .csv file = 3; SlimFast .txt file = 4).
##' @param ab.track a Logical indicating if absolute coordinates should be 
##' used.
##' @param cores Number of cores used for parallel computation. This can be the
##'   cores on a workstation, or on a cluster.

##' @return
##' \itemize{
##' \item{trackll} A list of tracks, each item of a track list correspond to a
##' folder. This list can then be fed into other functions for comparison. }

##' @examples
##' folder1=system.file("extdata","SWR1",package="sojourner")
##' folder2=system.file("extdata","HTZ1",package="sojourner")
##' trackll=compareFolder(folders=c(folder1,folder2), input=3)
##' str(trackll,max.level=1)

##' @import reshape2
##' @export compareFolder
##'
## FUTURE: maybe plot on dt
###############################################################################


compareFolder=function(folders, input=1, ab.track=FALSE,cores=1){

    # the number of folder to compare can be extended using ... statement
    # folder.list=list(folder1,folder2,folder3,folder4,folder5)
    # remove null folders by subsetting un-null folders
    #     null.folder=sapply(folder.list,is.null)
    #     folder.list=folder.list[!null.folder]

    folder.list=sapply(folders,list,simplify=TRUE)

    names(folder.list)=sapply(folder.list,basename)

    sample.list=list()


    if (ab.track == TRUE){

        for (i in seq_along(folder.list)) {
            sample.list[i] = mergeTracks(
                folder=folder.list[[i]], 
                createTrackll(folder=folder.list[[i]],input = input, 
                            ab.track=TRUE,cores=cores))
            cat("\n...\n") # seperator makes ouput clearer
            names(sample.list)[i]=names(folder.list)[i]
        }

    }else{

        for (i in seq_along(folder.list)) {
            # i=1
            sample.list[i] = mergeTracks(
                folder=folder.list[[i]], 
                createTrackll(folder=folder.list[[i]],
                            input = input, ab.track=FALSE,cores=cores))
            cat("\n...\n") # seperator makes ouput clearer
            names(sample.list)[i]=names(folder.list)[i]
        }
    }

    #names(sample.list)=names(folder.list)
    return(sample.list)

}











