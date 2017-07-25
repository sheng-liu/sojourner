
## readParticleTracker-methods
##
##
###############################################################################
##' @name readParticleTracker
##' @aliases readParticleTracker
##' @title readParticleTracker
##' @rdname readParticleTracker-methods
##' @docType methods
##'
##' @description Read output file (tracks/trajectories in csv format) from ParticleTracker (a program of ImageJ plugin MosaicSuit).

##' @usage
##' readParticleTracker(folder, merge = F, ab.track = F, mask = F, cores = 1, frameRecord = T)
##'
##' .readParticleTracker(file, interact = F, ab.track = F, frameRecord = F)
##'
## @method # this roxygen directive does not working
##' @param folder Full path to Diatrack output file.
##' @param merge An logical indicate if the output list should be merged into one. Default merge = FALSE, output list is divided by file names.
##' @param mask An logical indicate if image mask should be applied to screen tracks. Default False. Note the mask file should have the same name as the Diatrack output txt file with a "_MASK.tif" ending. Users can use plotMask() and plotTrackOverlay() to see the mask and its effect on screening tracks.
##' @param cores Number of cores used for parallel computation. This can be the cores on a workstation, or on a cluster. Tip: each core will be assigned to read in a file when paralelled.

##' @param frameRecord Add a fourth column to the track list after the xyz-coordinates for the frame that coordinate point was found (especially helpful when linking frames).

##' @return
##' \itemize{
##' \item{merge = F} Defult. A list of list of data.frames. First level is a list of file names in Diatrack output folder, second level is a list of data.frames from individual output file.

##' \item{merge = T} A list list of data.frames. First level is the folder name. second level is a list of data.frames from all Diatrack output files merged into one
##' }


## @section Usage : {
## readParticleTracker(folder,merge=F)
## }

##' @examples
##' # reading in tracks
##' folder=system.file("extdata","ImageJ",package="smt")
##' trackll=readParticleTracker(folder)
##' str(trackll,max.level=2)
##'
##' # masking with image mask
##' trackll.masked=readParticleTracker(folder=folder,merge=F,mask=T)
##' str(trackll.masked,1)
##'
##' # compare the masking effect
##' plotTrackOverlay(trackll)
##' plotTrackOverlay(trackll.masked)
##'
##' # if Nucclear image is available
##' plotNucTrackOverlay(folder=folder,trackll=trackll)
##' plotNucTrackOverlay(folder=folder,trackll=trackll.masked)
##'
##' # plot mask
##' plotMask(folder=folder)



##' @details
##' The usage of readParticleTracker() is equivalent to ReadDiatrack().
##' default merge = FALSE, so the researcher can assay variations between files. Keep both output as two level list is for simplicity of downstream analysis.
##'
##' Note: the folder name should not contain ".", as it is a key character for subsequent indexing of file names.
##'
##' trackID=fileID.frameID.duration.indexPerFile.indexPerTrackll
##'
##' This "indexPerFile" is the index within a diatrackFile.
##'
##' This "indexPerTrackll" is the index within a trackll, which is unique.
##'
##' The macro used for generating the csv file is also included in ImageJ folder of the package: folder=system.file("extdata","ImageJ",package="smt")
##'

##' @export readParticleTracker
##' @export .readParticleTracker

###############################################################################

##------------------------------------------------------------------------------
## .readParticleTracker
## a function to read ParticleTracker (a program of ImageJ plugin MosaicSuit) output .csv file and returns a list of tracks


.readParticleTracker=function(file,interact=F,ab.track=F, frameRecord=F){


    # interactively open window
    if (interact==T) {
        file=file.choose()
    }

    file.name=basename(file)
    cat("\nReading ParticleTracker file: ",file.name,"...\n")

    # reading in data
    data=read.csv(file=file,header=T)

    vars=c("Trajectory","Frame","x","y","z")
    track.data=dplyr::select(data,one_of(vars))

    track.list=split(track.data,f=track.data$Trajectory)
    # trackl.list=noquote(track.list)
    # reset row.names of data.frame
    track.list=lapply(track.list,function(df){row.names(df)=NULL;return(df)})

    ## name the track.list
    # the track list already has name, which is named in the particle tracker,
    # to be compatible with downstream operation, need to first remove the name,
    # and rename it.
    # names(track.list)=NULL

    # somehow flow thorugh the same code
    # the readParticleTracker output trackl has quotation
    # the readDiatrack output trackl doesn't have quotation
    # the quotation can't be removed, and all works

    # init.frame.id
    init.fm=function(df){df$Frame[1]}
    duration.fm=function(df){length(df$Frame)}

    frame.id=sapply(track.list,init.fm)
    duration=sapply(track.list,duration.fm)


    ##TODO, modify the macro to remove the last affix .tif in xxx.tif.csv
    # for now

    # file.subname=substr(file.name,
    #                     start=nchar(file.name)-12,
    #                     stop=nchar(file.name)-8)
    # #     file.subname=substr(file.name,
    # #                         start=nchar(file.name)-8,
    # #                         stop=nchar(file.name)-4)


    if (substr(file.name, start=nchar(file.name)-7,stop=nchar(file.name)-4) == ".tif"){
        file.subname=substr(file.name, start=nchar(file.name)-12,stop=nchar(file.name)-8)
    } else {
        file.subname=substr(file.name, start=nchar(file.name)-8, stop=nchar(file.name)-4)
    }

    # file.id
    file.id=rep(file.subname,length(duration))

    # indexPerFile
    indexPerFile=seq(from=1,to=length(duration))

    ## trackID=fileID.frameID.duration.indexPerFile
    track.name=paste(file.id,frame.id,duration,indexPerFile,sep=".")
    #names(track.list)=noquote(track.name)
    names(track.list)=track.name

    # remove columns

    if (frameRecord){
        subset.df=function(df){df[,c("x","y","z", "Frame")]}
        track.list=lapply(track.list,subset.df)
    } else {
        subset.df=function(df){df[,c("x","y","z")]}
        track.list=lapply(track.list,subset.df)
    }


    # convert normal trackll to ab.trackll for plotting
    # can be paralleled
    # this has to be a seperate function here, as reading in it is not going through loops.

    if (ab.track==T) {

        cat ("\nConverting to ab.trackl for plotting\n")
        abTrack=function(track){
            data.frame(x=track$x-min(track$x),
                       y=track$y-min(track$y))
        }
        ab.track.list=lapply(track.list,abTrack)

    }

    cat("\n", file.subname, "read and processed.\n")

    if (ab.track==T) return(ab.track.list) else return(track.list)

}


readParticleTracker=function(folder,merge= F,ab.track=F,mask=F,cores=1, frameRecord=T){

    trackll=list()
    track.holder=c()

    # getting a file list of Diatrack files in a directory
    file.list=list.files(path=folder,pattern=".csv",full.names=T)
    file.name=list.files(path=folder,pattern=".csv",full.names=F)
    folder.name=basename(folder)

    # read in mask
    mask.list=list.files(path=folder,pattern="_MASK.tif",full.names=T)

    if (mask==T & length(mask.list)==0){
        cat("No image mask file ending '_MASK.tif' found.\n")

    }

    # read in tracks
    # list of list of data.frames,
    # first level list of file names and
    # second level list of data.frames

    max.cores=parallel::detectCores(logical=F)

    if (cores==1){

        # TODO: if cores are not assigned and detected cores more than one
        # automatic assign 1/4 of max.cores
        # if (cores==1 & max.cores>1)
        # switch cores=c(1,n,auto)

        for (i in 1:length(file.list)){


            track=.readParticleTracker(file=file.list[i],ab.track=ab.track, frameRecord = frameRecord)

            # add indexPerTrackll to track name
            indexPerTrackll=1:length(track)
            names(track)=mapply(paste,names(track),indexPerTrackll,sep=".")

            trackll[[i]]=track
            names(trackll)[i]=file.name[i]
        }

    }else{

        # parallel this block of code
        # assign each file to a CPU for reading in using .readParticleTracker

        if (cores>max.cores)
            stop("Number of cores specified is greater than recomended maxium: ",max.cores)

        cat("Initiated parallel execution on", cores, "cores\n")
        # use outfile=" to display result on screen
        cl <- parallel::makeCluster(spec=cores,type="PSOCK",outfile="")
        # register cluster
        parallel::setDefaultCluster(cl)

        # pass environment variables to workers

        # parallel::clusterExport(cl,varlist=c(".readParticleTracker","ab.track"),envir=environment())
        #
        # trackll=parallel::parLapply(cl,file.list,function(fname){
        #     track=.readParticleTracker(file=fname,ab.track=ab.track)

        parallel::clusterExport(cl,varlist=c(".readParticleTracker","ab.track", "frameRecord"),envir=environment())

        trackll=parallel::parLapply(cl,file.list,function(fname){
            track=.readParticleTracker(file=fname,ab.track=ab.track, frameRecord = frameRecord)

            # add indexPerTrackll to track name
            indexPerTrackll=1:length(track)
            names(track)=mapply(paste,names(track),indexPerTrackll,sep=".")
            return(track)
        })

        # stop cluster

        cat("\nStopping clusters...\n")

        parallel::stopCluster(cl)

        names(trackll)=file.name

    }

    # cleaning tracks by image mask
    if (mask==T){
        trackll=maskTracks(trackll=trackll,maskl=mask.list)
    }

    # merge masked tracks
    # merge has to be done after mask


    if (merge==T){


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

        # trackll=track.holder
    }
    cat("\nProcess complete.\n")

    return(trackll)
}

# file="/Volumes/SDXC\ Disc/DoScience/Projects/SWR1/DotTracking/Result/2016-10-07/Traj_swc5AA_Swr1Halo_2hRAP_fr10ms_120mW_37-1.tif.csv"
# folder="/Volumes/SDXC\ Disc/DoScience/Projects/SWR1/DotTracking/Data/2016-10-07"



