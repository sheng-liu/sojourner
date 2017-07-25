
## readDiatrack-methods
##
##
###############################################################################
##' @name readDiatrack
##' @aliases readDiatrack
##' @title readDiatrack
##' @rdname readDiatrack-methods
##' @docType methods
##'
##' @description read output file (tracks/trajecotries) from Diatrack.
##' @usage
##' readDiatrack(folder, merge = F, ab.track = F, mask = F, cores = 1, frameRecord = T)
##'
##' .readDiatrack(file, interact = F, ab.track = F, frameRecord = F)
##' 
##' maskTracks(folder, trackll)
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
## readDiatrack(folder,merge=F)
## }

##' @examples
##' folder=system.file("extdata","SWR1",package="smt")
##' trackll=readDiatrack(folder)
##' str(trackll,max.level=2)
##'
##' # Not run:
##' # masking with image mask
##' # trackll.masked=readDiatrack(folder=track.folder,merge=F,mask=T)
##' # str(trackll.masked,1)
##'
##' # compare the masking effect
##' # plotTrackOverlay(trackll)
##' # plotTrackOverlay(trackll.masked)
##'
##' # plot mask
##' # mask.list=list.files(path=folder,pattern="_MASK.tif",full.names=T)
##' # plotMask(mask.file=mask.list[[1]])



##' @details
##' default merge = FALSE, so the researcher can assay variations between files. Keep both output as two level list is for simplicity of downstream analysis.
##'
##' Note: the folder name should not contain ".", as it is a key charactero for subsequent indexing of file names.
##'
##'
##' the absolute coordinates trajectory has moved
##'
##' trackID=fileID.frameID.duration.indexPerFile.indexPerTrackll
##'
##' This "indexPerFile" is the index within a diatrackFile.
##'
##' This "indexPerTrackll" is the index within a trackll, which is unique.
##'

##' @export readDiatrack
##' @export .readDiatrack
##' @export maskTracks

###############################################################################

##------------------------------------------------------------------------------
## .readDiatrack
## a function to read one diatrack txt file and returns a list of tracks


.readDiatrack=function(file, interact=F,ab.track=F, frameRecord = F){

    # interactively open window
    if (interact==T) {
        file=file.choose()
    }

    file.name=basename(file)
    cat("\nReading Diatrack file: ",file.name,"...\n")

    ## skip the first 'comment line'
    data=read.table(file=file, header=F, skip=1)

    ## read in frame number line (for future use)
    frame.num=data[1,]

    ## remove frame number line for computation
    data=data[-1,]


    # frame.id
    frame.num.mx=matrix(frame.num,ncol=3,nrow=length(frame.num)/3,byrow=T)
    frame.id=unlist(frame.num.mx[,1])

    ## process the data
    # store coordinates of track in track.list
    track.list=list()
    # store absolute coordinates of track for comparison plots
    ab.track.list=list()
    # store num.tracks.per.file
    num.tracks.per.file=c()

    # select 3 column at a time
    # can use frame number to do this, but this way makes the program more
    # robust with little to non decrease in efficiency

    for (i in 1:(dim(data)[2]/3)) {

        #i=2

        triple=i*3
        track=dplyr::select(data,(triple-3+1):triple)
        colnames(track)=c("x","y","z")
        track=dplyr::filter(track,x!=0,y!=0)

        if (frameRecord){
            track <- cbind(track, "Frame" = c(frame.id[[i]]:(frame.id[[i]]+nrow(track)-1)))
        }

        # the [[]] is important, otherwise only x is included
        track.list[[i]]=track

        # store num.tracks.per.file
        num.tracks.per.file[i]=dim(track)[1]


        ## preprocess to fix coordinates from 0 to max
        ## absolute value of trajectory movement

        abTrack=data.frame(x=track$x-min(track$x),
                            y=track$y-min(track$y))
        ab.track.list[[i]]=abTrack

    }

    ## name the tracks


    # frame.id
    frame.num.mx=matrix(frame.num,ncol=3,nrow=length(frame.num)/3,byrow=T)
    frame.id=unlist(frame.num.mx[,1])

    # duration


    # duration=table(frame.id)
    duration=num.tracks.per.file

    # file.id
    file.subname=substr(file.name,
           start=nchar(file.name)-8,
           stop=nchar(file.name)-4)

    file.id=rep(file.subname,length(duration))

    # indexPerFile
    indexPerFile=seq(from=1,to=length(duration))

    ## trackID=fileID.frameID.duration.indexPerFile
    track.name=paste(file.id,frame.id,duration,indexPerFile,sep=".")

    # name the track
    names(track.list)=track.name
    names(ab.track.list)=track.name


    cat("\n", file.subname, "read and processed.\n")

    if (ab.track==T) return(ab.track.list) else return(track.list)

}

##------------------------------------------------------------------------------
##
# read in mask and derive positive pix form a mask
maskPoint=function(mask.file,plot=F){

    mask.file.name=basename(mask.file)
    # read in tiff mask
    # library(rtiff)
    cat("Reading mask file    ",mask.file.name,"\n")
    mask=rtiff::readTiff(fn=mask.file)
    # plot(mask)

    pospt=which(mask@red!=0,arr.ind=T)
    pos.point=with(data.frame(pospt),data.frame(x=col,y=row))

    # horizontal is the same vertical is fliped as the pixels is counted from
    # upper left in the image, but is counted from lower left in the plot.
    if (plot) {
        plot(mask)
        plot(x=pos.point$x,y=pos.point$y)
        #ggplot(pos.point,aes(x=x,y=y))+geom_point()
    }


    return(pos.point)
}

##------------------------------------------------------------------------------
##
# Use each trajectory's geometric center as unit for clusterization.
# Each data point is an averaged trajectory.

trackCenter=function(trackll){

    # arithmetic mean for geometric center (centroid)
    track.center=list()
    length(track.center)=length(trackll)
    names(track.center)=names(trackll)

    for (i in 1:length(trackll)){
        track.center[[i]]=lapply(trackll[[i]],function(x){
            # round coords
            apply(x,2,function(coord){round(mean(coord))})})
    }

    return(track.center)
}


##------------------------------------------------------------------------------
##
# track.center and pos.point should be one to one cooresponding
posTracks=function(track.center,pos.point){

    # convert list to data.frame
    track.center.df=do.call(rbind.data.frame,track.center)
    names(track.center.df)=c("x","y","z")
    pos.tracks=plyr::match_df(track.center.df,pos.point,on=c("x","y"))
    return(pos.tracks)

}

##------------------------------------------------------------------------------
##

maskTracks=function(folder, trackll){
    
    # read in mask
    maskl=list.files(path=folder,pattern="_MASK.tif",full.names=T)
    
    if (length(maskl)==0){
        cat("No image mask file ending '_MASK.tif' found.\n")
        
    }
    
    mask.track.index=list()
    length(mask.track.index)=length(trackll)
    names(mask.track.index)=names(trackll)
    
    masked.tracks=list()
    length(masked.tracks)=length(trackll)
    names(masked.tracks)=names(trackll)


    for (i in 1:length(trackll)){
        track.center=trackCenter(trackll)[[i]]
        pos.point=maskPoint(maskl[[i]],plot=F)
        mask.track.index[[i]]=posTracks(track.center,pos.point)

        index=rownames(mask.track.index[[i]])

        masked.tracks[[i]]=lapply(trackll[i],function(x){x[index]})[[1]]

    }

    return(masked.tracks)

}


##------------------------------------------------------------------------------
## Note:the list can be named, this wil change the read.distrack.folder 's
## naming no need for naming it

# the mask file has to be named corresponding to its txt file name to work
# correspondingly. as it is read into two list, file.list, and mask.list. there
# is not direct comparison of file name function add in yet in v0.3.4


readDiatrack=function(folder,merge=F,ab.track=F,mask=F,cores=1, frameRecord = T){


    trackll=list()
    track.holder=c()

    # getting a file list of Diatrack files in a directory
    file.list=list.files(path=folder,pattern=".txt",full.names=T)
    file.name=list.files(path=folder,pattern=".txt",full.names=F)
    folder.name=basename(folder)


    # read in tracks
    # list of list of data.frames,
    # first level list of file names and
    # second level list of data.frames

    max.cores=parallel::detectCores(logical=F)

    if (cores==1){

        for (i in 1:length(file.list)){


            track=.readDiatrack(file=file.list[i],ab.track=ab.track, frameRecord = frameRecord)


            # add indexPerTrackll to track name
            indexPerTrackll=1:length(track)
            names(track)=mapply(paste,names(track),indexPerTrackll,sep=".")

            trackll[[i]]=track
            names(trackll)[i]=file.name[i]
        }

    }else{

        # parallel this block of code
        # assign reading in using .readDiatrack to each CPUs

        # detect number of cores
        # FUTURE: if more than one, automatic using multicore

        if (cores>max.cores)
            stop("Number of cores specified is greater than recomended maxium: ",max.cores)

        cat("Initiated parallel execution on", cores, "cores\n")
        # use outfile="" to display result on screen
        cl <- parallel::makeCluster(spec=cores,type="PSOCK",outfile="")
        # register cluster
        parallel::setDefaultCluster(cl)

        # pass environment variables to workers

        parallel::clusterExport(cl,varlist=c(".readDiatrack","ab.track", "frameRecord"),envir=environment())

        # trackll=parallel::parLapply(cl,file.list,function(fname){
        trackll=parallel::parLapply(cl,file.list,function(fname){
            track=.readDiatrack(file=fname,ab.track=ab.track, frameRecord = frameRecord)

            # add indexPerTrackll to track name
            indexPerTrackll=1:length(track)
            names(track)=mapply(paste,names(track),indexPerTrackll,sep=".")
            return(track)
        })

        # stop cluster
        cat("\nStopping clusters...\n")
        parallel::stopCluster(cl)

        names(trackll)=file.name
        # names(track)=file.name

    }

    # cleaning tracks by image mask
    if (mask==T){
        trackll=maskTracks(folder = folder, trackll=trackll)
    }

    # merge masked tracks
    # merge has to be done after mask


    # if (merge==T){
    #     for (i in 1:length(file.list)){
    #         trackll[[i]]=track[[i]]
    #         names(trackll)[i]=file.name[i]
    #     }
    # }


    # trackll naming scheme
    # if merge==F, list takes the name of individual file name within folder
    # file.name > data.frame.name
    # if merge==T, list takes the folder name
    # folder.name > data.frame.name
    
    # merge masked tracks
    # merge has to be done after mask

    if (merge==T){

        # trackll naming scheme
        # if merge==F, list takes the name of individual file name within folder
        # file.name > data.frame.name
        # if merge==T, list takes the folder name
        # folder.name > data.frame.name

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

    #     }else{
    #
    #         # list of list of data.frames,
    #         # first level list of folder names and
    #         # second level list of data.frames
    #
    #         for (i in 1:length(file.list)){
    #
    #             track=.readDiatrack(file=file.list[i],ab.track=ab.track)
    #             # concatenate tracks into one list of data.frames
    #             track.holder=c(track.holder,track)
    #
    #         }
    #
    #         # add indexPerTrackll to track name
    #         indexPerTrackll=1:length(track.holder)
    #
    #         names(track.holder)=mapply(paste,names(track.holder),
    #                                    indexPerTrackll,sep=".")
    #
    #         # make the result a list of list with length 1
    #         trackll[[1]]=track.holder
    #         names(trackll)[[1]]=folder.name
    #

    #
    #
    #         if (mask==T){
    #             trackll=maskTracks(trackll,mask.list)
    #         }
    #
    #     }

    cat("\nProcess complete.\n")
    return(trackll)
}

##-----------------------------------------------------------------------------
## Note:

## if want to keep the names of each data frame come from, use
## if (merge) do.call(c,trackll)




