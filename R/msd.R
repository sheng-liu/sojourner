## msd-methods
##
##
###############################################################################
##' @name msd
##' @aliases msd msd.perc msd.track.vecdt
##' @title msd
##' @rdname msd-methods
##' @docType methods
##'
##' @description calculate mean square displacement for individual trajectory or
##'   summarize on trajectories.

##' @usage
##' msd(trackll,dt=6,resolution=0.107,summarize=FALSE,cores=1,
##' plot=FALSE,output=FALSE,filter=c(min=7,max=Inf))
##'   msd.track.vecdt(trackll,vecdt=NULL,resolution=0.107,output=F)
##'   msd.perc(trackll,percentage=0.25,filter=c(min=7,max=Inf),
##' trimmer=c(min=1,max=31),resolution=0.107,output=FALSE)
##'   
##' @param dt Time intervals. Default 6.
##' @param resolution ratio of pixel to uM.
##' @param trackll Track list output from createTrackll().
##' @param summarize An logical indicate if MSD should be calculated on
##'   individual trajectories (Default) or summarized on all trajectories.
##' @param filter a vector specifies the minimum and max length of trajecotries
##'   to be analyzed. Take only trajectories that has number of frames greater
##'   than (>=) min and less than (<) max.
##' @param cores Number of cores used for parallel computation. This can be the
##'   cores on a workstation, or on a cluster. Tip: the computation on each file
##'   will be parallel assigned to each CPU core.
##' @param plot An logical indicate if plot should be generated. See Values for
##'   detail.
##' @param output An logical indicate if output should be generated. See Values
##'   for detail.
##' @param vecdt A list containing varying dt values.
##' @param percentage compute msd based on (tierd) percentage of its total length.
##' @param trimmer vector used for trimming via trimTrack()
##' @return \itemize{ \item{SummarizedMSD} MSD summarized over all trajectories
##' as a function of dt.
##'
##' \item{InidvidualMSD} MSD of individual trajectories at specified dt. Row
##' number corresponding to its dt. Notice only the trajectories that satisfies
##' the specified dt is output, trajectories that does not satisfy (i.e.
##' trajectories satisfies 1:(dt-1)) is not output here.
##'
##' \item{StandardError} Standard Error of the sample mean measures the
##' variations of sample mean to underlying mean, it is estimated as
##' SE=SD/sqrt(N).
##'
##'
##' \item{SampleSize} The sample size (number of tracks/trajectories) used for
##' calculating the msd and standard error.
##'
##' \item{Trackll} The msd function also returns the processed trackll. If passed
##' to a variable, one can then export the trackll with this variable.
##' }


##' @examples
##' # read in using createTrackll()
##' folder=system.file("extdata","SWR1",package="sojourner")
##' trackll=createTrackll(folder=folder, input=3)
##'
##' # filter track based length
##' trackll.flt=filterTrack(trackll,filter=c(min=5,max=Inf))
##' msd=msd(trackll.flt,dt=6,summarize=TRUE,plot=TRUE)
##' str(msd)
##'
##' # focus on a group of trajectory by setting filter greater than dt
##' trackll.flt2=filterTrack(trackll,filter=c(min=7,max=Inf))
##' msd2=msd(trackll.flt2,dt=6,summarize=FALSE,plot=TRUE) # individual
##' msd2=msd(trackll.flt2,dt=6,summarize=TRUE,plot=TRUE) # summarized
##' 

##' @details
##' msd() calculate track (/trajectory)'s mean square displacement as a function
##' of time (dt). For a track of N steps, at each dt, there are N-dt number of
##' sub-trajectory/sub-tracks, mean of dt-wise sub-trajectories/ step-wise sub
##' tracks average subtracks into one number at each dt.
##'
##' the dt number of su-btracks each contains N:N-dt steps. Because minimum step
##' is 1 (N-dt > = 1), so the maxium dt is N-1 (dt < = N-1). As dt increase, the
##' number of steps used to generate that mean decrease with the maxmum dt
##' (dt=N-1) generated from one step.
##'
##' if one wants to focus on a group of trajectory's evolution, he can simply
##' filter on a number that is bigger than the dt he wanted to plot MSD.
##'
##' by assinging cores, computation is paralelled on each each list of trackll
##' (corresponding to one movie file).

###############################################################################


##------------------------------------------------------------------------------
## msd.track

## calculate msd for tracks (data.frame)
## track, data.frame, xyz
## resolution 107nm=1 pixel

# msd() calculate trajectory(/track)'s msd as a function of time (dt). for a
# track that has N steps, at each dt, there are N-dt number of
# sub-trajectory/sub-tracks, mean of dt-wise sub-trajectories/ step-wise sub
# tracks average subtracks into one number at each dt

# the dt number of su-btracks each contains N:N-dt steps
# (N-dt > = 1, so dt < = N-1) because minimum step is 1 (N-dt > = 1),
# so the maxium dt is N-1 (dt < = N-1)
# as dt increase, the number of steps used to generate that mean decrease
# with the maxmum dt (dt=N-1) generated from one step

# compile track wise computations to bytecode for performance
# compiler::enableJIT(3)
# msd.track=compiler::cmpfun(msd.track,options=list(optimize=3))
# compiler::enableJIT(0)

# separate at.dt is to remove repeated calcualtion in msd.trackl

##' @export msd.track
msd.track=function(track,dt=6,resolution=0.107,at.dt=F){

    if (at.dt == F){

        # calculate msd for track at 1~dt
        msd_track=c()
        for (i in seq_len(dt)){

            # caculate msd for track at specified dt
            track.sqd=squareDisp(track,dt=i,resolution=resolution)
            # extract square.disp into a vector
            square.disp=do.call(rbind,track.sqd)$square.disp
            # remove NA and get the genuine mean
            msd_track[i]=mean(square.disp,na.rm=TRUE)
        }

    }else{

        # caculate msd for track at specified dt
        track.sqd=squareDisp(track,dt=dt,resolution=resolution)
        # extract square.disp into a vector
        square.disp=do.call(rbind,track.sqd)$square.disp
        # remove NA and get the genuine mean
        msd_track=mean(square.disp,na.rm=TRUE)
    }

    return(msd_track)
}


# alternative more vectorized implementation
# msd.dt.subtrack=sapply(track.sqd,function(x){
#    mean(x$square.disp,na.rm=TRUE)})
# msd.dt.track[i]=mean(msd.dt.subtrack)

# get the sum of all subtrajectories then average them
#         sum.square.disp=lapply(track.sqd,function(trk){
#           # subsetting list with [["colnames]]
#           sum(trk[["square.disp"]],na.rm=TRUE)})
#         # get the msd
#         msd=mean(do.call(rbind,sum.square.disp))


##-----------------------------------------------------------------------------
## msd.trackl

## calculate msd.track for a list of tracks, trackl

msd.trackl=function(trackl,dt=6,resolution=0.107){

    # extrack track length
    track.len=vapply(trackl,function(x) dim(x)[1], integer(1))

    # stop if dt greater than max track length
    if (dt>(max(track.len)-1)) {
        stop("\nmax track length:\t",max(track.len),"\ndt:\t\t\t",dt,
                "\nTime interval (dt) greater than max track length-1\n")
        }

    msd=list()  # msd list
    msd.individual=list() # msd for individual track at specified dt
    msd.summarized=c() # summarized msd over all tracks
    std.summarized=c() # standard deviation
    num.tracks.sel=c() # store number of tracks used for calculation

    # compute msd for 1:dt
    # note the number of tracks selected at each can vary
    for (i in seq_len(dt)){

        # select tracks longer than i to compute
        trackl.sel=trackl[(track.len-1)>=i]
        num.tracks.sel[i]=length(trackl.sel)

        # stop if no tracks satisfys, otherwise output the number of tracks
        if (num.tracks.sel[i] == 0){
            stop("no track satisfies dt =",i,"\n")
        }else{

            # add "\n" at the end makes it print all, when removed, "\r" it has
            # the effect of animation; add an extra "\n" outside of the loop to
            # make system output in a new line.
            cat("\r",num.tracks.sel[i],"\ttracks length > & =\t",i)
            }

        # calculate msd.track for dt=i
        msd.individual[[i]]=vapply(trackl.sel,function(x){
            msd.track(track=x,dt=i,resolution=resolution,at.dt=TRUE)
            },FUN.VALUE=double(1))

        # calculate summarized msd and se
        N=length(msd.individual[[i]])
        msd.summarized[i]=mean(msd.individual[[i]],na.rm=TRUE)
        std.summarized[i]=sd(msd.individual[[i]],na.rm=TRUE)/sqrt(N)

    }

    # add an extra "\n" outside of the loop to make system output in a new line.
    cat("\n")


    # output

    # collapse /binding  msd.individual list into matrix cooresponding to their
    # name for plotting. When reshape2::melt for plotting, reshape2::melt function
    # generates an index when pass in as matrix but not data.frame. format
    # msd.individual from vector to single row matrix. this is in line with tidy
    # data format where variables are columns and observations are rows.

    # format msd.individual from vector to single row matrix.
    msd.individual.mx=lapply(msd.individual,function(x){
        t(as.matrix(x))
        })

    # collapse /binding matrix of list with colnames
    msd.individual=do.call(plyr::rbind.fill.matrix,msd.individual.mx)

    # format msd.summarized into tidy data format it is so much easier to do the
    # formating here than down the line when need processing
    msd.summarized=cbind(SummarizedMSD=msd.summarized,
                        StandardError=std.summarized,
                        NumTracksAtDt=num.tracks.sel)

    msd$IndividualMSD=msd.individual
    msd$SummarizedMSD=msd.summarized

    return(msd)

}

# old version corrections
# when sapply(...,SIMPLIFY=TRUE)
# i=1, sapply returns a vecotr
# i>1, sapply returns a matrix

# # convert result into unified matrix format
# if (is.matrix(msd.individual)!=TRUE){
#     msd.individual=t(as.matrix(msd.individual))}
#
# # compute msd.summarized and std.summarized
# N=length(msd.individual[i,][!is.na(msd.individual[i,])])
# msd.summarized[i]=mean(msd.individual[i,],na.rm=TRUE)
# std.summarized[i]=sd(msd.individual[i,],na.rm=TRUE)/sqrt(N)


# output
#     msd$InidvidualMSD=msd.individual
#     msd$SummarizedMSD=msd.summarized
#     msd$StandardError=std.summarized
#     msd$NumTracksAtDt=num.tracks.sel

##-----------------------------------------------------------------------------
## msd.trackll

# This function applies msd.trackl() to a folder of list, trackll. If cores are
# more than one, it assigns lists to CPU cores to compute msd.trackl() in
# parallel.

# paralel track-wise computation, rather give it to CUDA (thousands); parallel
# trackl computation to CPUs (tens)

msd.trackll=function(trackll,dt=6,resolution=0.107,cores=1){

    # detect number of cores
    # FUTURE: if cores more than one, automatic using multicore
    max.cores=parallel::detectCores(logical=TRUE)

    if (cores == 1){

        msd_trackll=lapply(trackll,function(x){
            msd=msd.trackl(trackl=x,dt=dt,resolution=resolution)
            cat("\n...\n") # a seperator to make output clearer
            return(msd)
        })

        }else{

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
                                    varlist=c("msd.trackl","dt","resolution"),
                                    envir=environment())

            msd_trackll=parallel::parLapply(cl,trackll,function(x){
                msd=msd.trackl(trackl=x,dt=dt,resolution=resolution)
                cat("\n...\n") # a seperator to make output clearer
                return(msd)
            })

            # stop cluster
            cat("Stop clusters...\n")
            parallel::stopCluster(cl)

    return(msd_trackll)
        }

}

# when testing, need to pass all variables
# parallel::clusterExport(cl,varlist=c(
#     "squareDisp","msd.track","msd.trackl","dt","resolution"),
# envir=environment())


##------------------------------------------------------------------------------
## msd

# computes msd for single trajectory as well as summarized msd for all
# trajectories over a specified time lags (τ) of 1, 2, 3, … ,τ time steps (dt)
# in the system. It calls msd.trackll() function (can be paralleled), and does
# plotting and output msd csv files.


##' @export msd
##'
msd=function(trackll,dt=6,resolution=0.107,summarize=FALSE,cores=1,plot=FALSE,
             output=FALSE,
             filter=c(min=7,max=Inf)){

    ## keep this code here for backward compatability
    ## first remove documentation, code will be remove in next version
    ## filterTrack
    trackll=filterTrack(trackll,filter=filter)


    # compute MSD using msd.trackll() function
    MSD=msd.trackll(trackll,dt=dt,resolution=resolution,cores=cores)
    file.name=names(trackll)

    # if file.name has ".", replace it with "_", as it interference with Index
    # indentifier "."
    file.name=gsub('\\.', '_', file.name)

    # if summarize
    if (summarize == TRUE){

        # remove the IndividualMSD
        MSD.summarized=lapply(MSD,function(x){x$SummarizedMSD})

        # if file.name has ".", replace it with "_", as it interference with Index
        # indentifier "."
        names(MSD.summarized)=gsub('\\.', '_', names(MSD.summarized))


        # do.call(rbind,MSD.summarized) lost rownames
        n= do.call(rbind.data.frame,MSD.summarized)

        # extract Index from rownames

        # readParticleTracker introduces .csv into track.name
        # readDiatrack introduces .txt into track.name

        # using file name as part of index impose vonerability to the program,
        # unless rules of naming a file is strictly followed, problem may arise
        # sporadically. Rules: "." is conserved, do not use it in file name.

        # this time file name has ".", stopped the program
        # it is better to remove the file name part or do the exchange at the
        # begining of the program

        if (length(grep("txt",rownames(n)[1]))!=0){
            # if track.name has .txt
            Index=strsplit(rownames(n),".txt.")
        }else if (length(grep("csv",rownames(n)[1]))!=0){
            # if track.name has .csv
            Index=strsplit(rownames(n),".csv.")
        }else {
            # rest, track name doesn't have either .txt or .csv
            Index=strsplit(rownames(n),split="[.]")  # split="\\."
        }

        Index=do.call(rbind,Index)
        colnames(Index)=c("file.name","dt")

        rownames(n)=NULL
        p=cbind(n,Index)

        # tested alternative way using melt
        # m=reshape2::melt(m)
        # MSD.summarized has three variables, and each needs to be treated
        # differently, use data.frame itself will fit the purpuse; unlike
        # MSD.individual only has one variable (with 6 rows), use
        # reshape2::melt() would be appropriet.

        # tested removing last .txt or .csv in tracklist name during
        # readDiatrack and treadParticleTracker, the full name maybe used in
        # passing input file name into functions, for now maybe it is good to
        # keep the full name and as the input file type is only two, that
        # solves.

        # change dt from factor to integer/numeric
        # alternative as.numeric(levels(x))[x]
        p=transform(p,dt=as.integer(as.character(dt)))
        yMin=p$SummarizedMSD-p$StandardError
        yMax=p$SummarizedMSD+p$StandardError
        msd.plot=ggplot(

            # p,aes(x=as.integer(as.character(dt)), not work
            p,aes_string(x="dt",y="SummarizedMSD",group="file.name",
                         col="file.name"))+
            geom_line()+geom_point()+
            geom_errorbar(aes_string(ymin="yMin",ymax="yMax"), width=.1)+
            # this makes integer breaks
            scale_x_continuous(breaks=scales::pretty_breaks())+
            labs(x="Time intervals", y="SummarizedMSD (um^2)")+
            theme_bw()+
            theme(legend.title=element_blank())

        if (plot == TRUE) plot(msd.plot)

        if (output == TRUE){
            fileName=paste("MSD Summarized-",
                            .timeStamp(file.name[1]),"___.csv",sep="")
            cat("\nOutput MSD for summarized trajectories.\n")
            write.csv(file=fileName,p)
        }

        return(invisible(MSD.summarized))

    }else{
        # calculate msd for individual MSD at 1:dt and output without
        # summarize,this is useful for diffusion coefficient calculation

        # extract only IndividualMSD from list
        MSD.individual=lapply(MSD,function(x){x$IndividualMSD})

        if(plot == TRUE){

            # MSD.individual only has one variable (with 6 rows), use
            # reshape2::melt()

            p=reshape2::melt(MSD.individual,level=1,na.rm=TRUE)
            colnames(p)=c("index","track.name","msd","file.name")

            # note group needs to be on two variable (i.e. track.name and
            # file.name) as track.name starts over again from a new file, use
            # interaction() realize it

            
            p$inter=interaction(p$file.name,p$track.name)
            
            msd.plot=ggplot(p,aes_string(x="index",y="msd",
                                # this conventional way somehow does not work
                                # group=interaction("file.name","track.name"),
                                
                                    group="inter",
                                    col="file.name"))+
                geom_line()+
                # this makes integer breaks
                scale_x_continuous(breaks=scales::pretty_breaks())+
                labs(x="Time intervals", y="MSD (um^2)")+
                theme_bw()+
                theme(legend.title=element_blank())

            plot(msd.plot)
        }

        # output csv, each list a seperate file
        if (output == TRUE){
            # reformat output
            MSD.individual.output=lapply(MSD.individual,function(x){
                x=reshape2::melt(x,na.rm=TRUE)
                colnames(x)=c("dt","track.name","msd")
                return(x)
            })

            for (i in seq_along(MSD.individual.output)){
                fileName=paste("MSD individual-",
                                .timeStamp(file.name[i]),".csv",sep="")
                cat("\nOutput MSD for individual trajectories.\n")
                write.csv(file=fileName,MSD.individual.output[[i]])
            }
        }

    return(invisible(MSD.individual))
    }

}

##------------------------------------------------------------------------------
## msd.track.vecdt

## This function is a special case of msd.track(), where dt is of different
## values, corresponding to the first 25% of the total length, when calculate
## MSD for each track.

## It is specifically designed for diffusion coefficient calculation using
## percentage method Dcoef.perc(). It uses first 25% of MSD-dt plot positions
## for fitting and derivation of the coefficient.

## For all trajectories, because their length (N) can vary, the number of time
## intervals (dt) used for analysis (N/4) also varies; however, for each
## trajectory, dt is a fixed number corresponding to 1/4 of its length.


##' @export msd.track.vecdt
msd.track.vecdt=function(trackll,vecdt=NULL,resolution=0.107,output=F){

    # copy trackll's structure
    msd.lst=list()
    length(msd.lst)=length(trackll)
    names(msd.lst)=names(trackll)


    track.name=lapply(vecdt,names)

    # i folder name level
    for (i in seq_along(trackll)){

        # j data.frame level
        for (j in seq_along(trackll[[i]])){

            # add "\n" at the end makes it print all, when removed, "\r" it has
            # the effect of animation; add an extra "\n" outside of the loop to
            # make system output in a new line.
            cat("\rcalculating MSD for individual tracks...","folder ",i,
                " track ",j)
            msd.lst[[i]][[j]]=as.matrix(msd.track(track=trackll[[i]][[j]],
                                                  dt=vecdt[[i]][[j]],
                                                  resolution=resolution))

            # use matrix as output to be uniform with msd.inividual output
            # names(vecdt[[i]][j]) subsets to list; names(msd.lst[[i]][[j]])
            # subsets to element
            colnames(msd.lst[[i]][[j]])=names(vecdt[[i]][j])
        }
    }

    # add an extra "\n" outside of the loop to make system output in a new line.
    cat("\n")

    if(output == TRUE){

        p=reshape2::melt(msd.lst)
        colnames(p)=c("frame.index","track.name","msd","track.number",
                      "file.name")

        fileName=paste("MSD individual-",
                        .timeStamp("vecdt"),".csv",sep="")
        cat("\nOutput MSD for individual trajectories.\n")
        write.csv(file=fileName,p)
    }
    return(msd.lst)
}


## dt is in a list, track is in a list, multiple variables in the operation, 
## use for loop's i j system maybe better. mapply, similar as all other apply
## functions, works only for functions that require preferentially for one
## parameter. e.g. mapply(rep, 1:4, 4:1), repeat 1 4 times, 2, 3times, etc.

# mapply(msd.track,trackll,n,MoreArgs=list(resolution=0.107))
# mapply(msd.track,vectdt=n[[1]],trackll=trackll[[1]])

##------------------------------------------------------------------------------
## msd.track.vecdt

# compute msd based on (tierd) percentage of its total length, rather than
# specified dt. This is specially designed for coefficient calculation using
# percentage method.

##' @export msd.perc
msd.perc=function(trackll,percentage=0.25,filter=c(min=7,max=Inf),
                    trimmer=c(min=1,max=31),resolution=0.107,output=FALSE){


    cat("\napplying percentage,",percentage,"\n")

    # filter off tracks
    trackll=filterTrack(trackll,filter=filter)
    # trimm off tracks
    trackll=trimTrack(trackll,trimmer=trimmer)

    # determine the length of each trajectory N then compute first 25% N's
    # msd, manipulate N before hand, then pass in the N vector to dt
    N=list()

    for (i in seq_along(trackll)){
        N[[i]]=vapply(trackll[[i]],function(x){dim(x)[1]},integer(1))
    }
    names(N)=names(trackll)

    # n, number of frames used for MSD calculation
    n=N
    l=vapply(n,function(x) length(x), integer(1))


    # for loop, use i to traversing through list, j traversing though vector
    for (i in seq_along(n)){

        for (j in seq_len(l[i])) {

            # remove last frame to derive dt
            # remove last dt to be more accurate in slope calculation
            if (n[[i]][[j]]>=22) {
                n[[i]][[j]]=round(percentage*n[[i]][[j]]-1-1)

            }else if (21>=n[[i]][[j]]&n[[i]][[j]]>=15) {
                n[[i]][[j]]=round(0.4*n[[i]][[j]]-1-1)

            }else if (14>=n[[i]][[j]]&n[[i]][[j]]>=10) {
                n[[i]][[j]]=round(0.6*n[[i]][[j]]-1-1)

            }else if (9>=n[[i]][[j]]&n[[i]][[j]]>=7) {
                n[[i]][[j]]=round(1*n[[i]][[j]]-1-1)

                # remove last frame to derive dt
                # do not remove last dt
            }else if (6>=n[[i]][[j]]) {
                n[[i]][[j]]=round(1*n[[i]][[j]]-1)
            }
        }
    }

    # lapply(n,summary);lapply(N,summary)

    # calculate msd
    msd.lst=msd.track.vecdt(trackll,vecdt=n,resolution=resolution,output=output)

    return(msd.lst)


}

