
## selComponentTracks-methods
##
##
###############################################################################
##' @name selComponentTracks
##' @aliases selComponentTracks
##' @title selComponentTracks
##' @rdname selComponentTracks-methods
##' @docType methods
##'
##' @description Select trajectory based on component fitting on diffusion
##'   coefficient.
##' @usage selComponentTracks(trackll,fit,likelihood=0.9,dcoef,log.transformed=F,output=F)
##' @param fit Component fitting result form fitNormDistr() function.
##' @param likelihood The likelihood of a trajecotry to be in fitted group. This parameter specifies the strigency of selecting trajectories to be in the fitted group and therefore influence the number of trajectories been selected.
##' @param dcoef Diffusion coefficent calcualted by Dcoef, which provide the link between trajecotry index and diffusion coefficent.
##' @param log.transformed A flag indicating if the fitting is been log transformed, select TRUE if fitting was done in fitNormDistr(log.transform = T,..). This parameter will be removed in later version by directly read the info from the output of fitNormalDistr() function.
##' @param output A logical indicating if output of selected trajectory index, which can be used for plot individual trajectory using plotTrack.


##' @return
##' \itemize{
##' \item{combined list of trackll} The result is a combined list of selected trajectories. The list is one level higher than trackll, use subsetting to output trackll.e.g. trackll[[1]], or trackll[["SWR1"]]. }

##' @examples
##'
##' ## selComponentTracks() usage
##' # 1. select componentTracks per folder (cross movie) by using compareFolders
##' # 2. select componentTracks per movie base, use plotComponentTracks to plot
##'  component tracks back to initial Nuclei image.
##'
##' ## 1. select componentTracks per folder (cross movie) by using compareFolders
##' folder1=system.file("extdata","SWR1",package="sojourner")
##' folder2=system.file("extdata","HTZ1",package="sojourner")
##' trackll=createTrackll(folder=c(folder1,folder2),input=1)
##' MSD=msd(trackll=trackll)
##' dcoef=Dcoef(MSD,dt=6,plot=T,output=F)
##' # fit dcoef
##' # for replication purpose set seed to fix number
##' fit=fitNormDistr(dcoef,components=2,log.transform=T,combine.plot=F,output=F,seed=484)
##'
##' # select component tracks from fitting
##' trackll.sel=selComponentTracks(trackll=trackll,fit=fit,likelihood=0.9,
##' dcoef=dcoef,log.transformed=T,output=F)
##' # subset component tracks to further analyze msd, dcoef
##' trackll.swr1=trackll.sel[["SWR1_WT_140mW_image6.txt"]]
##' msd(trackll.swr1,plot=T)
##' msd(trackll.swr1,summarize=T,plot=T)
##' Dcoef(trackll=trackll.swr1,plot=T)
##' plotTrackOverlay(trackll.swr1)
##'
##' # plotNucTrackOverlay(folder=folder1, trackll=trackll.swr1)
##' dwellTime(trackll.swr1)
##'
##' # Output trajectory index to plot individually
##' # trackll.sel=selComponentTracks(trackll=trackll,fit=fit,likelihood = 0.9,
##' dcoef = dcoef,log.transformed = T,output = T)
##' # specify index file path.
##' index.file=system.file("extdata","INDEX","componentTrackID-SWR1.comp.1.csv",package="sojourner")
##' index.file2=system.file("extdata","INDEX","componentTrackID-SWR1.comp.2.csv",package="sojourner")
##' movie.folder=system.file("extdata","SWR1_2",package="sojourner")
##' plotTrackFromIndex(index.file=index.file,movie.folder = movie.folder)
##' plotTrackFromIndex(index.file=index.file2,movie.folder = movie.folder)
##'
##'
##' ## 2. select componentTracks per movie base, use plotComponentTracks to plot 
##' component tracks back to initial Nuclei image.
##' ## plotComponentTrackOverlay
##' folder3=system.file("extdata","SWR1_2",package="sojourner")
##' trackll=readDiatrack(folder3)
##'
##' ## use merge=T for per folder comparison, the analsyis result can't be plot back to original image
##' ## To see component tracks on original nuclei image, set merge=F (for per movie analysis)
##' ## may not make much sense to msd on individual movie, 
##' however for plot component track back to original nuclei image.
##'
##' ## compute MSD
##' MSD=msd(trackll=trackll,plot=T)
##' msd(trackll=trackll,summarize=T,plot=T)
##'
##' ## calculate Dcoef
##' dcoef=Dcoef(MSD=MSD,method="static",plot=TRUE)
##'
##' ## fit normal distribution to define component
##' ## set seed to reproduce results (see fitNormalDistr() for details on seed)
##' fit=fitNormDistr(dcoef,components=2,log.transform=T,combine.plot=F,output=F,seed=481)
##'
##' ## select component tracks based on fitting
##' trackll.sel=selComponentTracks(trackll=trackll,fit=fit,likelihood = 0.9,
##' dcoef = dcoef,log.transformed = T,output = F)
##' ## plot component tracks
##' plotComponentTrackOverlay(folder=folder3,trackll.sel=trackll.sel)
##'
##'



##' @export selComponentTracks
##'

###############################################################################


selComponentTracks=function(
    trackll,
    fit,likelihood=0.9,dcoef,log.transformed=F,output=F){

    # trackll is not used for calculation but for subsetting tracks
    comp=list()
    length(comp)=length(fit)
    names(comp)=names(fit)


    comp.track.index.lst=list()
    length(comp.track.index.lst)=length(fit)
    names(comp.track.index.lst)=names(fit)

    # secondary /internal list structure of comp.track.index.lst is created in a
    # loop below

    # lapply(comp.track.index.lst,function(){})

    #---------------------------------------------------------------------------
    # select tracks from fitNormDistr() outputs
    for (i in 1:length(fit)){

        comp[[i]]=data.frame(fit[[i]]$posterior)

        # specify secondary list structure (length and name)
        # length(comp.track.index.lst[i])=dim(comp[[i]])[2]

        # creates internal list structure of comp.track.index.lst
        comp.track.index.lst[[i]]=list()
        length(comp.track.index.lst[[i]])=dim(comp[[i]])[2]
        names(comp.track.index.lst[[i]])=colnames(comp[[i]])

        for (j in 1:dim(comp[[i]])[2]){

            # select likelihood > 0.9
            comp.track.index.lst[[i]][[j]]=which(comp[[i]][,j]>=likelihood)
        }
    }


    # print selected tracks on the console
    l.print=list()
    length(l.print)=length(comp.track.index.lst)
    names(l.print)=names(comp.track.index.lst)
    for (i in 1:length(comp.track.index.lst)){
        l.print[[i]]=lapply(comp.track.index.lst[[i]],length)
    }
    l.print.result=do.call(rbind,l.print)

    cat("\n","at likelihood of",likelihood,", the number of trajectories selected are:\n")
    print(l.print.result)


    # use comp.track.index.lst to extract track index from dcoef

    comp.trackID.lst=list()
    length(comp.trackID.lst)=length(comp.track.index.lst)
    names(comp.trackID.lst)=names(comp.track.index.lst)


    # log.transformed=fit$log.transformed
    # if log.transformed, remove negative
    # values to make the total index corresponding to..

    if (log.transformed == T){
        dcoef.log.trans=lapply(dcoef,function(x){x=x[(x[,"slope"]>0),]})
        dcoef=dcoef.log.trans
    }

    # dcoef is supposed to have the same structure as comp.track.index.list
    # need a sanity check here TODO

    for (i in 1:length(dcoef)){

        comp.trackID.lst[[i]]=list()
        length(comp.trackID.lst[[i]])=length(comp.track.index.lst[[i]])
        names(comp.trackID.lst[[i]])=names(comp.track.index.lst[[i]])

        for (j in 1:length(comp.track.index.lst[[i]])){

            # subset the trackID using comp.track.index.lst[[i]][[j]]
            comp.trackID.lst[[i]][[j]]= dcoef[[i]][comp.track.index.lst[[i]][[j]],]
        }
    }


    ## (optional) remove zero length list

    # remove empty list out of list
    # Filter(function(k) length(k)>0, mylist)

    # remove na list out of list
    # Filter(function(x) !is.na(x), l)

    comp.trackID.lst=lapply(comp.trackID.lst,function(x){
        Filter(function(k) length(k)>0, x)
    })


    #---------------------------------------------------------------------------
    # export as index file to plot individually

    if (output == T){
        for (i in 1:length(comp.trackID.lst)){
            for (j in 1:length(comp.trackID.lst[[i]])){
                fileName=paste("componentTrackID-",
                               paste(names(comp.trackID.lst)[[i]],".",
                                     names(comp.trackID.lst[[i]][j]),sep=""),
                               .timeStamp("-"),".csv",sep="")
                write.csv(file=fileName,comp.trackID.lst[[i]][[j]])
            }
        }

    }

    #---------------------------------------------------------------------------
    # subset selected component tracks from original trackll for msd, dwell time
    # and further analysis

    # extract trackID

    trackID=list()
    length(trackID)=length(comp.trackID.lst)
    names(trackID)=names(comp.trackID.lst)

    for (i in 1:length(comp.trackID.lst)){
        trackID[[i]]=lapply(comp.trackID.lst[[i]],function(x){rownames(x)})
    }

    componentTracks=list()
    length(componentTracks)=length(trackID)
    names(componentTracks)=names(trackID)


    for (i in 1:length(trackID)){

        ## creates internal list structure of componentTracks
        componentTracks[[i]]=list()
        length(componentTracks[[i]])=length(trackID[[i]])
        names(componentTracks[[i]])=names(trackID[[i]])


        for (j in 1:length(trackID[[i]])){
            componentTracks[[i]][j]=lapply(trackll[i],function(x){
                x[trackID[[i]][[j]]]})
        }

    }

    # remove empty list out of list
    # Filter(function(k) length(k)>0, mylist)

    # remove na list out of list
    # Filter(function(x) !is.na(x), l)

    # remove zero length list
    componentTracks=lapply(componentTracks,function(x){
        Filter(function(k) length(k)>0, x)
    })

    return(componentTracks)
}









