
## displacementCDF-methods
##
##
###############################################################################
##' @name displacementCDF
##' @aliases displacementCDF
##' @title displacementCDF
##' @rdname displacementCDF-methods
##' @docType methods
##'
##' @description calculate cumulative distribution function of all displacement for individual trajectories.

##' @usage displacementCDF(trackll,dt=6,resolution=0.107,plot=F,output=F)
##' @param dt Time intervals.
##' @param resolution ratio of pixel to µM.
##' @param trackll Track list output from readDiatrack().
##' @param plot An logical indicate if plot should be generated. See Values for
##'   detail.
##' @param output An logical indicate if output should be generated. See Values
##'   for detail.
##' @details The cumulative radial distribution function, P(r, i△t), is the probability of finding the diffusing particle within a radius r from the origin at time lag i△t:
##'
##' P(r,iΔt) = 1− e^(-r^2/4*D*(iΔt))
##'
##' the CDF and UniqueDisplacement in the output file is corresponding to P and r in this formula. If intend to generate the CDF plot from the output file, the CDF and UniqueDisplacement is corresponding to the y and x values in the CDF output plot.

##' @return
##' \itemize{
##' \item{ list of "stepwise.displacement" and "CDF.displacement",} {A list of stepwise.displacement" and "CDF.displacement". the name of the list is the Diatrack folder name.}
##'
##' \item{Output file,} {Displacement of individual trajectoreis at specified dt.}
##'
##'
##' \item{CDF plot,} {CDF plot of displacement for individual files.
##' }
##' }


##' @examples
##' folder1=system.file("extdata","SWR1",package="smt")
##' folder2=system.file("extdata","HTZ1",package="smt")
##' trackll=compareFolder(c(folder1,folder2))
##' displacementCDF(trackll,dt=1,plot=T)

##' @export displacementCDF

# do a cdf on single steps ( make step adjustable). already have all the steps, just need a CDF function. for each dt. better have a table of dt as first column, sqrt(square disp) as rows. then one can easily call CDF on this table's rows.


# it is of different length at different dt,
# for each dt for all trajectories, can use a matrix/data.frame as it is of same length, and dt number of such matrix/data.frame


## calculate displacement for tracks (data.frame)
## track, data.frame, xyz
## resolution 107nm=1 pixel

# focus only on displacement, however a one step displacement is also a track, for the consistancy in naming, also called displacement.track

##------------------------------------------------------------------------------
## displacement.track
## displacement for a single data.frame
##' @export displacement.track
displacement.track=function(track,dt=6,resolution=0.107,bivar=F){

    # validity check for dt less than track length
    if (dt >(dim(track)[1]-1)){
        stop("\ntrack length:\t",dim(track)[1],
             "\ndt:\t\t",dt,
             "\nTime interval (dt) greater than track length-1\n")
    }
    # instead of stop, calculate the maximum steps
    # cat("Time interval (dt) greater than track length-1, calculating the maximum time interval")

    # summarize displacement for track at all dt
    # note this function calculates only "at" all dt
    # displacement.dt.track=list()
    # for (i in 1:dt){

        # at each dt, there are dt number of sub-trajectory/sub-tracks
        # displacement of dt-wise sub-trajectories/ step-wise sub tracks

        # caculate displacement for track at specified dt

        track.sqd=squareDisp(track,dt=dt,resolution=resolution)
        track.disp=lapply(track.sqd,function(x){
            x["square.disp"]=sqrt(x["square.disp"])

            # specifically change the modified column name
            colnames(x)[which(colnames(x)=="square.disp")]="displacement"
            return(x)
        })

        # pull all the displacement at this dt together
        # displacement=do.call(rbind,track.disp)$displacement

        if(bivar==T){

            # get the displacement at that dt
            displacement=track.disp[[dt]]["displacement"]
            displacement=displacement[!is.na(displacement)]

            # at each dt, there is length(track)-dt number of displacement
            # for an 11 step track, dt=1, 10 displacment
            # dt=2, 9 displacement
            # dt=3, 8 displacement
            # ...
            # dt=10, 1 displacement

            # get the bivariate dx & dy at that dt
            bivariate=cbind(track.disp[[dt]]["dx"],track.disp[[dt]]["dy"])
            bivariate=bivariate[complete.cases(bivariate),]

            return(bivariate)

        }
        # put the displacement into a list as it has different length
        # displacement.dt.track[[i]]=displacement

        # to be compatible, get the displacement into a list
        displacement.dt.track=lapply(track.disp,function(x){
            d=x$displacement
            d=d[!is.na(d)]})
    #}

    return(displacement.dt.track)
}

##------------------------------------------------------------------------------
## displacement.trackl
## calculate displacement for a list of data.frame

## calculate displacement.track for trackl (list of data.frame) one level
displacement.trackl=function(trackl,dt=6,resolution=0.107,bivar=F){

    # validity check for max track length greater than dt
    track.len=sapply(trackl,function(x) dim(x)[1])
    if (dt>(max(track.len)-1)) {
        stop("\nmax track length:\t",max(track.len),
             "\ndt:\t\t\t",dt,
             "\nTime interval (dt) greater than max track length-1\n")
    }

    # subset trackl for each dt
    num.tracks=c()

    displacement.summarized=list()
    #length(displacement.summarized)=length(displacement.individual)
    #names(displacement.summarized)=names(displacement.individual)

    std.summarized=list()
    #length(std.summarized)=length(displacement.individual)
    #names(std.summarized)=names(displacement.individual)


    displacement=list()

    for (i in 1:dt){

        # select tracks longer than i (i.e. "satisfy" dt=i)
        trackl.dt=trackl[(track.len-1)>=i]

        # double check if tracks exist, calculate displacement for trackl
        if (length(trackl.dt)==0){
            stop("no track satisfies dt =",i,"\n")
        }

        cat("\n",length(trackl.dt),"tracks satisfy dt =",i,"\n")

        num.tracks[i]=length(trackl.dt)
        displacement.individual=sapply(trackl.dt,function(x){
            displacement.track(track=x,dt=i,resolution=resolution,bivar=bivar)},simplify = F)

        # as the result is of different length, the output is a list
    }




        # calculate mean and sd for displacement at each dt
        N=length(displacement.individual)
        for ( i in 1:N){
            displacement.summarized[i]=list(sapply(
                displacement.individual[[i]],mean,na.rm=T,simplify=F))

            std.summarized[i]=list(sapply(
                #displacement.individual[[i]],function(x){sd(x)/N},na.rm=T,simplify=F))
            displacement.individual[[i]],sd,na.rm=T,simplify=F))
        }

        names(displacement.summarized)=names(displacement.individual)
        names(std.summarized)=names(displacement.individual)


    # output
    displacement$InidvidualDisplacement=displacement.individual
    displacement$SummarizedDisplacement=displacement.summarized
    displacement$StandardError=std.summarized
    displacement$NumTracksAtDt=num.tracks

    return(displacement)

}

##------------------------------------------------------------------------------
## displacement.trackll

displacement.trackll=function(trackll,dt=6,resolution=0.107,bivar=F){


    displacement.trackll.lst=lapply(trackll,function(x){
        displacement=displacement.trackl(trackl=x,dt=dt,resolution=resolution,bivar=bivar)
        cat("\n...\n") # a seperator to make output clearer

        return(displacement)
    })

    return(displacement.trackll.lst)
}

##------------------------------------------------------------------------------
## plot CDF of individual displacement

# the minimum tracks to include for the dt is recommended as 50. Users can see the screen output, or NumTracksAtDt of displacement.trackll to decide the dt that suits.

displacementCDF=function(trackll,dt=6,resolution=0.107,plot=F,output=F,bivar=F){
    dp=displacement.trackll(trackll,dt=dt,resolution=resolution,bivar=bivar)

    # take "InidvidualDisplacement" out
    InidvidualDisplacement=list()
    length(InidvidualDisplacement)=length(dp)
    names(InidvidualDisplacement)=names(dp)
    for (i in 1:length(dp)){

        InidvidualDisplacement[[i]]=dp[[i]]["InidvidualDisplacement"]
    }

    # collapse all dp at dt
    dp.dt=list()
    length(dp.dt)=length(dp)
    names(dp.dt)=names(dp)
    for ( i in 1:length(dp)){
        # InidvidualDisplacement[[i]][[1]] # the [[1]] is to move it one level
        # off "IndividualDisplacement"
        dp.dt[[i]]=lapply(InidvidualDisplacement[[i]][[1]],function(x){
            x[dt]
        })

    }

    # reshape for plotting and output
    p=melt(dp.dt)

    ## plotting

        ecdf=ggplot(p,aes(x=value,group=L1,colour=L1))+stat_ecdf()+
            labs(x="Displacement (µm)",y="CDF")+
            theme_bw()+
            theme(legend.title=element_blank())
        # can use stat_ecdf(pad=F) to remove first -Inf and Inf dded on x in ggplot2::stat_cdf, however it seems not working in current version ggplot2 2.1.10, remove it manually in "preprocessing of data for output"

        histogram=ggplot(p,aes(x=value,group=L1,colour=L1))+geom_density()+
            labs(x="Displacement (µm)",y="Density")+
            theme_bw()+
            theme(legend.title=element_blank())
        if (plot==T){
        multiplot(ecdf,histogram,cols=1)
    }

    ## output
    file.name=names(trackll)
    stepwise.displacement=lapply(dp.dt,melt)

    # export the data within ecdf
    # only x, y is relevant
    # the CDF and UniqueDisplacement is corresponding to the y and x
    data=ggplot_build(ecdf)$data[[1]]

    ## preprocess the data before output

    dat=with(data,data.frame(y,x,group))

    # remove Inf
    dat=dat[!is.infinite(dat$x),]
    # make it a list
    CDF.displacement=split(dat,f=dat$group)
    # name the list
    # the order of the grouping in ggplot2 is sorted by factor,
    # so group 1 is levels(factor(c("SWR1","HTZ1")))[1]
    #     > levels(factor(c("SWR1","HTZ1")))
    #     [1] "HTZ1" "SWR1"
    name.CDF.displacement=levels(factor(names(stepwise.displacement)))
    names(CDF.displacement)=name.CDF.displacement

    # format output by removing group column and renaming x,y column
    CDF.displacement=lapply(CDF.displacement,function(x) {

        x$group=NULL
        names(x)=c("CDF","UniqueDisplacement")
        return(x)})



    for (i in 1:length(stepwise.displacement)){
        stepwise.displacement[[i]]["L2"]=NULL
        colnames(stepwise.displacement[[i]])=c("stepwiseDisplacement","trackIndex")
    }

    if (output==T){

        for (i in 1:length(stepwise.displacement)){
            fileName=paste("stepwiseDisplacement-",
                           .timeStamp(file.name[i]),"....csv",sep="")
            cat("\nOutput stepwiseDisplacement for",file.name[i],"\n")
            write.csv(file=fileName,stepwise.displacement[[i]],row.names = F)
        }

        for (i in 1:length(CDF.displacement)){
            fileName=paste("CDFDisplacement-",
                           .timeStamp(file.name[i]),"....csv",sep="")
            cat("\nOutput CDFDisplacement for",file.name[i],"\n")
            write.csv(file=fileName,CDF.displacement[[i]],row.names = F)
        }

    }

    output.lst=list(stepwise.displacement,CDF.displacement)
    names(output.lst)=c("stepwise.displacement","CDF.displacement")


    return(output.lst)

}

## -----

# DONE: output dat V
# DONE: base:ecdf 1228 seems to be more accurate as ggplot2::stat_ecdf also have two outside values, one is negative at zero, the other is extra 1 at the end, which should be removed.  V
