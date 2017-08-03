
## dwellTime-methods
##
##
###############################################################################
##' @name dwellTime
##' @aliases dwellTime
##' @title dwellTime
##' @rdname dwellTime-methods
##' @docType methods
##' @description Caclulate dwell time (/residence time) for trajecotries.

##' @usage dwellTime(trackll,t.interval=10,x.scale=c(min=0,max=250),plot=T,output=F)
##' @param trackll Track list output from readDiatrack().
##' @param t.interval t.interval time, default = 10ms.
##' @param plot An Logical indicate if plot should be generated. If plot = TRUE, the plot data will also be output.
##' @param output An Logical indicate if output should be generated. 1) dwell time of tracks in the track list output to csv file. Each item in the list will have an individual csv file. 2) Plot PDF and plot data will be saved.


##' @return
##' \itemize{
##' \item{dwell time list} A list of dwell time for every trajectory, separated by file names of the trajectory file in Diatrack file folder. If combined dewell time is intended, use readDiatrack(folder, merge=T) to generate a single length list, then apply this function.
##' \item{PDF} dwell time frequency plot in PDF format, when plot = TRUE.
##' \item{csv} dwell time output in csv format, when output = TRUE.
##' }

##' @examples
##' folder=system.file("extdata","SWR1",package="smt")
##' trackll=readDiatrack(folder)
##' dwellTime(trackll,plot=TRUE)

##' @import reshape2
##' @export dwellTime
###############################################################################


##------------------------------------------------------------------------------
## .dwellTime
## a function to calculate dwell time from a list of data.frame track (trackl). and returns a vector of dwell time.

## nomenclature
## track    data.frame with x,y,z coordinates
## trackl   list of data.frames with x,y,z coordinates, read from one track file
## trackll  list of list of data.frames with x,y,z coordinates, read from multiple track file

.dwellTime=function(trackl,t.interval=10){
    sapply(trackl,function(x){dim(x)[1]*t.interval})
}


dwellTime=function(trackll,t.interval=10,x.scale=c(min=0,max=250),plot=T,output=F){

    ## compute dwell time
    dwell.time=sapply(trackll,function(x){.dwellTime(x,t.interval)})
    file.name=names(trackll)

    ## reshape data for plot
    dwell.time.mlt=reshape2::melt(dwell.time)

    if (length(trackll)==1){

        colnames(dwell.time.mlt)=c("index","variable","value")
    }else{
        colnames(dwell.time.mlt)=c("value","variable")

    }

    histo.plot=ggplot(dwell.time.mlt,
                      aes(x=value,group=variable,fill=variable))+
        geom_histogram(binwidth=t.interval,position="dodge",colour="white")+ ##change from white to red?
        xlim(x.scale["min"],x.scale["max"])+
        theme_bw()+
        theme(legend.title=element_blank())+
        labs(x="Lifetime of trajectories (ms)", y="Number of trajecotries")

    density.plot=ggplot(dwell.time.mlt,
                        aes(x=value,group=variable,col=variable,fill=variable))+
        geom_density(alpha=0.2)+
        xlim(x.scale["min"],x.scale["max"])+
        theme_bw()+
        theme(legend.title=element_blank())+
        labs(x="Lifetime of trajectories (ms)", y="Frequency of trajectories")

    if (plot==T) multiplot(histo.plot,density.plot,cols=1)

    ## output
    if (output==T){

        # output csv
#         for (i in 1:length(trackll)){
#             fileName=paste("Dwell Time-",.timeStamp(file.name[i]),".csv",sep="")
#             write.csv(file=fileName,dwell.time[[i]])
#         }

        # output plot
        tStamp.plotName=paste(.timeStamp(file.name[1]),"...",sep="")
        plotName=paste("Dwell Time Plot-",tStamp.plotName,".pdf",sep="")
        ggsave(filename=plotName,plot=histo.plot,width=8,height=4)

        # output plot data
        plotData=ggplot_build(histo.plot)$data
        plotFile=paste("Dwell Time Plot-",tStamp.plotName,".csv",sep="")
        write.csv(file=plotFile,plotData)
    }
    return(invisible(dwell.time))
}


##-----------------------------------------------------------------------------
##

# freqpoly=ggplot(dwell.time.mlt,aes(x=value,color=variable)) + geom_freqpoly(binwidth=t.interval)+labs(x="Dwell time (ms)", y="Count")+theme_bw()+ theme(legend.title=element_blank())+xlim(0,200)

#     histodensity=ggplot(dwell.time.mlt,aes(x=value,color=variable,fill=variable))+
#         geom_histogram(binwidth=t.interval,position="dodge")+
#         geom_density(aes(y=10*..count..),alpha=0.2)+
#         theme_bw()+
#         theme(legend.title=element_blank())+
#         xlim(0,200)



