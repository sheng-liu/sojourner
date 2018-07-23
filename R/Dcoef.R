## Dcoef-methods
##
##
###############################################################################
##' @name Dcoef
##' @aliases Dcoef
##' @title Dcoef
##' @rdname Dcoef-methods
##' @docType methods
##' @description Caclulate diffusion coefficient (Dcoef) for trajecotries.
##'
##' @usage
##' Dcoef( MSD=NULL,trackll=NULL,dt=6,filter=c(min=7,max=Inf),rsquare=0.8,resolution=0.107,
##'        binwidth=NULL,method=c("static","percentage","rolling.window"),
##'        plot=F,output=F,t.interval=0.01,profile=NULL)
##' @param MSD Mean Square Displacement calculated using msd() function. Either MSD or trackll can be passed into Dcoef for calculation of diffusion coefficient.
##' @param trackll Track list output from readDiatrack().
##' @param dt Time intervals. Default 6.
##' @param filter a vector specifies the minimum and max length of trajecotries
##'   to be analyzed. Take only trajectories that has number of frames greater
##'   than (>=) min and less than (<) max.
##' @param rsquare rsquare filter on Dcoef results. Default to be 0.8. Set value
##'   to 0 if rsquare filter is not desired.
##' @param resolution ratio of pixel to µM.
##' @param plot A parameter for plotting. Default FALSE, no plot; If TRUE,
##'   automatically plots "histogram" with count information, binwidth can be
##'   set through parameter binwidth; as well as "density" with
##'   density/frequency.
##'
##' @param binwidth binwidth used for histogram. Default NULL, automatically
##'   assign binwidth.
##' @param method "static", uses time lags 2~5 to calculate diffusion
##'   coefficient; "percentage", uses (tierd) percentage (default 0.25) of time
##'   lags (see Details). "rolling.window", time lags uses for Dcoef follows a
##'   rolling window with specified window size (default 4).
##'
##' @param output An Logical indicate if output should be generated. See Values
##'   for detail.
##' @param t.interval time interval between frames, default 0.010 s (10ms).
##' @param profile Location of preference file. By default (NULL), it is stored
##'   at : system.file("extdata","PREF","profile.csv",package="sojourner"). User can
##'   provide preference file by specifying the location of the file, e.g.
##'   profile="/Users/shengliu/Desktop/profile.csv".
##'
##' @return
##' \itemize{
##' \item \emph{Dcoef} A list of Dcoef for each file in trackll.
##' \item \emph{PDF} Log.Dcoef histogram fitted with density curve, when plot = TRUE.
##' \item \emph{csv} Dcoef output in csv format, when output = TRUE.
##' }
##' @details Generic parameters (parameter applied to all methods, such as
##'   resolution etc) are set in the function. Method dependent parameters (such
##'   as lag.start, lag.end for method = "static"), are stored in profile.csv in
##'   PREF folder under extdata. To change preference parameter, can either
##'   programably or manually go to folder
##'   system.file("extdata","PREF","profile.csv",package="sojourner"), and change the
##'   profile.csv.
##'   
##'   lag.start: time lag used as start of dt for compute Dcoef. Default 2.
##'   lag.end: Time lag used as end of dt for compute Dcoef. Default 2.
##'
##'   method for calculating Dcoef:
##' \itemize{
##'     \item \bold{static} stabilize the number of time lags used for fitting
##'     using time lag 2~ 5 despite the total time lags measured.
##'    \item \bold{percentage} "percentage", uses (tierd) percentage (default
##'    0.25) of time lags.
##' \tabular{rlll}{
##'     [,1] \tab TrackLength \tab Percentage \tab TimeLagsForFitting\cr
##'     [,2] \tab 31~ \tab 0.25 \tab 2~5-2~7\cr
##'     [,3] \tab 22~30 \tab 0.25 \tab 2~5-2~7\cr
##'     [,4] \tab 15~21 \tab 0.4 \tab 2~5-2~7\cr
##'     [,5] \tab 10~15 \tab 0.6 \tab 2~5-2~7\cr
##'     [,6] \tab 7~9 \tab 1 \tab 2~5-2~7
##' }
##'
##'    \item \bold{rolling.window}  time lags uses for Dcoef follows a rolling window with specified window size (default 4).
##'}


##'
##' @examples
##' # compare files
##' folder=system.file("extdata","SWR1",package="sojourner")
##' trackll=readDiatrack(folder)
##' MSD=msd(trackll=trackll)
##' Dcoef(MSD=MSD,method="static",plot=TRUE)
##'
##' # compare folders
##' folder1=system.file("extdata","SWR1",package="sojourner")
##' folder2=system.file("extdata","HTZ1",package="sojourner")
##' trackll2=compareFolder(c(folder1,folder2))
##' Dcoef(trackll=trackll2,method="percentage",plot=TRUE)
##' Dcoef(trackll=trackll2,method="rolling.window",plot=TRUE)


##' @export Dcoef
###############################################################################

# Dcoef (Diffusion coefficient)


Dcoef=function(MSD=NULL,trackll=NULL,dt=6,filter=c(min=7,max=Inf),rsquare=0.8,
               resolution=0.107,binwidth=NULL,
               method=c("static","percentage","rolling.window"),
               plot=F,output=F,t.interval=0.01,profile=NULL){


##------------------------------------------------------------------------------
    # validity check for input

    # if neither MSD or trackll is provided
    if (length(MSD) == 0 & length(trackll) == 0){
        stop("\nPlease provide either MSD or trackll")
    }

    method=match.arg(method)

    # if select percentage method but original trackll is not provided
    if (method == "percentage" & length(trackll) == 0){
        stop("\nPlease provide 'trackll' when using percentage method")
    }

    ## set corresponding switches

    ## read in preference parameters
    ## these are some method dependent parameters, the generic parameters (parameter applied to all methods) are set in the function

    ## enable user provided preference file
    if (is.null(profile)){
        profile=system.file("extdata","PREF","profile.csv",package="sojourner")
    }

    PARAM=read.csv(file=profile,header=T,row.names="PARAMETER")
    lag.start=PARAM["lag.start",]
    lag.end=PARAM["lag.end",]
    perc=PARAM["percentage",]
    window.size=PARAM["window.size",]


    # dispatch on "method"
    switch(method,
           static={
               cat("\napplying static,lag.start=",
                   lag.start,"\t","lag.end=", lag.end,"\n")
               static=T
               lag.start=lag.start
               lag.end=lag.end


               # if MSD is not provided
               if (length(MSD) == 0){
                   # calculate MSD
                   MSD=msd(trackll,dt=dt,resolution=resolution,
                           filter=filter,summarize=F)
               }
               # default using MSD if trackll and MSD both present

               # calculate Dcoef using static
               D.coef=Dcoef.static(MSD,lag.start=lag.start,lag.end=lag.end,
                                   t.interval=t.interval)
           },
           rolling.window={


               cat("\napplying rolling window...\n")
               static=F
               window.size=window.size

               # if MSD is not provided
               if (length(MSD) == 0){
                   # calculate MSD
                   MSD=msd(trackll,dt=dt,resolution=resolution,
                           filter=filter,summarize=F)
               }
               # default using MSD if trackll and MSD both present


               # calculate Dcoef using rolling window
               D.coef=Dcoef.roll(MSD,window.size=window.size,t.interval=t.interval)

           },
           percentage={

               static=T

               D.coef=Dcoef.perc(trackll,percentage=perc,weighted=F,
                                 filter=filter, resolution=resolution,
                                 t.interval=t.interval)

           })

#     if (plot == "variance"){
#         ## currently set rollingwindow only for variance plot
#         cat("\nvariance = TRUE, applying rolling window, filter swtiched on\n")
#         rolling.window=T
#         filter=T
#     }else{
#             rolling.window=F
#         }



##------------------------------------------------------------------------------
## call corresponding functions

#     if (rolling.window == T){
#
#         D.coef=Dcoef.roll(MSD,dt=dt)
#         D.coef.subset=rsquare.filter(D.coef,static=F)
#         Log.D.coef=Dcoef.log(D.coef.subset,static=F)
#
#     }else{
#
#         D.coef=Dcoef.static(MSD)
#         D.coef.subset=rsquare.filter(D.coef,static=T)
#         Log.D.coef=Dcoef.log(D.coef.subset,static=T)
#
#     }

    #if (length(rsquare)!=0){
#         D.coef.subset=rsquare.filter(D.coef,rsquare=rsquare,static=static)
#     }else{
#         D.coef.subset=D.coef
#     }

    # further process dispatch on method

    if(method == "static"||method == "percentage"){

        # subset
        D.coef.subset=rsquare.filter(D.coef,rsquare=rsquare)
        D.coef.subset.slope=lapply(D.coef.subset,function(x){x[,"slope"]})

        # logorithm
        Log.D.coef=lapply(D.coef.subset.slope,log10)

    }else if(method == "rolling.window"){
        # subset
        D.coef.subset=rsquare.filter.roll(D.coef,rsquare=rsquare)
        D.coef.subset.slope=lapply(D.coef.subset,function(x){
            for (i in 1:length(x)){
                x[[i]]=x[[i]][,"slope"]
            }
            return(x)
        })

        ## logorithm
        Log.D.coef=lapply(D.coef.subset.slope,function(x){
            for (i in 1:length(x)){
                x[[i]]=log10(x[[i]])
            }
            return(x)
        })

        # simpler than reverse setup
        #                 Log.D.coef=list()
        #                 for (i in 1:length(D.coef.subset.slope)){
        #                     Log.D.coef[[i]]=lapply(D.coef.subset.slope[[i]],log10)
        #                 }
        #                 names(Log.D.coef)=names(D.coef.subset.slope)
        #                 return(Log.D.coef)

    }

##------------------------------------------------------------------------------
## plot

    if (plot == T){

        cat("\nPlotting histogram...\n")
        # see count inforamtion
        histogram=plotHistogram(Log.D.coef,binwidth=binwidth,method=method)

        cat("\nPlotting density...\n")
        # plot frequency so it is easier to compare groups
        density=plotDensity(Log.D.coef,binwidth=binwidth,method=method)

    }

#     plot=match.arg(plot)
#     switch(plot,
#            variance={
#                if (method == "static"||method == "percentage"){
#
#                    cat("\n\nvariance plot for method static and percentage not available for sojourner v0.2 \n\n")
#
# #                    cat("variance plot for method static and percentage does not use rsquare filter. \n")
# #                    Log.D.coef.nofilter=Dcoef.log(D.coef,static=T)
# #                    plotVariance(Log.D.coef.nofilter,method=method)
#
#                }else{ plotVariance(Log.D.coef,method=method)}
#
#               },
#
#            ## needs more work to deal with a list
#
#            # see count inforamtion
#            histogram=plotHistogram(Log.D.coef,binwidth=binwidth,method=method),
#
#            # plot frequency so it is easier to compare groups
#            density=plotDensity(Log.D.coef,binwidth=binwidth,method=method)
#            # else do nothing
#            )


##------------------------------------------------------------------------------
## output

    if (output == T){

        # output csv
        for (i in 1:length(trackll)){
            fileName=paste("Dcoef-",.timeStamp(names(trackll)[i]),".csv",sep="")
            write.csv(file=fileName,D.coef.subset[[i]])
        }


    }

    return(invisible(D.coef.subset))
    # if no subsetting is intended, select rsquare=0
}






