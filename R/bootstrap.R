## bootstrap-methods
##
##
###############################################################################
##' @name bootstrap
##' @aliases bootstrap
##' @title bootstrap
##' @rdname bootstrap-methods
##' @docType methods
##'
##' @description 
##' Bootstrap confidience intervals with standard errors. bootstrap resamples
##' dataset (e.g. diffusion coefficients) to calculate confidience intervals for
##' a statistic measure of dataset.

##' @usage
##'   bootstrap(fittedObj,n.reps=100)
##'   
##   plotBootstrap(d.boot,alpha=1/2) -- not available in current version
##'   
##' @param fittedObj output from fitNormDistr
##' @param n.reps number of replicates for bootstrapping
## @param d.boot bootstrapped data, or the output from bootstrap
## @param alpha transparency adjustment, between 0 to 1.
##' @return \itemize{ 
##' List of length 2 containing:
##' \item{Fit}
##' results from fitNormDistr that was used as the input for the function
##' \item{Bootstraps}
##' results from bootstrapping. mean values of each samples,
##' standard error for the mean, lambda(proportion)
##'  }

##' @examples
##' # read in using readDiatrack
##' folder=system.file("extdata","SWR1",package="sojourner")
##' trackll=readDiatrack(folder)
##'
##' # filter track based length
##' trackll.flt=filterTrack(trackll,filter=c(min=5,max=Inf))
##' MSD=msd(trackll.flt,dt=6,summarize=FALSE,plot=TRUE)
##' dcoef=Dcoef(MSD=MSD,method="static",plot=FALSE)
##' # fit the dcoef result
##' fittedObj=fitNormDistr(dcoef)
##' 
##' # bootstrap new datasets
##' d.boot=bootstrap(fittedObj)
##' # manually set the number of bootstrap samples to 50
##' d.boot=bootstrap(fittedObj, n.reps=50)
##' 

##' @details
##' A wrapper of boot::boot and mixtools::boot.se adapted for data format in
##' sojourner package. Also returns stderr information by running boot.se from
##' mixtools and an additional method for one-component distribution calculates
##' the stderr separately. For multi-component distributions, the boot.se
##' function from the mixtools package was used. For single-component
##' distributions, a separate function was used to calculate the stderr and
##' confidence interval.
##' 
##' @import boot
###############################################################################

##------------------------------------------------------------------------------
## bootstrap

##bootstrap and then return stderr for means 
##Ref: http://www.stat.wisc.edu/~larget/stat302/chap3.pdf
.boot.se.onecomp=function(fitResult, n.reps=100){
    data=fitResult$x
    #for choosing in bootstraps
    f=function(d,i){d[i]}
    bootstraps = boot::boot(data,f,n.reps)
    means=apply(bootstraps$t, 1, mean)
    list(mu.se=sd(means), lambda.se=0)
}

# basically a wrapper for stderr calculation
.boot.se.all=function(fitResult, B){
    if(length(fitResult$mu) == 1){
        boot.se.result = .boot.se.onecomp(fitResult,B)
    } else {
        capture.output({multi.stderr=mixtools::boot.se(fitResult, B)})
        boot.se.result = multi.stderr
    }
    return(c(boot.se.result, n.reps=B))
}

##' @export bootstrap
bootstrap=function(fittedObj, n.reps=100){
    d.boot=lapply(fittedObj,.boot.se.all, B=n.reps)
    return(list(Fit=fittedObj, Bootstraps=d.boot))
}

## not ready yet
plotBootstrap=function(d.boot,alpha=1/2){
    # transpose for plotting
    sample_names=names(d.boot)
    for (i in 1:length(d.boot)){
        bs=as.data.frame(t(d.boot[[i]]))
        bsf=reshape2::melt(bs)
    
    # "alpha impacts the line of stat_ based geoms"
    # https://github.com/tidyverse/ggplot2/issues/1371
        p=ggplot(bsf,aes_string(x="value",group="variable"))+
            geom_line(alpha=alpha,position="identity",stat="density")+
            ggtitle(sample_names[i])+
            theme(plot.title = element_text(hjust = 0.5))
        plot(p)
    }
}



# TODO: isolate the boot.se part from fitNormDistr(), so the bootstrap is all in
# one place.


# Note the calculation of means for the bootstrapped sample seems to be directly
# using fitted posterior distributions, rather than fitted 100 times.






