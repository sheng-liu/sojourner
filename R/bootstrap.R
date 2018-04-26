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
##' @description Bootstrap confidience intervals with standard errors. bootstrap resamples dataset (e.g. diffusion coefficients) to calculate confidience intervals for a statistic measure of dataset.  

##' @usage
##'   bootstrap(dat,n.rep=10,ind=1,type="ordinary")
##'   
##'   plotBootstrap(d.boot,alpha=1/2)
##'   
##' @param n.rep number of replicates
##' @param ind the index of the column to be analyzed
##' @param dat data to be passed into bootstrap function, it needs to be in a data.frame format
##' @param type the type of resampling, can be "ordinary" (the default), "parametric", "balanced", "permutation", or "antithetic". See boot::boot parameter sim for details. 
##' @param d.boot bootstrapped data
##' @param alpha transparency adjustment, between 0 to 1.
##' @return \itemize{ 
##'
##' \item{Inidvidualbootstrap} returns a matrix of bootstrapped samples, 
##' each row is a bootstrapped replicate of original dataset, and the length of columns are matching the size of original dataset. for example, original data has 40 data points and you want to bootstrap 10 such datasets, you will have a 10 x 40 (row x col) datasets. 
##' }

##' @examples
##' # read in using readDiatrack
##' folder=system.file("extdata","SWR1",package="sojourner")
##' trackll=readDiatrack(folder)
##'
##' # filter track based length
##' trackll.flt=filterTrack(trackll,filter=c(min=5,max=Inf))
##' MSD=msd(trackll.flt,dt=6,summarize=TRUE,plot=TRUE)
##' dcoef=Dcoef(MSD=MSD,method="static",plot=TRUE)
##' 
##' # bootstrap new datasets
##' d.boot=bootstrap(dcoef)
##' 
##' # plot bootstrapped data
##' plotBootstrap(d.boot)

##' @details
##' A wrapper of boot::boot function adapted for data format in sojourner package. 
##' 
###############################################################################


##------------------------------------------------------------------------------
## bootstrap

library(boot)


# file="/Users/shengliu/OneDrive\ -\ Johns\ Hopkins\ University/OneDrive/DoScience/Projects/SWR1/_ParticleTracking/Data/2018-04-26/H2A.Z_Dcoefs.csv"
# dat=read.csv(file=file,header = F)

.bootstrap=function(dat,n.rep=10,ind=1,type="ordinary"){
    
    # simple/ordinary resampling using sample()
    # the benefit is to do some customized manipulation, the down is the ways to resample is limited to ordinary
    
    # dat2=dat[,ind]
    # n.row=length(dat2) # n.row=nrow(dat) #l=length(dd)
    # boot.sample=matrix(data=sample(dat2,size=n.rep*n.row,replace=T),
    #                    nrow=n.rep,ncol=n.row)
    
    
    # boot package allows "ordinary" (the default), "parametric", "balanced",
    # "permutation", or "antithetic"

    # a do-nothing function
    f=function(d,i){
        dd=d[i,] # dat is data.frame/matrix
        # dd=d[i] # when dat is a vector
        # f=mean(dd)
        return(dd)
    }
    
    bt=boot::boot(dat,f,R=n.rep,sim=type)
    boot.sample=bt$t
    
    return(boot.sample)
}


bootstrap=function(d,n.rep=10,ind=1,type="ordinary"){
    d.boot=lapply(d,.bootstrap)
    return(d.boot)
}


plotBootstrap=function(d.boot,alpha=1/2){
    # transpose for plotting
    bs=as.data.frame(t(boot.sample[[1]]))
    bsf=reshape2::melt(bs)
    
    # "alpha impacts the line of stat_ based geoms"
    # https://github.com/tidyverse/ggplot2/issues/1371
    p=ggplot(bsf,aes(x=value,group=variable))+
        geom_line(alpha=alpha,position="identity",stat="density")
    
    # fill with no lines
    # ggplot(bsf,aes(x=value,group=variable,fill=variable))+
    #     geom_density(alpha=1/10,linetype="blank")

    plot(p)
    return(p)
}





