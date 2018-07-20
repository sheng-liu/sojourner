## fitNormDistr-methods
##
##
###############################################################################
##' @name fitNormDistr
##' @aliases fitNormDistr
##' @title fitNormDistr
##' @rdname fitNormDistr-methods
##' @docType methods
##' @description fit normal distributions to diffusion coefficient caclulated by Dcoef method.
##'
##' @usage
##' fitNormDistr(dcoef,components=NULL,log.transform=F,binwidth=NULL,combine.plot=F,output=F,
##'              seed=NULL,meanRange=NULL, proportion=NULL, means=NULL, sd=NULL)
##' @param dcoef diffusion coefficient calculated from Dcoef().
##' @param components parameter specifying the number of components to fit. If NULL (default), a components analysis is done to determine the most likely components and this number is then used for subsequent analysis.
##' @param log.transform logical indicate if log10 transformation is needed, default F.
##' @param binwidth binwidth for the combined plot. If NULL (default), will automatic assign binwidth.
##' @param combine.plot Logical indicate if all the plot should be combined into one, with same scale (/same axises breaks), same color theme, and same bin size for comparison.
##' @param output Logical indicate if output file should be generated.
##' @param seed Seed for random number generator. This makes each run easily repeatable. Seed will be automatically assigned if no seed is specified (default). The seed information is stored as an attribute of the returned object. The seed can also be output to a txt file when output=T.
##' @param meanRange list of vectors, each of which have the min&max that estimates a range where a peak of the distribution should lie.
##' @param proportion numeric vector with estimates of each component's proportion of the whole data.
##' @param means numeric vector with estimates of mean(mu) values for each component.
##' @param sd numeric vector with estimates of standard deviation(sigma) values for each component.
##' @details
##' components analysis uses the likelihood ratio test (LRT) to assess the number of mixture components.
##' Bad Random seed generation may cause normalmixEM to crash. Running the function again would be the quickest solution to this issue.
##'
##' @return
##' \describe{
##' \item{proportions}{The proportions of mixing components.}
##' \item{mean}{The Means of the components.}
##' \item{sd}{The Standard Deviations (SD) of components if not log transformed; if log transformed, it is then interpreted as Coefficient of Variation (CV).}
##' \item{loglik}{The log likelihood, useful for compare different fitting result, the bigger the better fit.}
##' }

##' @examples
##' # compare folders
##' folder1=system.file("extdata","SWR1",package="sojourner")
##' folder2=system.file("extdata","HTZ1",package="sojourner")
##' trackll=compareFolder(c(folder1,folder2))
##' MSD=msd(trackll=trackll)
##' dcoef=Dcoef(MSD,dt=6,plot=T,output=F)
##' # fit dcoef
##' a=fitNormDistr(dcoef,components=NULL,log.transform=F,combine.plot=F,output=F)
##' # to repeat a
##' b=fitNormDistr(dcoef,components=NULL,log.transform=F,combine.plot=F,output=F,seed=attr(a,"seed"))
##' # if a and b are the same
##' mapply(identical,a[[1]],b[[1]])
##' #try with log transformation
##' c=fitNormDistr(dcoef,components=2,log.transform=T,combine.plot=F,output=F)
##' #trying with some parameters provided(this will be applied to all dcoef results)
##' d=fitNormDistr(dcoef,meanRange=list(c(0.1,0.3), c(0.5,0.7)))


##' @export fitNormDistr
###############################################################################
# fit normal distribution to diffusion coefficient

#function for seed setting
.setSeed=function(seed=NULL){
    # set seed
    if (is.null(seed)){
        seed=sample(0:647,1)
        set.seed(seed)
    }else{set.seed(seed)}
    
    note=paste("\nRandom number generation seed",seed,"\n")
    cat(note)
    return(seed)
}

#function that deals with one-component normal distribution fitting
.singlecompFit=function(data){
    fit.info = fitdist(data=data, distr="norm")
    #This dummy matrix is used to make it work out with the gg.mixEM function in Plotting.Helpers.R
    dummy.posterior = matrix()
    colnames(dummy.posterior) = c("comp1")
    #a list with components similar to that of mixEM so that it can be used in gg.mixEM
    fit.output = list(x=as.vector(data), lambda = c(1), mu = fit.info$estimate[1], sigma = fit.info$sd[1],
                      loglik = fit.info$loglik, restarts=0, ft="normalmixEM", posterior= dummy.posterior)
    fit.se = list(mu.se= c(fit.info$estimate[2]), sigma.se = c(fit.info$sd[2]), lambda.se = c(0.00))
    return(list(fit.output,fit.se))
}

#uses mclust package's mclustBootstrapLRT() function to estimate the number of components in distribution
.getCompNum=function(data){
    cat("\ncomponents analysis\n")
    components.test=mclust::mclustBootstrapLRT(data,model = "V",verbose=T)
    print(components.test)
    components.int=length(components.test$p.value)
    cat("\n\nmost likely components",components.int,"at significant level 0.05\n\n")
    return(components.int)
}

.getSummaryResult=function(result, result.se, log.transform){
    out=list(
        data.frame(proportion=result$lambda,
                   proportion.se=result.se$lambda.se),
        
        #why returning it back from log?
        #data.frame(mean=if(log.transform) 10^(result$mu)
        #           else result$mu),
        data.frame(mean=result$mu),
        mean.se=as.vector(result.se$mu.se),
        
        
        data.frame(sd=result$sigma,
                   sd.se=result.se$sigma.se),
        
        log.lik=result$loglik
    )
    
    out=list(out)
    out=t(do.call(cbind.data.frame,out))
    colnames(out) = NULL
    return(out)
}

.fitNormDistr=function(dcoef,components=NULL,log.transform=F,binwidth=NULL,combine.plot=F,output=F,meanRange=NULL,proportion=NULL,means=NULL,sd=NULL, seed=NULL){
    # scale=1e3
    if(is.null(seed)){stop("NULL seed, cannot execute .fitNormDistr")}
    name=names(dcoef)
    len=length(dcoef)
    
    #parameters that help initial settings.
    params = list(meanRange, proportion, means, sd)

    mixmdl.lst=list()
    length(mixmdl.lst)=len
    names(mixmdl.lst)=name

    # store the standard error, currently only for components more than 2
    mixmdl.se.lst=list()
    length(mixmdl.se.lst)=len
    names(mixmdl.se.lst)=name
    
    ## This part deals with input related errors
    #collection of parameter vector lengths
    param.len = lapply(params, length)
    param.len=param.len[param.len!=0]
    if(!is.null(components) & (any(param.len != components))){
        stop("Invalid parameter size: number of components and parameter vector length does not match")
    }

    for (i in 1:len){
        data=dcoef[[name[i]]][,"slope"]

        # log transformation
        if (log.transform == T) {
            data=log10(data)
            data=data[!is.na(data)]
        }
        
        components.int=0
        #components test
        if (is.null(components)&length(param.len) == 0){
            components.int = .getCompNum(data)
        }else{
            if(length(param.len) == 0){
                components.int=components
            }else if(!all(as.vector(param.len) == param.len[[1]])){
                stop("Parameter vector lengths do not match")
            }
            if(is.null(components)){
                components.int=param.len[[1]]
            } else{components.int=components}
        }
        
        if(!is.null(meanRange)){
            if (length(param.len) == 1){
                #try to estimate/extract mu and sigma values from the given ranges
                means=unlist(lapply(meanRange, mean))
                diffs=lapply(meanRange, function(x){max(x)-mean(x)})
                sd=unlist(diffs)/1.6
            }
            else{stop("Either only meanRange or any of proportion/mean/sd should be provided")}
        }
        
        if (components.int == 1){ #The single component is a special case here since it cannot be a mixEM class object
            oneComp=.singlecompFit(data=data)
            fit.output=oneComp[[1]]
            fit.se=oneComp[[2]]
            plot.mixEM=gg.mixEM(fit.output,binwidth=binwidth,reorder=F)
            mixmdl=fit.output
            mixmdl.lst[[i]] = fit.output
            mixmdl.se=fit.se
            mixmdl.se.lst[[i]] = fit.se
        } else{
            mixmdl=normalmixEM(data,k=components.int,maxit=1e4,epsilon=1e-10, lambda=proportion, mu=means, sigma=sd)
            # reorder also for plotting and for output result
            # this is recordered on mean value
            mixmdl=reorderEM(mixmdl)
            print(summary(mixmdl))
            plot.mixEM=gg.mixEM(mixmdl,binwidth=binwidth,reorder=T)
            # approximate standard error using parametic bootstrap
            cat("\napproximating standard error by parametic bootstrap...\n\n")
            capture.output({mixmdl.se=boot.se(mixmdl, B = 1)})
            # file="/dev/null" # this only for mac
            mixmdl.lst[[i]]=mixmdl
            mixmdl.se.lst[[i]]=mixmdl.se
        }
        suppressMessages(plot(plot.mixEM))
    }

        result.lst=list()
        length(result.lst)=len
        names(result.lst)=name

        for (i in 1:len){
            result=.getSummaryResult(mixmdl.lst[[i]], mixmdl.se.lst[[i]], log.transform)
            result.lst[[i]]=result
        }

    print(result.lst)

    if (combine.plot == T){

        # same scale, same binwidth, same breaks
        ss=same.scale(mixmdl.lst)

        # auto binwidth, smaller of the two
        if (is.null(binwidth)) {
            binwidth.vec=sapply(mixmdl.lst,function(mdl){
                auto.binwidth(mdl$x)
            })

            binwidth=min(binwidth.vec)
            cat("\ncombined binwidth =",binwidth,"\n")
        }

        if(components.int == 1){reorder=F}
        else{reorder=T}
        plot.lst=lapply(mixmdl.lst,function(x){

            p=gg.mixEM(x,binwidth=binwidth,reorder=reorder)+
                # this only plot polygon but not histogram when ylim is added
                # +xlim(ss$scale.x)+ylim(ss$scale.y)
                coord_cartesian(xlim = ss$scale.x,ylim=ss$scale.y)+
                # makes integer breaks at the maxium n
                scale_x_continuous(breaks=scales::pretty_breaks(n=10))
            return(p)
        })

        # plot
        do.call(gridExtra::grid.arrange,plot.lst)
        # save
        cmb.plot=gridExtra::marrangeGrob(plot.lst, nrow=2, ncol=1)
    }

    # output
    if (output == T){
        result.df=do.call(rbind.data.frame,result.lst)
        logTrans=""
        if(log.transform){logTrans=".logtrans"}
        fileName=paste("FitNormDistr-",
                   .timeStamp(name[1]),".seed",seed,logTrans,"....csv",sep="")
        cat("\nOutput FitNormDistr.\n")
        write.csv(file=fileName,result.df)
        
        # output plot
        if (combine.plot == T){

            fileName=paste("FitNormDistr-combinePlot-",
                           .timeStamp(name[1]),".seed",seed,logTrans,"....pdf",sep="")
            cat("\nOutput FitNormDistr plot.\n")

            # gridExtra::marrangeGrob non-interactive use, multipage pdf
            # width=NULL,height=NULL default to graphic device size
            ggsave(filename=fileName,cmb.plot,width=8,height=8)
        }
    }
    # use invisible() so the user would not be overwhelmed by the numbers
    # at the same time programmers can assign the value and use it
    return(invisible(mixmdl.lst))
}


fitNormDistr=function(dcoef,components=NULL,log.transform=F,binwidth=NULL,combine.plot=F,output=F,seed=NULL, meanRange=NULL, proportion=NULL, means=NULL, sd=NULL){
    seed=.setSeed(seed=seed)

    # return
    structure(.fitNormDistr(dcoef=dcoef,components=components,log.transform=log.transform,binwidth=binwidth,combine.plot=combine.plot,output=output,
                            meanRange=meanRange,proportion=proportion,means=means,sd=sd,seed=seed),seed=seed)
}


# Comment:
# the calculation of means for the bootstrapped sample seems to be directly
# using fitted posterior distributions, rather than fitted 100 times.


# TODO: 
# isolate the boot.se part from fitNormDistr(), so the bootstrap is all in one
# place.
