## fitNormDistr-methods
##
##
###############################################################################
##' @name fitNormDistr
##' @aliases fitNormDistr
##' @title fitNormDistr
##' @rdname fitNormDistr-methods
##' @docType methods
##' @description fit normal distributions to diffusion coefficient caclulated by Dcoef method and saves seed state as a attribute of the result
##'
##' @usage
##' fitNormDistr(dcoef,components=NULL,log.transform=FALSE,binwidth=NULL,
##' combine.plot=FALSE,output=FALSE, proportion=NULL, means=NULL,
##' sd=NULL, constrain=FALSE)
##' @param dcoef diffusion coefficient calculated from Dcoef().
##' @param components parameter specifying the number of components to fit. If NULL (default), a components analysis is done to determine the most likely components and this number is then used for subsequent analysis.
##' @param log.transform logical indicate if log10 transformation is needed, default F.
##' @param binwidth binwidth for the combined plot. If NULL (default), will automatic assign binwidth.
##' @param combine.plot Logical indicate if all the plot should be combined into one, with same scale (/same axises breaks), same color theme, and same bin size for comparison.
##' @param output logical indicate if output file should be generated.
##' @param proportion numeric vector with estimates of each component's proportion of the whole data.
##' @param means numeric vector with estimates of mean(mu) values for each component.
##' @param sd numeric vector with estimates of standard deviation(sigma) values for each component.
##' @param constrain logical indicate if mean and std deviation are set to the given value. This will not work for the unimodal distribution.
##' @details
##' Components analysis uses the likelihood ratio test (LRT) to assess the number of mixture components.
##' Bad Random seed generation may cause normalmixEM to crash. Using another seed will solve the issue.
##' 
##' Note: Ensure that a random number generator seed has been manually set! The seed is stored as an attribute of the returned object of fitNormDistr() and using the same seed makes results repeatable (see examples).
##' 
##' @return
##' \describe{
##' \item{proportions}{The proportions of mixing components.}
##' \item{mean}{The Means of the components.}
##' \item{sd}{The Standard Deviations (SD) of components if not log transformed; if log transformed, it is then interpreted as Coefficient of Variation (CV).}
##' \item{loglik}{The log likelihood, useful for compare different fitting result, the bigger the better fit.}
##' }
##' 
##' @examples
##' # compare folders
##' folder1=system.file("extdata","SWR1",package="sojourner")
##' folder2=system.file("extdata","HTZ1",package="sojourner")
##' trackll=compareFolder(folders=c(folder1,folder2), input=3)
##' MSD=msd(trackll=trackll)
##' dcoef=Dcoef(MSD,dt=6,plot=TRUE,output=FALSE)
##' 
##' # set unique seed (use any number)
##' set.seed(123)
##' 
##' # fit dcoef (function automatically saves seed state as an attribute of the result)
##' a=fitNormDistr(dcoef,components=NULL,log.transform=FALSE,combine.plot=FALSE,output=FALSE)
##' 
##' # to repeat results of 'a', load seed attribute of a into RNG state
##' .Random.seed=attr(a,"seed")
##' # or, reset the seed with same unique number
##' # set.seed(123)
##' 
##' b=fitNormDistr(dcoef,components=NULL,log.transform=FALSE,
##' combine.plot=FALSE,output=FALSE)
##' 
##' # if 'a' and 'b' are the same
##' mapply(identical,a[[1]],b[[1]])
##' 
##' #try with log transformation
##' set.seed(234)
##' c=fitNormDistr(dcoef,components=2,log.transform=TRUE,combine.plot=FALSE,output=FALSE)
##' 
##' # trying with some parameters provided(this will be applied to all dcoef results). 
##' # with constrain = FALSE, this will be used as the starting values for the EM-algorithm
##' # normally we should deal with only one dataset when working with constrains, 
##' # since it will apply to all of them.
##' folder3=system.file("extdata","HSF", package="sojourner")
##' trackll=compareFolder(c(folder3),input=3)
##' MSD=msd(trackll=trackll)
##' dcoef=Dcoef(MSD,dt=6,plot=TRUE,output=FALSE)
##' 
##' # try with constrain =TRUE, the values will be forced to eqaul the provided ones.
##' set.seed(345)
##' e=fitNormDistr(dcoef,means=c(0.3,0.5), constrain=TRUE)

##' @export fitNormDistr
###############################################################################
# fit normal distribution to diffusion coefficient

#function that deals with one-component normal distribution fitting
.singlecompFit=function(data, mean=NULL, sd=NULL, constrain=FALSE){
    if (constrain == TRUE){warning("single component fitting cannot be constrained.")}
    #    if(!is.null(mean) & !is.null(sd)){stop("Please provide only one of the mean or sd constraints.")}
    #    fixlist = list(c(mean=mean, sd=sd))
    #    fit.info = fitdist(data,"norm", fix.arg=list(mean=mean))
    #} else{fit.info = fitdist(data=data, distr="norm")}
    #This dummy matrix is used to make it work out with the gg.mixEM function in Plotting.Helpers.R
    fit.info = fitdist(data=data, distr="norm")
    dummy.posterior = matrix()
    colnames(dummy.posterior) = c("comp1")
    #a list with components similar to that of mixEM so that it can be used in gg.mixEM
    fit.output = list(x=as.vector(data), lambda = c(1), mu = fit.info$estimate[1], sigma = fit.info$sd[1],
                        loglik = fit.info$loglik, restarts=0, ft="normalmixEM", posterior= dummy.posterior)
    return(list(fit.output))
}

#uses mclust package's mclustBootstrapLRT() function to estimate the number of components in distribution
.getCompNum=function(data){
    cat("\ncomponents analysis\n")
    components.test=mclust::mclustBootstrapLRT(data,model = "V",verbose=TRUE)
    print(components.test)
    components.int=length(components.test$p.value)
    cat("\n\nmost likely components",components.int,"at significant level 0.05\n\n")
    return(components.int)
}

.getSummaryResult=function(result, log.transform){
    out=list(
        data.frame(proportion=result$lambda),

        #why returning it back from log?
        #data.frame(mean=if(log.transform) 10^(result$mu)
        #           else result$mu),
        data.frame(mean=result$mu),
        
        data.frame(sd=result$sigma),

        log.lik=result$loglik
    )
    
    out=list(out)
    out=t(do.call(cbind.data.frame,out))
    colnames(out) = NULL
    return(out)
}

.fitNormDistr=function(dcoef,components=NULL,log.transform=FALSE,binwidth=NULL,combine.plot=FALSE,output=FALSE,proportion=NULL,means=NULL,sd=NULL,constrain=FALSE){
    cat("\nIMPORTANT: Ensure a seed has been manually set! See help docs for more info.\n")
    name=names(dcoef)
    len=length(dcoef)
    
    #parameters that help initial settings.
    params = list(proportion, means, sd)

    mixmdl.lst=list()
    length(mixmdl.lst)=len
    names(mixmdl.lst)=name

    # store the standard error, currently only for components more than 2
    
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
        if (log.transform == TRUE) {
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
        
        
        if (components.int == 1){ #The single component is a special case here since it cannot be a mixEM class object
            if(constrain==TRUE){
                oneComp = .singlecompFit(data=data, mean=means, sd=sd, constrain=TRUE)
            } else {oneComp=.singlecompFit(data=data)}
            fit.output=oneComp[[1]]
            plot.mixEM=gg.mixEM(fit.output,binwidth=binwidth,reorder=FALSE)
            mixmdl=fit.output
            mixmdl.lst[[i]] = fit.output
        } else{
            if(constrain == TRUE) {
                mixmdl=normalmixEM(data,k=components.int,maxit=1e4,epsilon=1e-10,lambda=proportion,sd.constr=sd,mean.constr=means)
            } else {
                mixmdl=normalmixEM(data,k=components.int,maxit=1e4,epsilon=1e-10,lambda=proportion,mu=means,sigma=sd)
            }
            # reorder also for plotting and for output result
            # this is recordered on mean value
            mixmdl=reorderEM(mixmdl)
            print(summary(mixmdl))
            plot.mixEM=gg.mixEM(mixmdl,binwidth=binwidth,reorder=TRUE)
            # approximate standard error using parametic bootstrap
            cat("\napproximating standard error by parametic bootstrap...\n\n")
            # file="/dev/null" # this only for mac
            mixmdl.lst[[i]]=mixmdl
        }
        suppressMessages(plot(plot.mixEM))
    }

        result.lst=list()
        length(result.lst)=len
        names(result.lst)=name

        for (i in 1:len){
            result=.getSummaryResult(mixmdl.lst[[i]], log.transform)
            result.lst[[i]]=result
        }

    print(result.lst)

    if (combine.plot == TRUE){

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

        if(components.int == 1){reorder=FALSE}
        else{reorder=TRUE}
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
    if (output == TRUE){
        result.df=do.call(rbind.data.frame,result.lst)
        logTrans=""
        if(log.transform){logTrans=".logtrans"}
        fileName=paste("FitNormDistr-",
                    .timeStamp(name[1]),logTrans,"....csv",sep="")
        cat("\nOutput FitNormDistr.\n")
        write.csv(file=fileName,result.df)
        
        # output plot
        if (combine.plot == TRUE){

            fileName=paste("FitNormDistr-combinePlot-",
                            .timeStamp(name[1]),logTrans,"....pdf",sep="")
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


fitNormDistr=function(dcoef,components=NULL,log.transform=FALSE,binwidth=NULL,combine.plot=FALSE,output=FALSE, proportion=NULL, means=NULL, sd=NULL, constrain=FALSE){
    
    # collects current seed (recommended that a unique seed is set beforehand, e.g. set.seed(123))
    my_seed=.Random.seed
    
    result = .fitNormDistr(dcoef=dcoef,components=components,log.transform=log.transform,binwidth=binwidth,combine.plot=combine.plot,output=output,
                            proportion=proportion,means=means,sd=sd,constrain=constrain)
    
    # saves seed as attribute of result
    attr(result,"seed")=my_seed
    
    return(result)
}


# Comment:
# the calculation of means for the bootstrapped sample seems to be directly
# using fitted posterior distributions, rather than fitted 100 times.


# TODO: 
# constrain for the single component case should be implemented.
