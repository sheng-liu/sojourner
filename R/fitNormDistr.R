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
##'              seed=NULL, proportion=NULL, means=NULL, sd=NULL, constrain=F)
##' @param dcoef diffusion coefficient calculated from Dcoef().
##' @param components parameter specifying the number of components to fit. If NULL (default), a components analysis is done to determine the most likely components and this number is then used for subsequent analysis.
##' @param log.transform logical indicate if log10 transformation is needed, default F.
##' @param binwidth binwidth for the combined plot. If NULL (default), will automatic assign binwidth.
##' @param combine.plot Logical indicate if all the plot should be combined into one, with same scale (/same axises breaks), same color theme, and same bin size for comparison.
##' @param output logical indicate if output file should be generated.
##' @param seed Seed for random number generator. This makes each run easily repeatable. Seed will be automatically assigned if no seed is specified (default). The seed information is stored as an attribute of the returned object. The seed can also be output to a txt file when output=T.
##' @param proportion numeric vector with estimates of each component's proportion of the whole data.
##' @param means numeric vector with estimates of mean(mu) values for each component.
##' @param sd numeric vector with estimates of standard deviation(sigma) values for each component.
##' @param constrain logical indicate if mean and std deviation are set to the given value. This will not work for the unimodal distribution.
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
##' dcoef=Dcoef(MSD,dt=6,plot=TRUE,output=FALSE)
##' # fit dcoef
##' a=fitNormDistr(dcoef,components=NULL,log.transform=FALSE,combine.plot=FALSE,output=FALSE)
##' # to repeat a
##' b=fitNormDistr(dcoef,components=NULL,log.transform=FALSE,
##' combine.plot=FALSE,output=FALSE,seed=attr(a,"seed"))
##' # if a and b are the same
##' mapply(identical,a[[1]],b[[1]])
##' #try with log transformation
##' c=fitNormDistr(dcoef,components=2,log.transform=TRUE,combine.plot=FALSE,output=FALSE)
##' #trying with some parameters provided(this will be applied to all dcoef results). 
##' #with constrain = F, this will be used as the starting values for the EM-algorithm
##' #normally we should deal with only one dataset when working with constrains, 
##' #since it will apply to all of them.
##' folder3=system.file("extdata","HSF", package="sojourner")
##' trackll=compareFolder(c(folder3),input=2)
##' MSD=msd(trackll=trackll)
##' dcoef=Dcoef(MSD,dt=6,plot=TRUE,output=FALSE)
##' #try with constrain =T, the values will be forced to eqaul the provided ones.
##' e=fitNormDistr(dcoef,means=c(0.3,0.5), constrain=TRUE)


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
.singlecompFit=function(data, mean=NULL, sd=NULL, constrain=F){
    if (constrain == T){warning("single component fitting cannot be constrained.")}
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
    components.test=mclust::mclustBootstrapLRT(data,model = "V",verbose=T)
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

.fitNormDistr=function(dcoef,components=NULL,log.transform=F,binwidth=NULL,combine.plot=F,output=F,proportion=NULL,means=NULL,sd=NULL,seed=NULL,constrain=F){
    # scale=1e3
    if(is.null(seed)){stop("NULL seed, cannot execute .fitNormDistr")}
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
        
        
        if (components.int == 1){ #The single component is a special case here since it cannot be a mixEM class object
            if(constrain==T){
                oneComp = .singlecompFit(data=data, mean=means, sd=sd, constrain=T)
            } else {oneComp=.singlecompFit(data=data)}
            fit.output=oneComp[[1]]
            plot.mixEM=gg.mixEM(fit.output,binwidth=binwidth,reorder=F)
            mixmdl=fit.output
            mixmdl.lst[[i]] = fit.output
        } else{
            if(constrain == T) {
                mixmdl=normalmixEM(data,k=components.int,maxit=1e4,epsilon=1e-10,lambda=proportion,sd.constr=sd,mean.constr=means)
            } else {
                mixmdl=normalmixEM(data,k=components.int,maxit=1e4,epsilon=1e-10,lambda=proportion,mu=means,sigma=sd)
            }
            # reorder also for plotting and for output result
            # this is recordered on mean value
            mixmdl=reorderEM(mixmdl)
            print(summary(mixmdl))
            plot.mixEM=gg.mixEM(mixmdl,binwidth=binwidth,reorder=T)
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


fitNormDistr=function(dcoef,components=NULL,log.transform=F,binwidth=NULL,combine.plot=F,output=F,seed=NULL, proportion=NULL, means=NULL, sd=NULL, constrain=F){
    seed=.setSeed(seed=seed)

    # return
    structure(.fitNormDistr(dcoef=dcoef,components=components,log.transform=log.transform,binwidth=binwidth,combine.plot=combine.plot,output=output,
                            proportion=proportion,means=means,sd=sd,constrain=constrain,seed=seed),seed=seed)
}


# Comment:
# the calculation of means for the bootstrapped sample seems to be directly
# using fitted posterior distributions, rather than fitted 100 times.


# TODO: 
# constrain for the single component case should be implemented.
