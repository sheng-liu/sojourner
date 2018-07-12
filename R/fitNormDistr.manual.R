## FitNormDistr.manual-methods
##
##
###############################################################################
##' @name fitNormDistr.manual
##' @aliases fitNormDistr.manual
##' @title fitNormDistr.manual
##' @rdname fitNormDistr.manual-methods
##' @docType methods
##' @description fit a normal distribution to diffusion coefficient caclulated with specified parameters
##'
##' @usage
##' fitNormDistr.manual(data, meanRange=NULL, proportion=NULL, means=NULL, 
##' sd=NULL, log.transform=F, output=F, components=NULL, binwidth=NULL, seed=NULL)
##' @param data single output from diffusion coefficient calculated from Dcoef(), like dcoef[[1]].
##' @param meanRange list of vectors, each of which have the min&max that estimates a range where a peak of the distribution should lie.
##' @param proportion numeric vector with estimates of each component's proportion of the whole data.
##' @param means numeric vector with estimates of mean(mu) values for each component.
##' @param sd numeric vector with estimates of standard deviation(sigma) values for each component.
##' @param log.transform logical indicate if log10 transformation is needed, default F.
##' @param output Logical indicaring if output file should be generated.
##' @param components parameter specifying the number of components to fit. If NULL (default), a components analysis is done to determine the most likely components and this number is then used for subsequent analysis.
##' @param binwidth binwidth for the combined plot. If NULL (default), will automatic assign binwidth.
##' @param seed Seed for random number generator. This makes each run easily repeatable. Seed will be automatically assigned if no seed is specified (default). The seed information is stored as an attribute of the returned object. The seed can also be output to a txt file when output=T.
##' @details
##' components analysis uses the likelihood ratio test (LRT) to assess the number of mixture components. This will be used if there are no estimates/information
##' regarding the data are provided. One should provide just the meanRange or the means, sd, and proportion. If both are provided, the function may raise an
##' error or run using the meanRange and ignoring other parameters.
##'
##' @return a single element from the list result that would result from FitNormDistr().
##' \describe{
##' \item{proportions}{The proportions of mixing components.}
##' \item{mean}{The Means of the components.}
##' \item{sd}{The Standard Deviations (SD) of components if not log transformed; if log transformed, it is then interpreted as Coefficient of Variation (CV).}
##' \item{loglik}{The log likelihood, useful for compare different fitting result, the bigger the better fit.}
##' }

##' @examples
##'
##' # compare folders
##' folder1=system.file("extdata","SWR1",package="sojourner")
##' folder2=system.file("extdata","HTZ1",package="sojourner")
##' trackll=compareFolder(c(folder1,folder2))
##' MSD=msd(trackll=trackll)
##' dcoef=Dcoef(MSD,dt=6,plot=T,output=F)
##' # fit first element of dcoef (SWR1)
##' a=fitNormDistr.manual(dcoef[[1]],components=NULL,log.transform=F,output=F)
##' # to repeat a
##' b=fitNormDistr.manual(dcoef[[2]],components=NULL,log.transform=F,output=F,seed=attr(a,"seed"))
##' # provide possible mean ranges assuming 2 component fit
##' c=fitNormDistr.manual(dcoef[[2]],meanRange=list(c(0.2,0.4),c(0.6,0.8)),log.transform=F,output=F)
##' 


##' @export fitNormDistr.manual
###############################################################################

#this will take in a matrix at a time instead of a list of matrices
fitNormDistr.manual=function(data, meanRange=NULL, proportion=NULL, means=NULL, sd=NULL, log.transform=F, output=F, components=NULL, binwidth=NULL, seed=NULL){
    note=.setSeed(seed=seed)
    # output seed
    if (output==T){
        name=names(dcoef)
        fileName=paste("FitNormDistr-",
                       .timeStamp(name[1]),"-Seed....txt",sep="")
        writeLines(text=note,con=fileName)
    }
    
    #parameters that help initial settings.
    params = list(meanRange, proportion, means, sd)
    slopes = data[,"slope"]
    # log transformation
    if (log.transform==T) {
        data=log10(data)
        data=data[!is.na(data)]
    }
    
    ## This part deals with input related errors
    #collection of parameter vector lengths
    param.len = lapply(params, length)
    param.len=param.len[param.len!=0]
    if(!is.null(components) & (any(param.len != components))){
        stop("Invalid parameter size: number of components and parameter vector length does not match")
    }
    
    components.int=0
    #components test
    if (is.null(components)&length(param.len)==0){
        components.int = .getCompNum(slopes)
    }else{
        if(length(param.len)==0){
            components.int=components
        }else if(!all(as.vector(param.len)==param.len[[1]])){
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
    
    if (components.int==1){ #The single component is a special case here since it cannot be a mixEM class object
        oneComp=.singlecompFit(data=slopes)
        fit.output=oneComp[[1]]
        fit.se=oneComp[[2]]
        plot.mixEM=gg.mixEM(fit.output,binwidth=binwidth,reorder=F)
        mixmdl=fit.output
        mixmdl.se=fit.se
    } else{
        mixmdl=normalmixEM(slopes,k=components.int,maxit=1e4,epsilon=1e-10, lambda=proportion, mu=means, sigma=sd)
        # reorder also for plotting and for output result
        # this is recordered on mean value
        mixmdl=reorderEM(mixmdl)
        print(summary(mixmdl))
        plot.mixEM=gg.mixEM(mixmdl,binwidth=binwidth,reorder=T)
        # approximate standard error using parametic bootstrap
        cat("\napproximating standard error by parametic bootstrap...\n\n")
        capture.output({mixmdl.se=boot.se(mixmdl, B = 100)})
        # file="/dev/null" # this only for mac
    }
    # supprress warnings
    suppressMessages(plot(plot.mixEM))
    
    #output to .csv file
    if (output==T){
        
        result.df=do.call(rbind.data.frame,list(.getSummaryResult(mixmdl, mixmdl.se, log.transform)))
        fileName=paste("FitNormDistr.manual-",
                       .timeStamp(name[1]),"....csv",sep="")
        cat("\nOutput FitNormDistr.manual.\n")
        write.csv(file=fileName,result.df)
    }
    return(invisible(mixmdl))
}