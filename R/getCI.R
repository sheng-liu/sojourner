## getCI-methods
##
##
########################################################################################################################
##' @name getCI
##' @aliases getCI
##' @title getCI
##' @rdname getCI-methods
##' @docType methods
##' @description given a dcoef result, return Confidence Interval information(mean, bounds, std err) for diffusion
##' coefficient and the proportions of each components.
##'
##' @usage
##' getCI(dcoef.result,confidence=0.95,output=F,n.reps=100,meanRange=NULL,proportion=NULL,means=NULL,sd=NULL,log.transform=F,components=NULL)
##' @param dcoef.result diffusion coefficient calculated from Dcoef().
##' @param confidence the level of confidence that is used to calculate the confidence interval.
##' @param output Logical indicate if output file should be generated.
##' @param n.reps number of replicates for bootstrapping
##' @param meanRange list of vectors, each of which have the min&max that estimates a range where a peak of the distribution should lie.
##' @param proportion numeric vector with estimates of each component's proportion of the whole data. Values in the vector would add up to 1.
##' @param means numeric vector with estimates of mean(mu) values for each component.
##' @param sd numeric vector with estimates of standard deviation(sigma) values for each component.
##' @param log.transform logical indicate if log10 transformation is needed, default F.
##' @param components parameter specifying the number of components to fit. If NULL (default), a components analysis is done to determine the most likely components and this number is then used for subsequent analysis.
##' @details
##' Supplied with a dcoef output, uses that to bootstrap and calculate the confidence range.
##' May hit an error for bootstrapping due to random number generation, but re-running the function will very likely resolve the problem.
##' The t-distribution/critical-t was used to calculate the Confidence Interval. 
##' For multi-component distributions, the boot.se function from the mixtools package was used.
##' For single-component distributions, a separate function was used to calculate the stderr and confidencet interval.
##' 
##'
##' @return
##' list of items, each of which will contain for each distribution component:
##' \describe{
##' \item{Estimate}{Mean estimated from sample}
##' \item{CI lower}{Lower bound of the confidence interval}
##' \item{CI upper}{Upper bound of the confidence interval}
##' \item{Std. Error}{Std. Error for given data}
##' }
##'
##' @examples
##' # compare folders
##' folder1=system.file("extdata","SWR1",package="sojourner")
##' folder2=system.file("extdata","HTZ1",package="sojourner")
##' trackll=compareFolder(c(folder1,folder2))
##' #get msd
##' MSD=msd(trackll=trackll)
##' #run Dcoef()
##' dcoef=Dcoef(MSD,dt=6,plot=T,output=F)
##' # get confidence intervals for this dcoef result which contains data from two different folders
##' a=getCI(dcoef)
##' # to manually set confidence to 80%
##' b=getCI(dcoef, confidence=0.8, output=F)
##' # get confidence interval with the some estimate values specified, same specification will be used for both data
##' e=getCI(dcoef, confidence=0.95, means=c(0.1, 0.4))
##' # If you want to apply to just one of this dcoef results
##' f=getCI(dcoef[1], confidence=0.95, means=c(0.1,0.4))
##' 
###TODO add plot with confidence interval

##small helper function that calculates the CI in terms of stderr. This returns the critical t-values, given the confidence.
.get.seRange=function(confidence, num.samples){
    diff=(1-confidence)/2
    range.vector=c(diff, 1-diff)
    range.se=qt(range.vector, num.samples-1)
    print(range.se)
    return(range.se)
}

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

.boot.se.all=function(fitResult, B){
    if(length(fitResult$mu) == 1){
        return(.boot.se.onecomp(fitResult,B))
    } else {
        capture.output({multi.stderr=boot.se(fitResult, B)})
        return(multi.stderr)
    }
}

##' @export getCI
getCI=function(dcoef.result, confidence=0.95, output=F, n.reps=100, meanRange=NULL, proportion=NULL, means=NULL, sd=NULL, log.transform=F, components=NULL){
    #from the desired confidecne lvl, calculate the crtitcal-t value
    SE.range=.get.seRange(confidence, n.reps)
    inputNames = names(dcoef.result)
    ciList = list()
    #fit and get multiple components
    fit.result=fitNormDistr(dcoef.result, meanRange=meanRange, proportion=proportion, means=means, sd=sd, log.transform=log.transform, components=components)
    #generate std.err value
    fit.se=lapply(fit.result, .boot.se.all, B=n.reps)
    for (i in 1:length(fit.result)){
        out.df = data.frame()
        
        #Extract values form the fitting result and then apply the calculated critical-t values stored in SE.range
        #generate the rownames in meannames and proportionnames each time we do this
        meannames=c()
        for (j in 1:length(fit.result[[i]]$mu)) {
            out.df = rbind(out.df, c(Estimate=fit.result[[i]]$mu[j], CI.lower=SE.range[1]*fit.se[[i]]$mu.se[j]+fit.result[[i]]$mu[j], 
                                     CI.upper=SE.range[2]*fit.se[[i]]$mu.se[j]+fit.result[[i]]$mu[j], Std.Error=fit.se[[i]]$mu.se[j]))
            meannames=c(meannames, paste(j,"-compCI",sep=""))
        }
        proportionnames=c()
        for (j in 1:length(fit.result[[i]]$lambda)) {
            out.df = rbind(out.df, c(Estimate=fit.result[[i]]$lambda[j], CI.lower=SE.range[1]*fit.se[[i]]$lambda.se[j]+fit.result[[i]]$lambda[j], 
                                     CI.upper=SE.range[2]*fit.se[[i]]$lambda.se[j]+fit.result[[i]]$lambda[j], Std.Error=fit.se[[i]]$lambda.se[j]))
            proportionnames=c(proportionnames, paste(j,"-ProportionCI", sep=""))
        }
        
        #add column/row names to the dataframe and then add the dataframe to list
        colnames(out.df) = c("Estimate", "CI.lower", "CI.upper", "Std.Error")
        rownames(out.df) = c(meannames,proportionnames)
        ciList[i] = list(out.df)
    }
    
    
    names(ciList) = inputNames
    #output .csv file has similar format as that of fitNormDistr
    if (output==T){
        whole.df = data.frame()
        fileName=paste("ConfidenceInterval-",
                       .timeStamp(paste(names(ciList),sep = "")),".csv",sep="")
        for (i in 1:length(ciList)){
            single.df = data.frame(ciList[[i]])
            rownames(single.df) = paste(inputNames[[i]],rownames(single.df), sep = ".")
            whole.df = rbind(whole.df, single.df)
        }
        write.csv(whole.df, file = fileName)
    }
    return(ciList)
}