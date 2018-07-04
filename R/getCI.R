## getCI-methods
##
##
########################################################################################################################
##' @name getCI
##' @aliases getCI
##' @title getCI
##' @rdname getCI-methods
##' @docType methods
##' @description given number of samples, return Confidence Interval information(mean, bounds, std err) for diffusion
##' coefficient and the proportions of each components.
##'
##' @usage
##' getCI(dcoef,plot=F, output=F)
##' @param dcoef diffusion coefficient calculated from Dcoef().
##' @param plot Logical indicating if plot should be generated.
##' @param output Logical indicating if output file should be generated.
##' @details
##' Supplied with number of samples as dcoef, uses that to calculate the confidence range.
##' Each entry in dcoeflist will consist of number of samples regarding that data.
##' May hit an error for bootstrapping due to random number generation, but re-running the function will do the job.
##' 
##'
##' @return
##' list of items, each of which will contain:
##' \describe{
##' \item{Estimate}{Mean estimated from sample}
##' \item{CI lower}{Lower bound of the confidence interval}
##' \item{CI upper}{Upper bound of the confidence interval}
##' \item{Std. Error}{Std. Error for given data}
##' }

##' @examples
##'
##' # compare folders
##' folder1=system.file("extdata","SWR1",package="sojourner")
##' folder2=system.file("extdata","HTZ1",package="sojourner")
##' trackll=compareFolder(c(folder1,folder2))
##' #get msd
##' MSD=msd(trackll=trackll)
##' #run Dcoef()
##' dcoef=Dcoef(MSD,dt=6,plot=T,output=F)
##' #resample data with bootstrapping(get 10 resampled data)
##' bootstraps=bootstrap(dcoef, attribute="slope")
##' # get confidence interval for this bootstrapped samples
##' a=getCI(bootstraps)
##' # to manually set confidence to 80% and generate output csv file
##' b=getCI(bootstraps, confidence=0.8, output=T)
##' 
##' @importFrom gmodels ci
## we deal with list of samples(possibly resampled)
###TODO add plot with confidence interval
require(gmodels)

##helper fucntion for calculating confidence interval for given attribute
.extractCI=function(fitResult, data, component, confidence){
    wanted=lapply(fitResult, get, x=data)
    return(ci(t(as.data.frame(wanted))[,component], confidence=confidence))
}

##' @export getCI
getCI=function(samples, confidence=0.95, output=F){
    inputNames = names(samples)
    ciList = list()
    for (i in 1:length(samples)){
        distances=list()
        #reformat the samples prior to running fitNormDistr
        boot.list=as.list(data.frame(t(samples[[i]])))
        boot.matricies = lapply(boot.list, as.matrix)
        listlen=length(boot.matricies)
        for (j in 1:listlen) {
            #fitNormDistr has to recognize this column name.
            colnames(boot.matricies[[j]]) <- "slope"
        }
        fit.result = fitNormDistr(boot.matricies, components = 2, combine.plot = T)
        
        #calculate confidence intervals for each components.
        firstcompCI=.extractCI(fit.result, "mu", 1, confidence)
        secondcompCI=.extractCI(fit.result, "mu", 2, confidence)
        firstProportionCI = .extractCI(fit.result, "lambda", 1, confidence)
        secondProportionCI = .extractCI(fit.result, "lambda", 2, confidence)
        
        cidf = data.frame(firstcompCI, secondcompCI, firstProportionCI, secondProportionCI)
        #append to list
        ciList = c(ciList,list(t(cidf)))
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