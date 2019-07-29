## getCI-methods
##
##
###############################################################################
##' @name getCI
##' @aliases getCI
##' @title getCI
##' @rdname getCI-methods
##' @docType methods
##' @description given a dcoef result, return Confidence Interval
##' information(mean, bounds, std err) for diffusion coefficient and the
##' proportions of each components.
##'
##' @usage
##' getCI(bootstrap.result,confidence=0.95,output=FALSE)
##' @param bootstrap.result diffusion coefficient calculated from Dcoef().
##' @param confidence the level of confidence that is used to calculate the
##' confidence interval.
##' @param output Logical indicate if output file should be generated.
##' @details Supplied with a bootstrap output, it calculates the confidence
##' range. The t-distribution/critical-t was used to calculate the Confidence
##' Interval.
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
##' trackll=compareFolder(folders=c(folder1,folder2), input=3)
##' #get msd
##' MSD=msd(trackll=trackll)
##' #run Dcoef()
##' dcoef=Dcoef(MSD,dt=6,plot=TRUE,output=FALSE)
##' #fit the dcoef result
##' normalFit=fitNormDistr(dcoef)
##' # perform bootstrapping for this dcoef result
##' d.boot = bootstrap(normalFit, n.reps=100)
##' # get confidence intervals for this dcoef result which contains data from 
##' # two different folders
##' a=getCI(d.boot)
##' # to manually set confidence to 80%
##' b=getCI(d.boot, confidence=0.8, output=FALSE)
##' 
###TODO add plot with confidence interval

##small helper function that calculates the CI in terms of stderr. This 
##returns the critical t-values, given the confidence.
.get.seRange=function(confidence, num.samples){
    diff=(1-confidence)/2
    range.vector=c(diff, 1-diff)
    range.se=qt(range.vector, num.samples-1)
    return(range.se)
}

##' @export getCI
getCI=function(bootstrap.result, confidence=0.95, output=FALSE){
    inputNames = names(bootstrap.result$Fit)
    ciList = list()
    #get fit component
    fit.result=bootstrap.result$Fit
    #get std.err component
    fit.se=bootstrap.result$Bootstraps
    for (i in seq_along(fit.result)){
        #from the desired confidecne lvl, calculate the crtitcal-t value
        n.reps=fit.se[[i]]$n.reps
        SE.range=.get.seRange(confidence, n.reps)
        out.df = data.frame()
        
        #Extract values form the fitting result and then apply the calculated
        #critical-t values stored in SE.range generate the rownames in 
        #meannames and proportionnames each time we do this
        meannames = c()
        for (j in seq_along(fit.result[[i]]$mu)) {
            out.df = rbind(
                out.df,
                c(
                    Estimate = fit.result[[i]]$mu[j],
                    CI.lower = SE.range[1] * fit.se[[i]]$mu.se[j] + 
                        fit.result[[i]]$mu[j],
                    CI.upper = SE.range[2] * fit.se[[i]]$mu.se[j] +
                        fit.result[[i]]$mu[j],
                    Std.Error = fit.se[[i]]$mu.se[j]
                )
            )
            meannames = c(meannames, paste(j, "-compCI", sep = ""))
        }
        proportionnames = c()
        for (j in seq_along(fit.result[[i]]$lambda)) {
            out.df = rbind(
                out.df,
                c(
                    Estimate = fit.result[[i]]$lambda[j],
                    CI.lower = SE.range[1] * fit.se[[i]]$lambda.se[j] + 
                                fit.result[[i]]$lambda[j],
                    CI.upper = SE.range[2] * fit.se[[i]]$lambda.se[j] +
                        fit.result[[i]]$lambda[j],
                    Std.Error = fit.se[[i]]$lambda.se[j]
                )
            )
            proportionnames = c(proportionnames, 
                                paste(j, "-ProportionCI", sep =
                                                            ""))
        }
        
        #add column/row names to the dataframe and add the dataframe to list
        colnames(out.df) = c("Estimate", "CI.lower", "CI.upper", "Std.Error")
        rownames(out.df) = c(meannames, proportionnames)
        ciList[i] = list(out.df)
    }
    
    names(ciList) = inputNames
    #output .csv file has similar format as that of fitNormDistr
    if (output == TRUE){
        whole.df = data.frame()
        fileName=paste("ConfidenceInterval-",
                        .timeStamp(paste(names(ciList),sep = "")),".csv",
                       sep="")
        for (i in seq_along(ciList)){
            single.df = data.frame(ciList[[i]])
            rownames(single.df) = paste(inputNames[[i]],rownames(single.df), 
                                        sep = ".")
            whole.df = rbind(whole.df, single.df)
        }
        write.csv(whole.df, file = fileName)
    }
    print(ciList)
    return(ciList)
}
