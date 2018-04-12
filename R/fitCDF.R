## fitCDF-methods
##
##
###############################################################################
##' @name fitCDF
##' @aliases fitCDF
##' @title fitCDF
##' @rdname fitCDF-methods
##' @docType methods
##' @description Caclulate apparent diffusion coefficient (Dcoef) for
##'   trajecotries by fitting displacementCDF.
##'
##' @usage
##'
##' fitCDF((cdf, components=c("one","two","three"),
##'         start=list(
##'             oneCompFit=list(D=c(0,2)),
##'             twoCompFit=list(D1=c(0,2),D2=c(0,2),alpha=c(0,1)),
##'             threeCompFit=list(D1=c(0,2),D2=c(0,2),D3=c(0,2),
##'                               alpha=c(0,1),beta=c(0,1))),
##'         t.interval=0.01,
##'         maxiter.search=1000,
##'         maxiter.optim=1000,
##'         output=F,
##'         seed=NULL)
##'
##' @param cdf cdf calculated from displacementCDF().
##' @param components parameter specifying the number of components to fit.Currently support one to three components fit.
##' @param start the start value for fitting.
##' @param t.interval time interval for image aquisition. Default 0.01 sec.
##' @param maxiter.search maximum iteration in random search start value process. defual to 1000.
##' @param maxiter.optim maximum iteration in local optimization process. Default ot 1000.
##' @param output Logical indicaring if output file should be generated.
##' @param seed Seed for random number generator. This makes each run easily repeatable. Seed will be automatically assigned if no seed is specified (default). The seed information is stored as an attribute of the returned object. The seed can also be output to a txt file when output=T.
##' @return
##' \itemize{
##' \item{on screen output and file} Result and parameters of goodness of the fit.
##' \item{Plot,} fiting plot.
##' }
##' @details calculating Dceof by fitting displacementCDF.
##'
##' Reducing the range can greatly increase the precision of the searching; alternatively, if the range are unavailable, increase the maxiter.search so more points will be searched through with the cost of computation time. maxiter.optim barely need to change, if it does not converge with default setting maxiter=1000, most likely the problem is in the initial values.

##'
##' @examples
##'
##' # compare folders
##' folder1=system.file("extdata","SWR1",package="sojourner")
##' folder2=system.file("extdata","HTZ1",package="sojourner")
##' trackll=compareFolder(c(folder1,folder2))
##' cdf=displacementCDF(trackll,dt=1,plot=F,output=F)
##'
##' # specify ranges of parameter value of interest
##' fitCDF(cdf,components="two",
##'       start=list(
##'                 twoCompFit=list(D1=c(0,2),D2=c(0,2),alpha=c(0,1)))
##'                 )
##'
##' # repeat a fit
##' a=fitCDF(cdf,components="two",output=F)
##' b=fitCDF(cdf,components="two",output=F,seed=attr(a,"seed"))
##'
##' # if result are identical
##' x=summary(a[[1]])
##' y=summary(b[[1]])
##' # formula records environment, exclude from the comparison
##' mapply(identical,x[names(x)!="formula"],y[names(y)!="formula"])


# the only difference is the function environment, of course it is runned in two functions
# > mapply(identical,summary(a[[1]]),summary(b[[1]]))
# formula    residuals        sigma           df cov.unscaled         call
# FALSE         TRUE         TRUE         TRUE         TRUE         TRUE
# convInfo      control    na.action coefficients   parameters
# TRUE         TRUE         TRUE         TRUE         TRUE
# > summary(b[[1]])$formula
# P ~ p3(r, D1, D2, alpha)
# <environment: 0x11a893640>
#     > summary(a[[1]])$formula
# P ~ p3(r, D1, D2, alpha)
# <environment: 0x11e79f208>

##
###############################################################################

## dt needs to be avariable in this equation so it is flexible and has a meaning to its unit um/s

# ------------------------------------------------------------------------------
# one component fit
# library(truncnorm)

one.comp.fit=function(r,P,start=list(D=c(0,3)),t.interval=0.01,maxiter.optim=1000,name){
    # with one parameter, D
    p1 = function(r,D){1 - exp(-r^2/(4*D*t.interval))}

    title=paste("One component fit -",name)
    cat("\n\n","==>",title,"\n")

    # it seems for one parameter fit, any nls program does really well, in
    # another words, the fitting seems independent of the intial value.
    # maybe one parameter makes the relationship linear, so it always reaches to
    # one conclusion

    D=truncnorm::rtruncnorm(1,a=data.frame(start)[1,],b=data.frame(start)[2,])
    # this defaults normal distribution with mean = 0, sd = 1

    start=list(D=D)

    # fit equation 2 to data P
    ocfit=nls(P ~ p1(r,D),start=start,control = nls.control(maxiter = maxiter.optim))

    print(ocfit);cat("\n")

    # plot
    pdf(paste(Sys.time(), ".pdf", sep = ""))
    plot(r,P,main=title,cex=0.3)
    curve(p1(x,D=coef(ocfit)),add=TRUE,col="red")
    dev.off()
    
    plot(r,P,main=title,cex=0.3)
    curve(p1(x,D=coef(ocfit)),add=TRUE,col="red")
    
    #coef(summary(ocfit))

    return(ocfit)
}

# ------------------------------------------------------------------------------
# two components fit

two.comp.fit=function(r,P,start=list(D1=c(0,2),D2=c(0,2),alpha=c(0,1)),
                      t.interval=0.01,maxiter.search=1000,maxiter.optim=1000,name){

    ## equation
    p3 =function(r,D1,D2,alpha){
        1 - (alpha*exp(-r^2/(4*D1*t.interval)) + 
                 (1-alpha)*exp(-r^2/(4*D2*t.interval)))}

    title=paste("Two components fit -",name)
    cat("\n\n","==>",title,"\n")

    cat("\nBrute force random search start value...\n\n")
    r.search.tcfit=nls2(P ~ p3(r,D1,D2,alpha),start=data.frame(start),
                       # algorithm="brute-force",
                       algorithm="random-search",
                       control = nls.control(maxiter = maxiter.search))

    print(coef(r.search.tcfit))

    ## local optimization using minpack.lm::nlsLM
    cat("\nLocal optimization...\n\n")

    tcfit=minpack.lm::nlsLM(P ~ p3(r,D1,D2,alpha),
                start=coef(r.search.tcfit),
                # three paramert ~D1, D2, alpha
                lower=c(start$D1[[1]],start$D2[[1]],start$alpha[[1]]),
                upper=c(start$D1[[2]],start$D2[[2]],start$alpha[[2]]),
                
                # barely need control of maxiter for minpack.lm::nlsLM, if it does not
                # converge with default setting maxiter=1024, most likely the
                # problem is in the initial values
                control = nls.control(maxiter = maxiter.optim)
                )

    print(tcfit);cat("\n")
    
    pdf(paste(Sys.time(), ".pdf", sep = ""))
    ## plotting
    plot(r,P,main=title,cex=0.3)
    curve(p3(x,
             coef(tcfit)["D1"],
             coef(tcfit)["D2"],
             coef(tcfit)["alpha"]),
          add=T,col="red"
    )
    dev.off()
    plot(r,P,main=title,cex=0.3)
    curve(p3(x,
             coef(tcfit)["D1"],
             coef(tcfit)["D2"],
             coef(tcfit)["alpha"]),
          add=T,col="red"
    )

    return(tcfit)
}
# ------------------------------------------------------------------------------
# three components fit

three.comp.fit=function(r,P,start=list(D1=c(0,2),D2=c(0,2),D3=c(0,2),
                                       alpha=c(0,1),beta=c(0,1)),
                        t.interval=0.01,maxiter.search=1000,maxiter.optim=1000,name){


    ## equation
    p5=function(r,D1,D2,D3,alpha,beta){
        1 - (
            alpha*exp(-r^2/(4*D1*t.interval)) +
                beta*exp(-r^2/(4*D2*t.interval)) +
                (1-alpha-beta)*exp(-r^2/(4*D3*t.interval))
        )}

    title=paste("Three components fit -",name)
    cat("\n\n","==>",title,"\n")

    cat("\nBrute force random search start value...\n\n")
    r.search.thcfit=nls2(P ~ p5(r,D1,D2,D3,alpha,beta),start=data.frame(start),

                        # the virtue of using random search is it reduces the
                        # chance that one parameters is fixed at the edge value,
                        # which is barely the case.

                        # algorithm="brute-force",
                        algorithm="random-search",
                        control = nls.control(maxiter = maxiter.search))

    print(coef(r.search.thcfit))

    ## local optimization using minpack.lm::nlsLM
    ## from minpack.lm for low or zero noise data
    cat("\nLocal optimization...\n\n")
    thcfit=minpack.lm::nlsLM(P ~ p5(r,D1,D2,D3,alpha,beta),
                 start=coef(r.search.thcfit),
                 # lower=c(0,0,0,0,0),
                 # upper=c(Inf,Inf,Inf,1,1),
                 
                 lower=c(start$D1[[1]],start$D2[[1]],start$D3[[1]],start$alpha[1],start$beta[[1]]),
                 upper=c(start$D1[[2]],start$D2[[2]],start$D3[[2]],start$alpha[[2]],start$beta[[2]]),

                 # barely need control of maxiter for minpack.lm::nlsLM, if it does not
                 # converge with default setting maxiter=1024, most likely the
                 # problem is in the initial values
                 control = nls.control(maxiter = maxiter.optim)
                 )

    print(thcfit);cat("\n") # needed for print when compiled as pacakge
    
    pdf(paste(Sys.time(), ".pdf", sep = ""))
    ## plot
    plot(r,P,main=title,cex=0.3)
    curve(p5(x,coef(thcfit)["D1"],coef(thcfit)["D2"],coef(thcfit)["D3"],
             coef(thcfit)["alpha"],coef(thcfit)["beta"]),add=T,col="red")
    dev.off()
    
    plot(r,P,main=title,cex=0.3)
    curve(p5(x,coef(thcfit)["D1"],coef(thcfit)["D2"],coef(thcfit)["D3"],
             coef(thcfit)["alpha"],coef(thcfit)["beta"]),add=T,col="red")

    return(thcfit)

}




# ------------------------------------------------------------------------------
# fitCDF
##' @export fitCDF
.fitCDF=function(cdf, components=c("one","two","three"),
                start=list(
                    oneCompFit=list(D=c(0,2)),
                    twoCompFit=list(D1=c(0,2),D2=c(0,2),alpha=c(0,1)),
                    threeCompFit=list(D1=c(0,2),D2=c(0,2),D3=c(0,2),
                                      alpha=c(0,1),beta=c(0,1))),
                t.interval=0.01,
                maxiter.search=1000,
                maxiter.optim=1000,
                output=F){

    # use lapply to do it for all folders
    cdf.displacement=cdf$CDF.displacement
    name=names(cdf.displacement)

    method=match.arg(components)

    fit=list()
    length(fit)=length(cdf.displacement)
    names(fit)=names(cdf.displacement)

    for (i in 1:length(cdf.displacement)){

        r=cdf.displacement[[i]]$UniqueDisplacement
        P=cdf.displacement[[i]]$CDF


        fit[[i]]=switch(method,
               one={one.comp.fit(r,P,start=start$oneCompFit,
                                 t.interval=t.interval,name[i],
                                 maxiter.optim=maxiter.optim)},

               two={two.comp.fit(r,P,start=start$twoCompFit,
                                 t.interval=t.interval,name[i],
                                 maxiter.search=maxiter.search,
                                 maxiter.optim=maxiter.optim)},

               three={three.comp.fit(r,P,start=start$threeCompFit,
                                     t.interval=t.interval,name[i],
                                     maxiter.search=maxiter.search,
                                     maxiter.optim=maxiter.optim)}
               )

    }

    result.lst=lapply(fit,function(x) coef(summary(x)))



#     lapply(fit,function(x) coef(x))
#     lapply(fit,function(x) summary(x))
#     lapply(fit,function(x) logLik(x))
#
#     lapply(fit,function(x) sum(residuals(x)^2)) # residual sum-of-squares
#     lapply(fit,function(x) deviance(x)) # is/equals to residual sum-of-squares

    result.lst=lapply(fit,function(x) {
        estimate=coef(summary(x))[,"Estimate"]
        std.error=coef(summary(x))[,"Std. Error"]
        res.sum.of.squares=deviance(x)
        re=cbind(estimate,std.error,res.sum.of.squares)}
        )


    print(result.lst)

    # output
    if (output==T){

        result.df=do.call(rbind.data.frame,result.lst)
        fileName=paste("FitCDF-",
                       .timeStamp(name[1]),"....csv",sep="")
        cat("\nOutput FitCDF.\n")
        write.csv(file=fileName,result.df)
    }

    return(invisible(fit))

}

fitCDF=function(cdf, components=c("one","two","three"),
                start=list(
                    oneCompFit=list(D=c(0,2)),
                    twoCompFit=list(D1=c(0,2),D2=c(0,2),alpha=c(0,1)),
                    threeCompFit=list(D1=c(0,2),D2=c(0,2),D3=c(0,2),
                                      alpha=c(0,1),beta=c(0,1))),
                t.interval=0.01,
                maxiter.search=1000,
                maxiter.optim=1000,
                output=F,
                seed=NULL){

    # set seed
    result=seedIt(expr=.fitCDF(cdf=cdf, components=components,
                               start=start,
                               t.interval=t.interval,
                               maxiter.search=maxiter.search,
                               maxiter.optim=maxiter.optim,
                               output=output),seed=seed)

    # output seed
    if (output==T){

        name=names(result)
        fileName=paste("FitCDF-",
                       .timeStamp(name[1]),"-Seed....txt",sep="")

        note <- paste("\nRandom number generation seed",attr(result,"seed"),"\n")
        writeLines(text=note,con=fileName)
    }

    return(result)

}







# ------------------------------------------------------------------------------
# DONE:

# a better density plot than ggplot2 there, or adjust it to be better /professiona looking
#   plot(density(r),main="distribution of displacement r")

## output summary into csv files


# the range of diffusion coefficient
# In dilute aqueous solutions the diffusion coefficients of most ions are
# similar and have values that at room temperature are in the range of 0.6 ×
# 10−9 to 2 × 10−9 m2/s. For biological molecules the diffusion coefficients
# normally range from 10−11 to 10−10 m2/s.

# https://en.wikipedia.org/wiki/Fick%27s_laws_of_diffusion
