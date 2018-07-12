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
##' fitNormDistr(dcoef,components=NULL,log.transform=F,binwidth=NULL,
##'              combine.plot=F,output=F,seed=NULL)
##' @param dcoef diffusion coefficient calculated from Dcoef().
##' @param components parameter specifying the number of components to fit. If NULL (default), a components analysis is done to determine the most likely components and this number is then used for subsequent analysis.
##' @param log.transform logical indicate if log10 transformation is needed, default F.
##' @param binwidth binwidth for the combined plot. If NULL (default), will automatic assign binwidth.
##' @param output Logical indicaring if output file should be generated.
##' @param combine.plot Logical indicating if all the plot should be combined into one, with same scale (/same axises breaks), same color theme, and same bin size for comparison.
##' @param seed Seed for random number generator. This makes each run easily repeatable. Seed will be automatically assigned if no seed is specified (default). The seed information is stored as an attribute of the returned object. The seed can also be output to a txt file when output=T.
##' @details
##' components analysis uses the likelihood ratio test (LRT) to assess the number of mixture components.
##' Bad Random seed generation may cause normalmixEM to crash. running the function again would be the quickest solution to this issue.
##'
##' @return
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
##' # fit dcoef
##' a=fitNormDistr(dcoef,components=NULL,log.transform=F,combine.plot=F,output=F)
##' # to repeat a
##' b=fitNormDistr(dcoef,components=NULL,log.transform=F,combine.plot=F,output=F,seed=attr(a,"seed"))
##' # if a and b are the same
##' mapply(identical,a[[1]],b[[1]])
##' #try with log transformation
##' c=fitNormDistr(dcoef,components=2,log.transform=T,combine.plot=F,output=F)


##' @export fitNormDistr
###############################################################################
# fit normal distribution to diffusion coefficient

# library(mixtools)

#function for seed setting
.setSeed=function(seed=NULL){
    # set seed
    if (is.null(seed)){
        seed=sample(0:647,1)
        set.seed(seed)
    }else{set.seed(seed)}
    
    note=paste("\nRandom number generation seed",seed,"\n")
    cat(note)
    return(note)
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
    fit.se = list(mu.se= c(fit.info$estimate[2]), sigma.se = c(fit.info$sd[2]), lambda.se = c(0))
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
        
        data.frame(mean=if(log.transform) 10^(result$mu)
                   else result$mu),
        mean.se=as.vector(result.se$mu.se),
        
        
        data.frame(sd=result$sigma,
                   sd.se=result.se$sigma.se),
        
        log.lik=result$loglik
    )
    
    print(class(out))
    print(list(out))
    print(lapply(out, length))
    print(out)
    out=list(out)
    out=t(do.call(cbind.data.frame,out))
    colnames(out) = NULL
    return(out)
}

.fitNormDistr=function(dcoef,components=NULL,log.transform=F,binwidth=NULL,combine.plot=F,output=F){
    auto=FALSE
    # scale=1e3
    name=names(dcoef)
    len=length(dcoef)

    mixmdl.lst=list()
    length(mixmdl.lst)=len
    names(mixmdl.lst)=name

    # store the standard error, currently only for components more than 2
    mixmdl.se.lst=list()
    length(mixmdl.se.lst)=len
    names(mixmdl.se.lst)=name


    for (i in 1:len){

        data=dcoef[[name[i]]][,"slope"]

        # log transformation
        if (log.transform==T) {
            data=log10(data)
            data=data[!is.na(data)]
        }

        # components test
        if (is.null(components)){
            components.int = .getCompNum(data)
            auto=T
        }else{
            components.int=components
        }


        if (components.int==1){#The single component is a special case here since it cannot be a mixEM class object
            oneComp=.singlecompFit(data=data)
            fit.output=oneComp[[1]]
            fit.se=oneComp[[2]]
            plot.mixEM=gg.mixEM(fit.output,binwidth=binwidth,reorder=F)
            
            mixmdl.lst[[i]] = fit.output
            mixmdl.se.lst[[i]] = fit.se
        } else{
            # convergence creteria epsilon = 1e-10 is to filter out false positives
            mixmdl=normalmixEM(data,k=components.int,maxit=1e4,epsilon = 1e-10, lambda = c(0.7,0.3))

            # reorder also for plotting and for output result
            # this is recordered on mean value
            mixmdl=reorderEM(mixmdl)

            print(summary(mixmdl))
            # plot(mixmdl,which=2)
            plot.mixEM=gg.mixEM(mixmdl,binwidth=binwidth,reorder=T)

            mixmdl.lst[[i]]=mixmdl

            # approximate standard error using parametic bootstrap
            cat("\napproximating standard error by parametic bootstrap...\n\n")
            capture.output({mixmdl.se=boot.se(mixmdl, B = 100)})
            # file="/dev/null" # this only for mac

            mixmdl.se.lst[[i]]=mixmdl.se
        }
        # supprress warnings
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

    if (combine.plot==T){

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

        if(components.int==1){reorder=F}
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
    if (output==T){
        if(auto){
            for (i in 1:length(result.lst)){
                result.df=do.call(rbind.data.frame,list(result.lst[[i]]))
                fileName=paste("FitNormDistr-",
                               .timeStamp(name[i]),"....csv",sep="")
                cat("\nOutput FitNormDistr for multiple files\n")
                write.csv(file=fileName,result.df)
            }
        }else{
            result.df=do.call(rbind.data.frame,result.lst)
            fileName=paste("FitNormDistr-",
                       .timeStamp(name[1]),"....csv",sep="")
            cat("\nOutput FitNormDistr.\n")
            write.csv(file=fileName,result.df)
        }
        
        # output plot
        if (combine.plot==T){

            fileName=paste("FitNormDistr-combinePlot-",
                           .timeStamp(name[1]),"....pdf",sep="")
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


fitNormDistr=function(dcoef,components=NULL,log.transform=F,binwidth=NULL,combine.plot=F,output=F,seed=NULL){
    note=.setSeed(seed=seed)

    # output seed
    if (output==T){

        name=names(dcoef)
        fileName=paste("FitNormDistr-",
                       .timeStamp(name[1]),"-Seed....txt",sep="")
        writeLines(text=note,con=fileName)
    }
    
    # return
    structure(.fitNormDistr(dcoef=dcoef,components=components,log.transform=log.transform,binwidth=binwidth,combine.plot=combine.plot,output=output),seed=seed)
}

# TODO:
# 1 check how the se is calculated for the 1 component case and for the multi-component cases
# 2 combine the output of seed into the df