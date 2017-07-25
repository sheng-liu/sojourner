## fitNormDistr-methods
##
##
###############################################################################
##' @name fitNormDistr
##' @aliases fitNormDistr
##' @title fitNormDistr
##' @rdname fitNormDistr-methods
##' @docType methods
##' @description fit normal distribution to diffusion coefficient caclulated by Dcoef method.
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
##'@details
##'components analysis uses the likelihood ratio test (LRT) to assess the number of mixture components.
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
##' folder1=system.file("extdata","SWR1",package="smt")
##' folder2=system.file("extdata","HTZ1",package="smt")
##' trackll=compareFolder(c(folder1,folder2))
##' MSD=msd(trackll=trackll)
##' dcoef=Dcoef(MSD,dt=6,plot=T,output=F)
##' # fit dcoef
##' a=fitNormDistr(dcoef,components=NULL,log.transform=F,combine.plot=F,output=F)
##' # to repeat a
##' b=fitNormDistr(dcoef,components=NULL,log.transform=F,combine.plot=F,output=F,seed=attr(a,"seed"))
##' # if a and b are the same
##' mapply(identical,a[[1]],b[[1]])


##' @export fitNormDistr
###############################################################################
# fit normal distribution to diffusion coefficient

# library(mixtools)

.fitNormDistr=function(dcoef,components=NULL,log.transform=F,binwidth=NULL,combine.plot=F,output=F){

    # random.seed=NULL
    # set seed for random number generation
    # according to stata How to choose a seed
#     rseed= if (is.null(random.seed)) sample(0:647,1) else random.seed
#     cat("\nset.seed (",rseed, ") for random number generation\n")

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

            cat("\ncomponents analysis\n")
            # verbose=T, add progress bar, it takes long
            components.test=mclust::mclustBootstrapLRT(data,model = "V",verbose=T)
            print(components.test)

            components.int=length(components.test$p.value)
            cat("\n\nmost likely components",components.int,"at significant level 0.05\n\n")
        }else{
            components.int=components
        }


        if (components.int==1){

            cat("\nSMT v0.3 currently does not deal with fitting one component,\ncomponents sets to 2 automatically.\n")
            components.int=2


            # this fit ignores the negative values
            # data=data[data>=0]
            # mixmdl=fitdist(data,"norm")
            # mixmdl=fitdist(data*scale,"norm",method="mle")
            # small values need scaling, mme doesn't and generate the same fit.

            # these code needs to reformat to the same as components==2
#             mixmdl=fitdistrplus::fitdist(data,"norm",method="mme")
#             print(summary(mixmdl))
#
#             # denscomp(fitg,demp=T) or plot(mixmdl)
#             fitdistrplus::denscomp(mixmdl, addlegend=FALSE)
#             mixmdl.lst[[i]]=mixmdl

        }

            # convergence creteria epsilon = 1e-10 is to filter out false positives
            mixmdl=normalmixEM(data,k=components.int,maxit=1e4,epsilon = 1e-10)

            # reorder also for plotting and for output result
            # this is recordered on mean value
            mixmdl=reorderEM(mixmdl)

            print(summary(mixmdl))
            # plot(mixmdl,which=2)
            plot.mixEM=gg.mixEM(mixmdl,binwidth=binwidth,reorder=T)
            # supprress warnings
            # `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
                suppressMessages(plot(plot.mixEM))


            # mixmdl[c("mu","sigma","lambda")]
            mixmdl.lst[[i]]=mixmdl

            # approximate standard error using parametic bootstrap
            cat("\napproximating standard error by parametic bootstrap...\n\n")
            capture.output({mixmdl.se=boot.se(mixmdl, B = 100)})
            # file="/dev/null" # this only for mac

            mixmdl.se.lst[[i]]=mixmdl.se


    }


    # abstract result from the fitting model
    if (components.int==1){

        result.lst=lapply(mixmdl.lst,function(x){
            s=summary(x)
            result=list(mean=s$estimate[1],sd=s$estimate[2])
            if (log.transform) result$mean=10^(result$mean)
            result=do.call(cbind.data.frame,result)
            print(result)
            return(result)
        })

    }else{


        #         result.lst=lapply(mixmdl.lst,
        #                          function(x){
        #                              result=list(proportion=x$lambda,
        #                                          mean=x$mu,
        #                                          sd=x$sigma,
        #                                          log.lik=x$loglik)
        #                              if (log.transform) {
        #                                  result$mean=10^(result$mean)
        #
        #                                  }
        #                              result=do.call(cbind.data.frame,result)
        #                              return(result)})


        #names(result)[which(names(result)=="sd")]="cv"

        result.lst=list()
        length(result.lst)=len
        names(result.lst)=name

        for (i in 1:len){

            result=list(
                data.frame(proportion=mixmdl.lst[[i]]$lambda,
                           proportion.se=mixmdl.se.lst[[i]]$lambda.se),

                #  data.frame(mean=ifelse(log.transform,
                #               10^(mixmdl.lst[[i]]$mu),mixmdl.lst[[i]]$mu),

                data.frame(mean=if(log.transform) 10^(mixmdl.lst[[i]]$mu)
                           else mixmdl.lst[[i]]$mu),
                mean.se=as.vector(mixmdl.se.lst[[i]]$mu.se),


                data.frame(sd=mixmdl.lst[[i]]$sigma,
                           sd.se=mixmdl.se.lst[[i]]$sigma.se),

                log.lik=mixmdl.lst[[i]]$loglik
            )

            result=t(do.call(cbind.data.frame,result))
            result.lst[[i]]=result
        }

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


        plot.lst=lapply(mixmdl.lst,function(x){

            p=gg.mixEM(x,binwidth=binwidth,reorder=T)+
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

        result.df=do.call(rbind.data.frame,result.lst)
        fileName=paste("FitNormDistr-",
                       .timeStamp(name[1]),"....csv",sep="")
        cat("\nOutput FitNormDistr.\n")
        write.csv(file=fileName,result.df)

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

    # set seed
    if (is.null(seed)){
        seed=sample(0:647,1)
        set.seed(seed)
    }else{set.seed(seed)}

    note <- paste("\nRandom number generation seed",seed,"\n")
    cat(note)

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
# 1 components=1 format
# 2 combine the output of seed into the df


