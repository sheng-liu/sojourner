## plotting helpers

##------------------------------------------------------------------------------
## automate binwidth
auto.binwidth=function(x){
    x=x[is.na(x) == F]
    xrange=range(x)
    nbins=median(c(nclass.Sturges(x),nclass.scott(x),nclass.FD(x)))
    binwidth=(xrange[2]-xrange[1])/nbins
    cat("auto binwidth =",binwidth,"\n")
    return(binwidth)
}

##------------------------------------------------------------------------------
## plotHistogram
plotHistogram=function(Log.D.coef,binwidth=0.5, method){
    p=reshape2::melt(Log.D.coef)


    if (method == "static"||method == "percentage"){

        colnames(p)=c("Log.D.coef","file.name")

        # auto binwidth
        if (is.null(binwidth)) binwidth=auto.binwidth(p$Log.D.coef)

        # overlay histogram and density plot without changing count as y axies
        Dcoef.plot=ggplot(p,aes_string(x="Log.D.coef",group="file.name",col="file.name"))+
            geom_histogram(aes_string(y = "..count..",fill="file.name"),
                           binwidth=binwidth,position="dodge")+

            # geom_density(aes(y=0.5*..count..,fill=file.name),alpha=0.2)+
            # geom_density(aes(y=binwidth*..count..,fill=file.name),alpha=0.2)+
            theme_bw()+
            theme(legend.title=element_blank())

        # add visual aid for actual D.coef
        xbreaks=scales::cbreaks(range=c(min(p$Log.D.coef,na.rm=T),
                                max(p$Log.D.coef,na.rm=T)))

        #xbreaks$labels=paste(xbreaks$breaks,10^(xbreaks$breaks),sep="\n")

        lab=paste("(",round(10^(xbreaks$breaks),digits=2),")",sep="")
        xbreaks$labels=paste(xbreaks$breaks,lab,sep="\n")

        Dcoef.plot= Dcoef.plot + scale_x_continuous(breaks=xbreaks$breaks,
                                                    labels=xbreaks$labels)

        plot(Dcoef.plot)
    }else if (method == "rolling.window"){

        colnames(p)=c("Log.D.coef","window.name","file.name")

        # auto binwidth
        if (is.null(binwidth)) binwidth=auto.binwidth(p$Log.D.coef)

        facet.plot=ggplot(p,aes_string(x="Log.D.coef",group="file.name",col="file.name"))+
            geom_histogram(aes_string(y = "..count..",fill="file.name"),
                           binwidth=binwidth,position="dodge")+

            # geom_density(aes(y=0.5*..count..,fill=file.name),alpha=0.2)+
            # this dynamic binwidth does not work
            # geom_density(aes(y=binwidth*..count..,fill=file.name),alpha=0.2)+
            theme_bw()+
            theme(legend.title=element_blank())+
            facet_grid(window.name ~ .)

        merged.plot=ggplot(p,aes_string(x="Log.D.coef",group="file.name",col="file.name"))+
            geom_histogram(aes_string(y = "..count..",fill="file.name"),
                           binwidth=binwidth,position="dodge")+
            geom_density(aes(y=0.5*..count..,fill=file.name),alpha=0.2)+
            theme_bw()+
            theme(legend.title=element_blank())

        multiplot(facet.plot,merged.plot,cols=2)

    }
}
## TODO:
## change the 0.5 to binwidth, so it is dynamic, it is not recognized somehow.


##------------------------------------------------------------------------------
## plotDensity
plotDensity=function(Log.D.coef,binwidth=0.5,method){

    p=reshape2::melt(Log.D.coef)

    if (method == "static"||method == "percentage"){


        colnames(p)=c("Log.D.coef","file.name")

        # auto binwidth
        if (is.null(binwidth)) binwidth=auto.binwidth(p$Log.D.coef)

        Dcoef.plot=ggplot(p,
                          aes_string(x="Log.D.coef",
                              group="file.name",
                              # col=file.name,
                              fill="file.name"))+
            geom_histogram(aes_string(y = "..density..",fill="file.name"),
                           colour="white",
                           binwidth=binwidth,
                           position="dodge")+
            geom_density(aes_string(y = "..density..",col="file.name"),alpha = 0.2)+
            theme_bw()+
            theme(legend.title=element_blank())

        xbreaks=scales::cbreaks(range=c(min(p$Log.D.coef,na.rm=T),
                                        max(p$Log.D.coef,na.rm=T)))

         #xbreaks$labels=paste(xbreaks$breaks,10^(xbreaks$breaks),sep="\n")

        lab=paste("(",round(10^(xbreaks$breaks),digits=2),")",sep="")
        xbreaks$labels=paste(xbreaks$breaks,lab,sep="\n")

        Dcoef.plot= Dcoef.plot + scale_x_continuous(breaks=xbreaks$breaks,
                                                    labels=xbreaks$labels)

        plot(Dcoef.plot)


    }else if (method == "rolling.window"){

        colnames(p)=c("Log.D.coef","window.name","file.name")

        # auto binwidth
        if (is.null(binwidth)) binwidth=auto.binwidth(p$Log.D.coef)

        # a perfect case for faceting
        facet.plot=ggplot(p,
               aes_string(x="Log.D.coef",group="file.name",
                   col="file.name",fill="file.name"))+
            geom_density(alpha = 0.2)+
            theme_bw()+
            theme(legend.title=element_blank())+
            facet_grid(window.name ~ .)

        merged.plot=ggplot(p,
               aes_string(x="Log.D.coef",group="file.name",
                   col="file.name",fill="file.name"))+
            geom_density(alpha = 0.2)+
            theme_bw()+
            theme(legend.title=element_blank())

        ## could also add a merged without 1234

        multiplot(facet.plot,merged.plot,cols=2)


    }


}

# ggplot(p,
#        aes(x=Log.D.coef))+
#     geom_density()

##------------------------------------------------------------------------------
## plotVariance

## not used but keep for now.

plotVariance=function(Log.D.coef,method){

        cat("Generating variance plot \n")

    ## plot data preparation
    ## plot mean of Log.D.coef of each individual trajectory, against variance of each individual trajectory

    # when the list have same length, it maybe easier to work with when converted into data.frame

    Log.D.coef.df=do.call(rbind.data.frame,Log.D.coef)
    MEAN=data.frame(apply(Log.D.coef.df,1,mean,na.rm=T))

    folder=c()
    for (i in 1:dim(MEAN)[1])
        folder[i]=unlist(strsplit(rownames(MEAN)[i],split ="[.]"))[1]
    MEAN=cbind(MEAN,folder)
    colnames(MEAN)=c("mean","folder")


    SD=data.frame(apply(Log.D.coef.df,1,sd,na.rm=T))
    colnames(SD)=c("sd")

    data=cbind(MEAN,SD)

    # plotting
    #     scatter=ggplot(data,aes(x=mean,y=sd,col=folder))+
    #         geom_point(alpha=0.8)+
    #         theme_bw()+
    #         theme(legend.title=element_blank())
    #
    #     mean.density=ggplot(data,aes(x=mean,col=folder))+
    #         geom_density()+
    #         theme_bw()+
    #         theme(legend.title=element_blank())
    #
    #     sd.density=ggplot(data,aes(x=sd,col=folder))+
    #         geom_density()+
    #         theme_bw()+
    #         theme(legend.title=element_blank())
    #
    #     multiplot(scatter,mean.density,sd.density, cols=1)

    # another implementation using gridExtra::grid.arrange
    scatter=ggplot(data,aes(x=mean,y=sd,col=folder))+
        geom_point(alpha=1,shape=21)+
        theme_bw()+
        theme(legend.title=element_blank())+
        theme(legend.position=c(1,1),legend.justification=c(1,1))

    mean.density=ggplot(data,aes(x=mean,fill=folder,col=folder))+
        geom_density(alpha=0.5)+
        theme(legend.title=element_blank())+
        theme(legend.position = "none")

    sd.density=ggplot(data,aes(x=sd,fill=folder,col=folder))+
        coord_flip()+
        geom_density(alpha=0.5)+
        theme(legend.title=element_blank())+
        theme(legend.position = "none")

    empty <- ggplot()+geom_point(aes(1,1), colour="white")+
        theme(
            plot.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank()
        )


    gridExtra::grid.arrange(mean.density, empty, scatter, sd.density, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))



    # this would have worked if mapply takes na.rm=T
    # mapply(mean,Log.D.coef[[1]][[1]],Log.D.coef[[1]][[2]],na.rm=T)


    # or this one
    # mapply(function(x){mean(x,na.rm=T)},Log.D.coef[[1]][[1]],Log.D.coef[[1]][[2]])

    # this means mapply is not taking the elements of each list into one vector, but used them as seperate


    # used alternative, concatanate, then apply
    #             C=mapply("c",Log.D.coef[[1]][[1]],Log.D.coef[[1]][[2]])
    #
    #
    #
    #             apply(C,2,mean,na.rm=T)
    #             mapply(mean,a,b)
    #
    #             lapply(Log.D.coef[[1]])
    #
    #         }
    #
    #         mapply(mean,Log.D.coef[[1]][[1]],Log.D.coef[[1]][[2]],na.rm=T)
    #         mapply(function(x){mean(x,na.rm=T)},
    #                      Log.D.coef[[1]][[1]],Log.D.coef[[1]][[2]])
    #
    #         mapply()
    #
    #         # collapse sublist rolling windowns into uper level list
    #         Log.D.coef=lapply(Log.D.coef,unlist)

}



##------------------------------------------------------------------------------
## .
## from Rcookbook
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    #library(grid)

    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)

    numPlots = length(plots)

    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
        # Make the panel
        # ncol: Number of columns of plots
        # nrow: Number of rows needed, calculated from # of cols
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                         ncol = cols, nrow = ceiling(numPlots/cols))
    }

    if (numPlots == 1) {
        print(plots[[1]])

    } else {
        # Set up the page
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

        # Make each plot, in the correct location
        for (i in 1:numPlots) {
            # Get the i,j matrix positions of the regions that contain this subplot
            matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

            print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                            layout.pos.col = matchidx$col))
        }
    }
}

##------------------------------------------------------------------------------
## .
## from StackOverflow

reorderEM=function(EM){

    order.mu=order(EM$mu)
    EM$mu=EM$mu[order.mu]
    EM$lambda=EM$lambda[order.mu]
    EM$sigma=EM$sigma[order.mu]
    #colnames(EM$posterior)=colnames(EM$posterior)[order.mu]

    # rather than change its name, should subset its content corder
    EM$posterior=EM$posterior[,colnames(EM$posterior)[order.mu]]
    # change to its name
    colnames(EM$posterior)=sort(colnames(EM$posterior))

    return(EM)

}

## the polygon approach
##'@export gg.mixEM
gg.mixEM <- function(EM,binwidth=NULL,reorder=T) {

    # To make multiple ggmixEM plot have the same color theme,  reorder EM's
    # posterior based on mu (how human eyes specify the order), also reorder
    # EM's lambda, mu, sigma, making it consistant within.

#     reorderEM=function(EM){
#
#         order.mu=order(EM$mu)
#         EM$mu=EM$mu[order.mu]
#         EM$lambda=EM$lambda[order.mu]
#         EM$sigma=EM$sigma[order.mu]
#         #colnames(EM$posterior)=colnames(EM$posterior)[order.mu]
#
#         # rather than change its name, should subset its content corder
#         EM$posterior=EM$posterior[,colnames(EM$posterior)[order.mu]]
#         # change to its name
#         colnames(EM$posterior)=sort(colnames(EM$posterior))
#
#         return(EM)
#
#     }

    if (reorder == T) EM=reorderEM(EM)

    # reconstruct x, may use sample
    x       <- with(EM,seq(min(x),max(x),len=1000))
    # parameters holder
    pars    <- with(EM,data.frame(comp=colnames(posterior), mu, sigma,lambda))

    em.df   <- data.frame(x=rep(x,each=nrow(pars)),pars)

    # reconstruct normal distribution with parameters
    em.df$y <- with(em.df,lambda*dnorm(x,mean=mu,sd=sigma))


    # auto binwidth
    if (is.null(binwidth)) binwidth=auto.binwidth(EM$x)

    ggplot(data.frame(x=EM$x),aes_string(x="x",y="..density..")) +
        geom_histogram(fill=NA,color="black",binwidth=binwidth)+
        # when distribution is truncated, it plots flippers
        # geom_polygon(data=em.df,aes(x,y,fill=comp),color="grey50", alpha=0.5)+
        geom_area(data=em.df,aes_string(x="x",y="y",fill="comp"),color="grey50", alpha=0.5,position = "identity")+
        scale_fill_discrete("Component\nMeans",labels=format(em.df$mu,digits=3))+
        theme_bw()
}

# gg.mixEM(EM)
