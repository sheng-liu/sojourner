
## plotTrack-methods
##
##
###############################################################################
##' @name plotTrack
##' @aliases plotTrack plotTrackFromIndex plotTrackOverlay plotNucTrackOverlay 
##' plotComponentTrackOverlay plotMask trackOverlayData
##' @title plotTrack
##' @rdname plotTrack-methods
##' @docType methods
##' @description Plot track/trajectory from track list. either randomly or
##'   specified.

##' @usage
##'
##' plotTrack(ab.trackll,resolution=0.107,frame.min=8,frame.max=100,
##' frame.start=1,frame.end=500)
##'
##' plotTrackFromIndex(index.file, movie.folder,
##' resolution=0.107,frame.min=1,frame.max=100,
##' frame.start=1,frame.end=500,input=1)
##'
##' plotTrackOverlay(trackll,max.pixel=128,nrow=2,ncol=2,width=16,height=16)
##'
##' plotNucTrackOverlay(folder,trackll=NULL,cores=1,
##' max.pixel=128,nrow=2,ncol=2,width=16,height=16)
##'
##' plotComponentTrackOverlay(folder,trackll.sel=NULL,
##' max.pixel=128,nrow=2,ncol=2,width=16,height=16)
##'
##' plotMask(folder,max.pixel=128,nrow=2,ncol=2,width=16,height=16)
##'
##' trackOverlayData(trackl)
##'
##' @param ab.trackll absolute coordinates for plotting, generated from
##'   readDiatrack(folder,ab.track=T).
##' @param trackl Track list
##' @param trackll Track list output from readDiatrack().
##' @param folder folder containing desired input data.
##' @param resolution ratio of pixel to uM.
##' @param frame.min minimum frame number for plotting.
##' @param frame.max max frame number for plotting.
##' @param frame.start the first frame to plot. Default 1.
##' @param frame.end last frame to plot. Default 500.
##' @param index.file a csv file that contains index of tracks in the first
##'   column. Leave a header line when preparing such a file.
##' @param movie.folder the path to the folder which contains Diatrack output
##'   txt files (presumably it is the same folder with movie files).
##' @param max.pixel Number of pixels of imaging regime.
##' @param nrow Number of rows in the final plot.
##' @param ncol Number of colums in the final plot.
##' @param width Width of the page for plotting.
##' @param height Height of the page for plotting.
##' @param cores Number of cores to be used.
##' @param trackll.sel Selected component trajectory output by
##'   selComponentTracks().
##'
##' @param input Input file type (Diatrack .txt file = 1; Diatrack .mat session
##'   file = 2; ImageJ .csv file = 3; SlimFast .txt file = 4).

##' @return
##' \itemize{

##' \item{PDF} One PDF file with all the frames satisfy the creteria. If trackll
##' has multiple items, it ouptus mutiple PDF files each corresponding to one
##' item.
##'
##' \item{csv} Outputs csv file of the coordiantes of the trajectory, which
##' users can use other plotting software (e.g. Prism or Excel) to plot tracks
##' in their favor.

##' }
##' @details
##' \itemize{
##' \item{plotTrackFromIndex:} if user provide a csv file with first column
##' listing the index of trajectories, this program will plot the tracks isted
##' in the csv file. It is useful after manipulating with the output from Dceof,
##' to plot the tracks that of interest to the user (e.g. highest Dcoef). User
##' need to provide the indexFile.csv, and specify the movie folder which
##' contains the movies where specified trajectories are tracked.
##'
##' \item{plotTrackOverlay:} plot all tracks in trackll overlaid on one plot.
##'
##' \item{plotNucTrackOverlay:} plot tracks in a movie overlayed with nuclei
##' image. The nuclei image file must end with "_Nuclei.tif" to be recognized.
##' If trackll is NULL (default), program will read in trackll from specified
##' folder and return trackll, otherwise it will take the specified trackll
##' directly.
##'
##' \item{plotComponentTrackOverlay:} plot tracks base on component fitting of
##' diffusion coefficient. Combined with selComponentTracks() function, together
##' it allows select and plot tracks based on component fitting of track
##' diffusion coefficient.
##'
##' \item{plotMask:} plot image mask. The mask file name must ended with
##' _MASK.tiff to be recognized.
##'
##' }

##' @examples
##' folder=system.file("extdata","SWR1",package="sojourner")
##' trackll.ab=readDiatrack(folder,ab.track=TRUE)
##' plotTrack(trackll.ab)
##'
##' ## plot from index file
##' index.file=system.file("extdata","INDEX","indexFile.csv",package="sojourner")
##' movie.folder=system.file("extdata","SWR1",package="sojourner")
##' plotTrackFromIndex(index.file=index.file,movie.folder = movie.folder)
##'
##' ## index file contain trajectories from multiple movie folders
##' folder1=system.file("extdata","SWR1",package="sojourner")
##' folder2=system.file("extdata","HTZ1",package="sojourner")
##' index.file2=system.file("extdata","INDEX","indexFile2.csv",package="sojourner")
##' plotTrackFromIndex(index.file=index.file2,movie.folder = c(folder1,folder2))
##'
##' ## masking with image mask
##' track.folder=system.file("extdata","SWR1_2",package="sojourner")
##' trackll=readDiatrack(folder=track.folder)
##' trackll.masked=readDiatrack(folder=track.folder)
##' str(trackll,1)
##' str(trackll.masked,1)
##'
##' ## compare the masking effect
##' plotTrackOverlay(trackll,nrow=1,ncol=1,width=8,height=8)
##' plotTrackOverlay(trackll.masked,nrow=1,ncol=1,width=8,height=8)
##'
##' ## compare masking effect with nuclei image
##' plotNucTrackOverlay(folder=track.folder,trackll,
##'                     nrow=1,ncol=1,width=8,height=8)
##' plotNucTrackOverlay(folder=track.folder,trackll.masked,
##'                     nrow=1,ncol=1,width=8,height=8)
##'
##' ## plot mask
##' plotMask(track.folder,nrow=1,ncol=1,width=8,height=8)
##'
##'
##' ## plotComponentTrackOverlay (see selComponentTracks() for more details)
##' folder2=system.file("extdata","SWR1_2",package="sojourner")
##' trackll=readDiatrack(folder2)
##'
##' ## use merge=T for per folder comparison, the analsyis result can't be plot

##' ##back to original image. To see component tracks on original nuclei image,
##' ##set merge=F, for per movie analysis.

##'
##' ## compute MSD
##' MSD=msd(trackll=trackll,plot=T)
##' msd(trackll=trackll,summarize=T,plot=T)
##'
##' ## calculate Dcoef
##' dcoef=Dcoef(MSD=MSD,method="static",plot=TRUE)
##'
##' ## fit normal distribution to define component
##' ## set seed to reproduce results (see fitNormalDistr() for details on seed)
##' fit=fitNormDistr(
##'        dcoef,components=2,log.transform=T,combine.plot=F,output=F,seed=481)
##'
##' ## select component tracks based on fitting
##' trackll.sel=selComponentTracks(trackll,
##'       fit=fit,likelihood = 0.9,dcoef = dcoef,log.transformed = T,output = F)
##'
##' ## plot component tracks
##' plotComponentTrackOverlay(folder2,trackll.sel=trackll.sel)
##'

## @import reshape2
##' @export plotTrack
##' @export .plotTrack
##' @export plotTrackFromIndex
## @import animation
## FUTURE: maybe plot on dt

##' @export plotTrackOverlay
##' @export plotMask
##' @export plotNucTrackOverlay
##' @export plotComponentTrackOverlay

## .plotMask(mask.list[1])

## FOR SHINY USAGE:
##' @export trackOverlayData
##' @export .plotNucTrackOverlay



##
## useful when the file is too big to include into the package.
## masking with image mask
## Not run:
## download Diatrack output txt file and mask
## save file to a location on local computer corresponding to your operation system

# mask.url="https://www.dropbox.com/s/6472kweldoh96xd/SWR1_halo_140mW3_20151025_SWR1WT_Substack%20%28400-5800%29_MASK.tif?dl=0"
# data.url="https://www.dropbox.com/s/lav6f0y4gd4j4lt/SWR1_halo_140mW3_20151025_SWR1WT_Substack%20%28400-5800%29.txt?dl=0"
# nuclei.url="https://www.dropbox.com/s/t6erxm3wze9pz0r/Temp_Placeholder_Nuclei.tif?dl=0"

# dir.create("~/masking_test/")
# download.file(mask.url, "~/masking_test/_MASK.tif")
# download.file(data.url, "~/masking_test/_DATA.txt")
# download.file(nuclei.url, "~/masking_test/_Nuclei.tif")
# track.folder="~/masking_test/"



###############################################################################

.plotTrack=function(ab.trackl,file.name="TrajectoryPlot",
                    resolution=0.107,frame.min=8,frame.max=100,frame.start=1,frame.end=500){
    
    # trackl is just a list of trajectories, with no upper level indicating folder
    ab.trackl.res=lapply(ab.trackl,function(x) x*resolution)
    m=max(sapply(ab.trackl.res,max))
    
    
    # add names to each plot
    name=names(ab.trackl)
    
    
    fileName=paste(.timeStamp(file.name),".pdf",sep="")
    
    pdf(file=fileName)
    
    par(mfrow=c(2,2))
    len=length(ab.trackl.res)
    
    frame.end=ifelse(len<frame.end,len,frame.end)
    
    for (i in frame.start:frame.end){
        
        p=ab.trackl.res[[i]]
        frame.len=dim(p)[1]
        if (frame.len>frame.min & frame.len<frame.max)
            plot(p$x,p$y,type="l",xlim=c(0,m),ylim=c(0,m),xlab="X (uM)",
                 ylab="Y (uM)",main=name[[i]])
    }
    
    # sub = name[[i]]
    dev.off()
    
    return(ab.trackl.res)
    
}

## TODO: this is good and simple, to be consistant with the plotting scheme, making it ggplot2 style if not require enormous amount of work.


##------------------------------------------------------------------------------
##

plotTrack=function(ab.trackll,resolution=0.107,
                   frame.min=8,frame.max=100,frame.start=1,frame.end=500){
    
    file.name=names(ab.trackll)
    
    
    for (i in 1:length(file.name)){
        
        # output plot
        cat("\nOutput track plot...\n")
        
        plot.coords=.plotTrack(
            ab.trackll[[i]],file.name[i],resolution=resolution,
            frame.min=frame.min,frame.max=frame.max,frame.start=frame.start,
            frame.end=frame.end)
        
        # output csv of the plot
        cat("\nOutput csv file for track plot...\n")
        
        plot.coords.df=do.call(rbind,plot.coords)
        # melt(plot.coords)
        
        
        fileName=paste("Track Coords-",.timeStamp(file.name[i]),".csv",sep="")
        write.csv(file=fileName,plot.coords.df)
        
    }
    
    # lapply(ab.trackll,function(ab.trackl,file.name){.plotTrack(ab.trackl,file.name,resolution=resolution,frame.min=frame.min,frame.max=frame.max,frame.start=frame.start,frame.end=frame.end)})
    # lapply can only take one input
    
}

##------------------------------------------------------------------------------
##

## plot trajectory according to index
## user need to put the corresponding movie files into a folder

plotTrackFromIndex=function(index.file, movie.folder,resolution=0.107,
                            frame.min=1,frame.max=100,frame.start=1,
                            frame.end=500,input=1){
    
    ## read trajectory index from the index.file
    index.df=read.csv(file=index.file,header=T)
    index=as.character(index.df[,1])
    
    
    #     # DONE: the number of folder to compare can be extended using ... statement
    
    #     folder.list=list()
    #     for (i in 1:length(movie.folder)){
    #         folder.list[[i]]=movie.folder[i]
    #     }
    #
    #
    #     # remove null folders by subsetting un-null folders
    #     null.folder=sapply(folder.list,is.null)
    #     folder.list=folder.list[!null.folder]
    #
    #     names(folder.list)=sapply(folder.list,basename)
    #
    #     ab.trackll=list()
    #
    #     for (i in 1:length(folder.list)) {
    #
    #         ## read in tracks in movie.folder with absolute coords,
    #         ## merge them as the input is merged csv files
    #         ab.trackll[i]=readDiatrack(folder=folder.list[[i]],merg=T,ab.track=T)
    #         cat("\n...\n") # seperator makes ouput clearer
    #         names(ab.trackll)[i]=names(folder.list)[i]
    #
    #
    #     }
    
    ab.trackll=compareFolder(folders=movie.folder,ab.track=T,input=input)
    
    # as it is only for one folder
    # trackl.plot=ab.trackll[[1]][index]
    
    ## because it is a merged list, use .plotTrack
    ## .plotTrack(trackl.plot,resolution=resolution,frame.min=frame.min,
    ## frame.max=frame.max,frame.start=frame.start,frame.end=frame.end)
    
    
    # or one can maintain the ab.trackl's structure, which has the folder name
    # because of list, this subsetting using index become amazingly simple
    # no further manipulation or matching needed. no assembly required.
    trackll.plot=lapply(ab.trackll,function(x){x[index]})
    
    
    
    
    ## remove NA (or combine it into one list)
    trackll.plot.narm=lapply(trackll.plot, function(x) x[!is.na(names(x))])
    
    ## if the list is length 0, remove it
    ## remove zero length list
    for (i in 1:length(trackll.plot.narm)){
        
        if (length(trackll.plot.narm[[i]]) == 0) trackll.plot.narm[[i]]=NULL
    }
    
    plotTrack(trackll.plot.narm,resolution=resolution,frame.min=frame.min,
              frame.max=frame.max,frame.start=frame.start,frame.end=frame.end)
    
    
}

# return(trackll.plot)
##------------------------------------------------------------------------------
##

# plot all movies in a folder

# becasue "indexPerTrackll"
# must be list from one movie
# TrackOverlay only plots the first object in trackll list,
# use trackll[n] to specifically nth object in the trackll list

# collaps tracks of a single file (i.e. trackl)
trackOverlayData=function(trackl){
    
    
    # plot.title=names(trackl)
    
    cat("\nProcessing",names(trackl))
    
    track.df=do.call(rbind.data.frame,trackl[[1]])
    
    # must use [[1]]
    #     > str(trackll[[1]],0)
    #     List of 908
    #     [list output truncated]
    #     > str(trackll[1],0)
    #     List of 1
    
    
    # split rownames
    n=track.df
    if (length(grep("txt",rownames(n)[1])) == 0){
        Index=strsplit(rownames(n),"\\.")
    }else{
        Index=strsplit(rownames(n),".txt.")
    }
    
    # trackID=fileID.frameID.duration.indexPerFile.indexPerTrackll
    Index.df=data.frame(do.call(rbind,Index))
    #do.call(rbind.data.frame,Index)
    
    colnames(Index.df)=c(
        "fileID","frameID","duration","indexPerFile","indexPerTrackll","SN")
    
    track.plot.df=cbind(track.df,Index.df)
    return(track.plot.df)
    
}

# for a single file (i.e. trackl)


.plotTrackOverlay=function(trackl,max.pixel=128){
    
    
    # get names of the trackll (/video) to put it on each graph
    plot.title=names(trackl)
    
    
    track.overlay.data=trackOverlayData(trackl)
    
    # note group needs to be on two /three variable (i.e. frameID and
    # indexPerFile) as indexPerFile sometimes repeat when multiple files are
    # merged. interaction() realize it
    
    track.overlay.data$inter=interaction(
        track.overlay.data$frameID,track.overlay.data$indexPerFile)
    
    p=ggplot(track.overlay.data,
             aes_string(x="x",y="y",
                        # group=interaction("frameID","duration","indexPerFile")))+
                        group="inter"))+
        geom_path()+
        
        #     > interaction("frameID","duration","indexPerFile")
        # [1] frameID.duration.indexPerFile
        # Levels: frameID.duration.indexPerFile
        
        # using interaction within aes_string (to pass RCMD check)
        # what you want is a column of actual values, not the name
        # the data has to be processed beforhand
        # https://stackoverflow.com/questions/19410781/problems-when-using-ggplot-aes-string-group-and-linetype/19415464#19415464
        
        scale_x_continuous(
            name="Pixel",
            breaks=seq(from=0, to=max.pixel,by=20),
            limits=c(0,max.pixel))+
        scale_y_continuous(
            name="Pixel",
            breaks=seq(from=0, to=max.pixel,by=20),
            limits=c(0,max.pixel))+
        
        ggtitle(plot.title)+
        
        # this makes integer breaks
        #        scale_x_continuous(breaks=scales::pretty_breaks(n=5))+
        #        scale_y_continuous(breaks=scales::pretty_breaks(n=5))+
        #labs(x="Pixel", y="Pixel")+
        
        theme_bw()
    
    return(p)
    
}


# mask and unmask needs to be done before plotting to save time on reading in data
# plot multiple trackl (videos)

plotTrackOverlay=function(trackll,max.pixel=128,nrow=2,ncol=2,width=16,height=16){
    
    # get plot.lst
    plot.lst=list()
    
    for (i in 1:length(trackll)) plot.lst[[i]]=.plotTrackOverlay(
        trackl=trackll[i],max.pixel=max.pixel)
    
    # output
    cat("\nOutput combined plot...")
    cmb.plot=gridExtra::marrangeGrob(plot.lst,nrow=nrow,ncol=ncol)
    
    fileName=paste(.timeStamp("TrackOverlay"),".pdf",sep="")
    # TODO: add folder name in .timeStamp
    # fileName=paste("TrackOverlay-",.timeStamp("folder"),".pdf",sep="")
    
    ggplot2::ggsave(filename=fileName,cmb.plot,width=width,height=height)
    
    cat("\nDone!")
    
}


##------------------------------------------------------------------------------
##

# Plot mask
.plotMask=function(mask.file,max.pixel=128){
    
    title=basename(mask.file)
    # read in tiff mask
    # library(rtiff)
    cat("\nReading mask file",title)
    mask=rtiff::readTiff(fn=mask.file)
    # plot(mask)
    
    # TODO 
    # EBImage::readImage(files)
    
    
    # img=EBImage::readImage(image.file)
    # d=img@.Data
    
    
    
    # 
    pospt=which(mask@red!=0,arr.ind=T)
    pos.point=with(data.frame(pospt),data.frame(x=col,y=row))
    
    # horizontal is the same vertical is fliped as the pixels is counted from
    # upper left in the image, but is counted from lower left in the plot.
    
    # shape 22 is little square, when squeezed, they have misconsumption of
    # smaller squres
    
    # shape 46 is the smallest dot, even when squeezed
    
    pp=ggplot(pos.point,aes_string(x="x",y="y"))+geom_point(alpha=1,shape=22)+
        scale_x_continuous(
            name="Pixel",
            breaks=seq(from=0, to=max.pixel,by=20),
            limits=c(0,max.pixel))+
        scale_y_continuous(
            name="Pixel",
            breaks=seq(from=0, to=max.pixel,by=20),
            limits=c(0,max.pixel))+
        ggtitle(title)+
        theme_bw()
    
    plot(pp)
    
    
    
    # return plot for better manipulation downstream
    # return(invisible(pos.point))
    return(pp)
    
}


plotMask=function(folder,max.pixel=128,nrow=2,ncol=2,width=16,height=16){
    
    mask.lst=list.files(path=folder,pattern="_MASK.tif",full.names=T)
    
    mask.plot.lst=lapply(mask.lst,.plotMask,max.pixel=max.pixel)
    
    cmb.plot.mask=gridExtra::marrangeGrob(mask.plot.lst,nrow=nrow,ncol=ncol)
    
    fileName=paste(.timeStamp("TrackMask"),".pdf",sep="")
    
    ggsave(filename=fileName,cmb.plot.mask,width=width,height=height)
    
    cat("\nDone!")
}

##------------------------------------------------------------------------------
##

# .plotNucTrackOverlay

# color="red"
# color=track.overlay.data$component

# process one movie at a time, trackl or component.trackl (which corresponding to
# one movie with two components)

.plotNucTrackOverlay=function(trackl=NULL,component.trackl=NULL,image.file,
                              max.pixel=128,color="red"){
    
    # double input, trackl or component.trackl to allow more flexible
    # manipulation of data before pass in to plot. trackl can be null,
    # component.trackl must have value
    
    img=EBImage::readImage(image.file)
    d=img@.Data
    
    # 2018-09-29
    # when trackl, color is charactor, aes_string(... shQuote(color))
    # when component.trackl, color become "factors", aes_string(... color)
    # when use aes, it reads both, when use aes_string, it read it as string.
    
    p=  ggplot()+
        geom_raster(data=reshape2::melt(d), 
                    aes_string(x="Var1",y="Var2",fill="value"),interpolate=FALSE)+
        scale_fill_gradient(low = "black", high = "white")+ guides(fill=FALSE)
    
    
    if (
        (is.null(component.trackl)&&is.null(trackl))||
        (!is.null(component.trackl)&&(!is.null(trackl)))
    ){stop(
        "Please provide either trackl or component.trackl")} 
    
    
    # if pass in trackl, i.e. color is character "red"
    
    # alternative check on the class of color
    # if (!is.factor(color))
    
    if (!is.null(trackl)){
        # get names of the trackll (/video) to put it on each graph the good
        # about trackl over data.frame component.trackl is it also passes in name
        # change trackOverlayData() to take in trackl
        plot.title=names(trackl)
        track.overlay.data=trackOverlayData(trackl)
        
        
        # p=  ggplot()+
        #     geom_raster(data=reshape2::melt(d), 
        #                 aes_string(x="Var1",y="Var2",fill="value"),interpolate=FALSE)+
        #     scale_fill_gradient(low = "black", high = "white")+ guides(fill=FALSE)+
        
        # ggplot()+
        
        p=p+
            geom_path(data=track.overlay.data,
                      aes_string(x="x",y="y",group="indexPerFile",color=shQuote(color)))
        ##aes_string(x="x",y="y",group="indexPerFile",color=color))+
        
        
    }else{
        
        # if pass in component.trackl, i.e. color is factor
        
        #plot.title=""
        plot.title=names(component.trackl)
        track.overlay.data=cmpOverlayData(component.trackl)
        
        
        p=p+
            geom_path(data=track.overlay.data,
                      # aes_string(x="x",y="y",group="indexPerFile",color=shQuote(color)))+
                      aes_string(x="x",y="y",group="indexPerFile",color=color))
        
    }
    
    
    p=p+
        scale_x_continuous(
            name="Pixel",
            breaks=seq(from=0, to=max.pixel,by=20),
            limits=c(0,max.pixel),
            # remove the blank between image and axies
            expand = c(0, 0))+
        scale_y_continuous(
            name="Pixel",
            breaks=seq(from=0, to=max.pixel,by=20),
            limits=c(0,max.pixel),
            # remove the blank between image and axies
            expand = c(0, 0))+
        ggtitle(plot.title)+
        
        # this makes integer breaks
        #        scale_x_continuous(breaks=scales::pretty_breaks(n=5))+
        #        scale_y_continuous(breaks=scales::pretty_breaks(n=5))+
        #labs(x="Pixel", y="Pixel")
        
        theme_bw()+
        theme(legend.title=element_blank())
    
    # removes line and text of axis, only picture
    #     +theme(line=element_blank(),
    #            text=element_blank())
    
    
    
    
    #     p=  ggplot()+
    #         geom_raster(data=reshape2::melt(d), 
    #                     aes_string(x="Var1",y="Var2",fill="value"),interpolate=FALSE)+
    #         scale_fill_gradient(low = "black", high = "white")+ guides(fill=FALSE)+
    # 
    #         # ggplot()+
    #         geom_path(data=track.overlay.data,
    #                   # aes_string(x="x",y="y",group="indexPerFile",color=shQuote(color)))+
    #                   aes_string(x="x",y="y",group="indexPerFile",color=color))+
    #         
    #         scale_x_continuous(
    #             name="Pixel",
    #             breaks=seq(from=0, to=max.pixel,by=20),
    #             limits=c(0,max.pixel),
    #             # remove the blank between image and axies
    #             expand = c(0, 0))+
    #         scale_y_continuous(
    #             name="Pixel",
    #             breaks=seq(from=0, to=max.pixel,by=20),
    #             limits=c(0,max.pixel),
    #             # remove the blank between image and axies
    #             expand = c(0, 0))+
    #         ggtitle(plot.title)+
    # 
    #         # this makes integer breaks
    #         #        scale_x_continuous(breaks=scales::pretty_breaks(n=5))+
    #         #        scale_y_continuous(breaks=scales::pretty_breaks(n=5))+
    #         #labs(x="Pixel", y="Pixel")
    # 
    #         theme_bw()+
    #         theme(legend.title=element_blank())
    # 
    #     # removes line and text of axis, only picture
    # #     +theme(line=element_blank(),
    # #            text=element_blank())
    
    if (!is.factor(color)) {p=p+theme(legend.position="none")}
    plot(p)
    
    return(p)
}

##------------------------------------------------------------------------------
##

# plotNucTrackOverlay

plotNucTrackOverlay=function(folder,trackll=NULL,cores=1,
                             max.pixel=128,
                             nrow=2,ncol=2,width=16,height=16){
    
    nuclei.lst=list.files(path=folder,pattern="_Nuclei.tif",full.names=T)
    
    if (is.null(trackll)){
        cat("\ntrackll not specified, read in Diatrack file\n")
        trackll=readDiatrack(folder=folder,cores=cores)
    }
    
    
    plot.lst=list()
    for (i in 1:length(trackll)) plot.lst[[i]]=
        
        suppressWarnings(
        .plotNucTrackOverlay(
        trackl=trackll[i],component.trackl=NULL,image.file=nuclei.lst[[i]],
        max.pixel=max.pixel,color="red")
        )
    
    # output
    cat("\nOutput combined plot...")
    cmb.plot=gridExtra::marrangeGrob(plot.lst,nrow=nrow,ncol=ncol)
    
    fileName=paste(.timeStamp("NucTrackOverlay"),".pdf",sep="")
    # TODO: add folder name in .timeStamp
    # fileName=paste("TrackOverlay-",.timeStamp("folder"),".pdf",sep="")
    
    ggplot2::ggsave(filename=fileName,cmb.plot,width=width,height=height)
    # tip: save as png help keep the raster pixel better than pdf
    
    cat("\nDone!")
    
    return(invisible(trackll))
    
}

##------------------------------------------------------------------------------
##

# plotComponentTrackOverlay

plotComponentTrackOverlay=function(folder,trackll.sel=NULL,
                                   max.pixel=128,
                                   nrow=2,ncol=2,width=16,height=16){
    
    nuclei.lst=list.files(path=folder,pattern="_Nuclei.tif",full.names=T)
    
    # no need to read them in, as only plot trackll.sel to the images
    # the folder is only to provide image location
    
    # if (is.null(trackll)){
    #     cat("\ntrackll not specified, read in Diatrack file\n")
    #     trackll=readDiatrack(folder=folder,merge=F,mask=mask,cores=cores)
    # }
    
    plot.lst=list()
    
    # track.overlay.data.lst=lapply(trackll.sel,function(x){cmpOverlayData(x)})
    
    track.overlay.data.lst=list()
    length(track.overlay.data.lst)=length(trackll.sel)
    names(track.overlay.data.lst)=names(trackll.sel)
    
    for (i in 1:length(trackll.sel)){
        
        track.overlay.data.lst[[i]]=cmpOverlayData(trackll.sel[i])
    }
    
    
    color.lst=lapply(track.overlay.data.lst,function(x){x$component})
    
    # trackll or trackll.sel 
    # plotComponentTrackOverlay: no visible binding for global variable "trackll"
    # print(track.overlay.data.lst)
    # print(color.lst)
    # print(class(color.lst[[1]]))
    
    ### trackll -> trackll.sel
    for (i in 1:length(trackll.sel)) plot.lst[[i]]=
        
        # Warning message:
        # In readTIFF(x, all = all, ...) :
        # TIFFReadDirectory: Unknown field with tag 65531 (0xfffb) encountered
        # these are unnecesary private tags can be suppressed
        # https://stackoverflow.com/questions/27608124/imagemagick-how-to-get-rid-of-tiffwarnings-768-message-about-unknown-field-wh
        
        suppressWarnings(.plotNucTrackOverlay(
            trackl=NULL,component.trackl=trackll.sel[i],
            image.file=nuclei.lst[[i]],
            max.pixel=max.pixel,color=color.lst[[i]]))
    
    
    # output
    cat("\nOutput combined plot...")
    cmb.plot=gridExtra::marrangeGrob(plot.lst,nrow=nrow,ncol=ncol)
    
    fileName=paste(.timeStamp("NucTrackOverlay"),".pdf",sep="")
    # TODO: add folder name in .timeStamp
    # fileName=paste("TrackOverlay-",.timeStamp("folder"),".pdf",sep="")
    
    ggplot2::ggsave(filename=fileName,cmb.plot,width=width,height=height)
    # tip: save as png help keep the raster pixel better than pdf
    
    cat("\nDone!")
    
    return(invisible(trackll.sel))
    
}


##------------------------------------------------------------------------------
##

# cmpOverlayData

# pass in component.trackl with name, rather than just comp1 comp2
# str(trackll.sel[1],2)
# $ 120mW_10ms1.txt:List of 2
# ..$ comp.1:List of 12
# ..$ comp.2:List of 26

# ratehr than
# > str(trackll.sel[[1]],1)
# List of 2
# $ comp.1:List of 12
# $ comp.2:List of 26

# component.trackl=trackll.sel[1]
# > str(component.trackl,2)
# List of 1
# $ 120mW_10ms1.txt:List of 2
# ..$ comp.1:List of 12
# ..$ comp.2:List of 26

# used to be pass in two component list
# now pass in one list of two component list
# to have name of the movie included

cmpOverlayData=function(component.trackl){
    
    # keep this, although this is no need as when do.call(rbind, list) convert
    # names into the names
    # add one more column into data.frame as identifier
    # column, then collaps tracks for (i in 1:length(component.trackl)){
    # component.trackl[[i]]=lapply(component.trackl[[i]],function(x,cmp.id){
    # component=rep(cmp.id,dim(x)[1]) cmp.id.x=cbind(x,component)
    # return(cmp.id.x) },cmp.id=names(component.trackl)[i]) }
    
    
    plot.title=names(component.trackl)
    
    cat("\nProcessing",names(component.trackl))
    # cat("\nProcessing",names(component.trackl[[1]]))
    
    # combine components into cmp.lst
    cmp.lst=list()
    length(cmp.lst)=length(component.trackl[[1]])
    cmp.name=names(component.trackl[[1]])
    
    # replace file.name ".",with "_",
    # as it interference with Index indentifier "."
    cmp.name=gsub('\\.', '_', cmp.name)
    names(cmp.lst)=cmp.name
    
    # for (i in 1:length(component.trackl)){
    #     cmp.lst[[i]]=do.call(rbind.data.frame,component.trackl[[i]])
    # }
    
    for (i in 1:length(component.trackl[[1]])){
        for (j in 1:length(component.trackl[[1]][i])){
            cmp.lst[[i]]=do.call(rbind.data.frame,component.trackl[[1]][[i]])
        }
    }
    
    cmp.df=do.call(rbind.data.frame,cmp.lst)
    
    # TODO: these name manipulations is used so frequently should be put into a
    # function
    
    # split rownames
    n=cmp.df
    
    if (length(grep("txt",rownames(n)[1])) == 0){
        Index=strsplit(rownames(n),"\\.")
    }else{
        Index=strsplit(rownames(n),".txt.")
    }
    
    # trackID=fileID.frameID.duration.indexPerFile.indexPerTrackll
    Index.df=data.frame(do.call(rbind,Index))
    # do.call(rbind.data.frame,Index)
    
    colnames(Index.df)=c("component","fileID","frameID","duration",
                         "indexPerFile","indexPerTrackll","SN")
    
    track.plot.df=cbind(cmp.df,Index.df)
    
    return(track.plot.df)
    
}

## TODO:
# replace or combine readTiff to eBIage
# TODO: max.pixel=128 can be removed


##------------------------------------------------------------------------------
##

#

# these below function can be used to add image as background

# .nucleiGrob=function(nuclei.file){
#
#     # read in tiff nuclei file
#     # library(EBImage)
#
#     title=basename(nuclei.file)
#     cat("\nReading nuclei file",title,"\n")
#     nuclei=EBImage::readImage(nuclei.file)
#     # display(nuclei)
#     # interpolate =F to make pixel clearer
#     nucleiGrob= grid::rasterGrob(nuclei, interpolate=F)
#
#     #     ggplot(data=pos.point,aes(x=x,y=y),shape=22)+
#     #         annotation_custom(grob=g,xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
#     #         geom_point(color="red")
#
#     #     ggplot(data.frame(x=1:10,y=1:10),aes(x=x,y=y))+
#     #         annotation_custom(grob=g,xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
#     #         geom_point()
#
#     return(nucleiGrob)
# }

# for one file (i.e. trackl)
# .plotPicTrackOverlay=function(trackl,nucleiGrob,max.pixel=128){
#
#     # get names of the trackll (/video) to put it on each graph
#     plot.title=names(trackl)
#
#     track.overlay.data=trackOverlayData(trackl)
#
#     p=ggplot(track.overlay.data,aes(x=x,y=y,group=indexPerFile))+
#         annotation_custom(grob=nucleiGrob,xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
#
#         geom_path(color="red")+
#         scale_x_continuous(
#             name="Pixel",
#             breaks=seq(from=0, to=max.pixel,by=20),
#             limits=c(0,max.pixel))+
#         scale_y_continuous(
#             name="Pixel",
#             breaks=seq(from=0, to=max.pixel,by=20),
#             limits=c(0,max.pixel))+
#         ggtitle(plot.title)+
#
#         # this makes integer breaks
#         #        scale_x_continuous(breaks=scales::pretty_breaks(n=5))+
#         #        scale_y_continuous(breaks=scales::pretty_breaks(n=5))+
#         #labs(x="Pixel", y="Pixel")+
#
#         theme_bw()
#
#     plot(p)
#     return(p)
#
# }

# plotPicTrackOverlay=function(folder,mask=F,cores=1,
#                              max.pixel=128,
#                              nrow=2,ncol=2,width=16,height=16){
#
#     nuclei.lst=list.files(path=folder,pattern="_Nuclei.tif",full.names=T)
#
#     nucleiGrob.lst=lapply(nuclei.lst,.nucleiGrob)
#     names(nucleiGrob.lst)=sapply(nuclei.lst,basename,
#                                  simplify = TRUE, USE.NAMES = F)
#
#
#     trackll=readDiatrack(folder=folder,merge=F,mask=mask,cores=cores)
#
#     plot.lst=list()
# #     for (i in 1:length(trackll)){
# #         .plotNucTrackOverlay(trackll,nrow=1,ncol=1,width=8,height=8)
# #
# #     }
#
#     for (i in 1:length(trackll)) plot.lst[[i]]=.plotNucTrackOverlay(
#         trackl=trackll[i],nucleiGrob.lst[[i]],max.pixel=max.pixel)
#
#
#     # output
#     cat("\nOutput combined plot...")
#     cmb.plot=gridExtra::marrangeGrob(plot.lst,nrow=nrow,ncol=ncol)
#
#     fileName=paste(.timeStamp("NucTrackOverlay"),".pdf",sep="")
#
#     # fileName=paste("TrackOverlay-",.timeStamp("folder"),".pdf",sep="")
#
#     ggplot2::ggsave(filename=fileName,cmb.plot,width=width,height=height)
#
#     cat("\nDone!")
#
#
# }
