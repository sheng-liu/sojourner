
## dwellTime-methods
##' @name dwellTime
##' @aliases dwellTime
##' @title dwellTime
##' @rdname dwellTime-methods
##' @docType methods
##' @description Caclulate dwell time (/residence time) for trajecotries.

##' @usage dwellTime(trackll,t.interval=10,x.scale=c(min=0,max=250),plot=TRUE,
##' output=FALSE)
##' @param trackll Track list output from readDiatrack().
##' @param t.interval t.interval time, default = 10ms.
##' @param x.scale x-scale min and max range.
##' @param plot An Logical indicate if plot should be generated. If 
##' plot = TRUE, the plot data will also be output.
##' @param output An Logical indicate if output should be generated. 1) dwell 
##' time of tracks in the track list output to csv file. Each item in the list 
##' will have an individual csv file. 2) Plot PDF and plot data will be saved.


##' @return
##' \itemize{
##' \item{dwell time list} A list of dwell time for all trajectories, 
##' separated by file names of the trajectory file in Diatrack file folder. 
##' If combined dwell time is intended, use readDiatrack(folder, merge=TRUE) 
##' to generate a  single length list, then apply this function.
##' \item{PDF} dwell time frequency plot in PDF format, when plot = TRUE.
##' \item{csv} dwell time output in csv format, when output = TRUE.
##' }

##' @examples
##' folder=system.file('extdata','SWR1',package='sojourner')
##' trackll=createTrackll(folder=folder, input=3)
##' dwellTime(trackll,plot=TRUE)

##' @import reshape2
##' @export dwellTime
############################################################################### 


##-----------------------------------------------------------------------------
## .dwellTime a function to calculate dwell time from a list of data.frame
## track (trackl) and returns a vector of dwell time.

## nomenclature track data.frame with x,y,z coordinates trackl list of
## data.frames with x,y,z coordinates, read from one track file trackll list
## of list of data.frames with x,y,z coordinates, read from bmultiple track
## file

.dwellTime = function(trackl, t.interval = 10) {
    vapply(trackl, function(x) {
        dim(x)[1] * t.interval
    }, FUN.VALUE=double(1))
}


dwellTime = function(trackll, t.interval = 10, x.scale = c(min = 0, max = 250), 
    plot = TRUE, output = FALSE) {
    
    ## compute dwell time
    
    #dwell.times = sapply(trackll, function(x) {
    #    .dwellTime(x, t.interval)
    #})
    
    dwell.time = lapply(trackll, function(x) {
        .dwellTime(x, t.interval)
    })
    # lapply() replacement with sapply() requires special case
    if (length(trackll)==1){
        dwell.time = matrix(unlist(dwell.time), 
            dimnames = list(names(unlist(dwell.time))))
        
        colnames(dwell.time) <- names(trackll)
        rownames(dwell.time) = lapply(rownames(dwell.time),function(x){
            gsub(paste(names(trackll),".",sep=""), "",x)})
    }
    
    file.name = names(trackll)
    
    ## reshape data for plot
    dwell.time.mlt = reshape2::melt(dwell.time)
    
    if (length(trackll) == 1) {
        
        colnames(dwell.time.mlt) = c("index", "variable", "value")
    } else {
        colnames(dwell.time.mlt) = c("value", "variable")
        
    }
    
    histo.plot = ggplot(dwell.time.mlt, aes_string(x = "value", 
        group = "variable", fill = "variable")) + 
        geom_histogram(binwidth = t.interval, position = "dodge", 
        colour = "white") + ## change from white to red?
    xlim(x.scale["min"], x.scale["max"]) + theme_bw() + 
        theme(legend.title = element_blank()) + 
        labs(x = "Lifetime of trajectories (ms)", y = "Number of trajecotries")
    
    density.plot = ggplot(dwell.time.mlt, aes_string(x = "value", 
        group = "variable", col = "variable", fill = "variable")) + 
        geom_density(alpha = 0.2) + xlim(x.scale["min"], x.scale["max"]) + 
        theme_bw() + theme(legend.title = element_blank()) + 
        labs(x = "Lifetime of trajectories (ms)", 
        y = "Frequency of trajectories")
    
    if (plot == TRUE) 
        multiplot(histo.plot, density.plot, cols = 1)
    
    ## output
    if (output == TRUE) {
        
        # output csv for (i in seq_along(trackll)){ fileName=paste('Dwell
        # Time-',.timeStamp(file.name[i]), '.csv',sep='')
        # write.csv(file=fileName,dwell.time[[i]]) }
        
        # output plot
        tStamp.plotName = paste(.timeStamp(file.name[1]), "___", sep = "")
        plotName = paste("Dwell Time Plot-", tStamp.plotName, ".pdf", sep = "")
        ggsave(filename = plotName, plot = histo.plot, width = 8, height = 4)
        
        # output plot data
        plotData = ggplot_build(histo.plot)$data
        plotFile = paste("Dwell Time Plot-", tStamp.plotName, ".csv", sep = "")
        write.csv(file = plotFile, plotData)
    }
    return(invisible(dwell.time))
}


##-----------------------------------------------------------------------------
## 

# freqpoly=ggplot(dwell.time.mlt,aes(x=value,color=variable)) +
# geom_freqpoly(binwidth=t.interval)+labs(x='Dwell time (ms)',y='Count')+
# theme_bw()+ theme(legend.title=element_blank())+xlim(0,200)

# histodensity=ggplot(dwell.time.mlt,aes(x=value,color=variable,
# fill=variable))+ geom_histogram(binwidth=t.interval,position='dodge')+
# geom_density(aes(y=10*..count..),alpha=0.2)+ theme_bw()+
# theme(legend.title=element_blank())+ xlim(0,200)



