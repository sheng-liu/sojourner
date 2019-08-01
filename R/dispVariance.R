## dispVariance-methods
##' @name dispVariance
##' @aliases dispVariance
##' @title dispVariance
##' @rdname dispVariance-methods
##' @docType methods
##' 
##' @description calculate square displacements for all tracks in a trackll 
##' datatype, and return the variances for the dispacements of each 
##' trajectories.
##' 
##' @usage
##' dispVariance(trackll, min=7, plot=FALSE, limits=c(), log=FALSE, 
##' output=FALSE)
##' @param trackll a list of track lists.
##' @param min minimum points on trajectory, should be at least 3 to work.
##' @param plot default: False, if true, show density plot for variances.
##' @param limits vector of size2, variance cut-off range that one wants to 
##' plot.
##' This will not affect the returned result.
##' @param log default: False, if true, apply log10 to variance value for new 
##' spread. like limits, this will only affect the plot, not the returned 
##' value.
##' @param output if True, generate a csv output for each tracklist files that 
##' are in trackll.
##' @return \itemize{
##' \item{Variances} calculated variacne for all trakcs in trackll}
##' 
##' @examples
##' folder=system.file('extdata','SWR1',package='sojourner')
##' trackll=createTrackll(folder=folder, input=3)
##'
##' # run dispVariance with default minimum tracklength (min=7)
##' dispVars = dispVariance(trackll)
##' 
##' # run dispVariance by setting the minimum tracklength to 3
##' dispVars = dispVariance(trackll, min=3)
##' 
##' # display plot only within certain range
##' dispVariance(trackll, plot=TRUE, limits = c(0,0.002))
##' 
##' # display plot with log-scale applied
##' dispVars = dispVariance(trackll, min=3, plot=TRUE, log=TRUE)
##' 
##' # display plot. Could get csv files if output = TRUE
##' dispVars = dispVariance(trackll, min=3, plot=TRUE, output=FALSE)
##' 
##' @details 
##' dispVariance applies the squareDisp function to each dataframe
##' containing trajectories. the tracks somehow had to be converted
##' into dataframes although they were expected to be in dataframes
##' in the first place.
##' 
##' Since the tracks of shorter length are filtered out in the process,
##' there is no guqrante that the length of tracklists equal that of 
##' the input.
##' 
##' The tracks should have length of at least 3, in order to have a
##' valid displacement variance. If min argument is less than 3, the
##' function will not be executed.
##' 
##' Generally, when plotting, you would want to use only one of limits
##' or log. Although you may use both, using only one of the two would
##' do the job.
##' 
##' @export dispVariance

############################################################################### 

##-----------------------------------------------------------------------------
## dispVariance_track

## calculate variance for displacements single track

##' @export dispVariance_track
dispVariance_track = function(track) {
    withDisp = squareDisp(track)  # get displacement
    sqDisp = withDisp[[1]]$square.disp  #[[1]] is necessary b/c 
    # squaredisplacement give a list in result
    clean = sqrt(sqDisp[!is.na(sqDisp)])  #get rid of NA and do the sqrt since 
    # it is currently dx^2 + dy^2
    var(clean)  # return variance value
}

##-----------------------------------------------------------------------------
## dispVariance.trackl

## calculate variance for displacements single tracklist
dispVariance.trackl = function(trackl) {
    
    lapply(trackl, dispVariance_track)
    # track.df = lapply(trackl, data.frame) track.withdisp =
    # sapply(track.df, squareDisp) sqdisp =
    # lapply(track.withdisp,function(x) {x[7]})#column with square-disp
    # sqdisp.clean = lapply(sqdisp, function(x) {x[!is.na(x)]}) disp.var =
    # lapply(sqdisp.clean, var) disp.var
}

##-----------------------------------------------------------------------------
## dispVariance

## calculate displacement variance for all tracks in given trackll also
## filter based on tracklength(must be greater than 2) plotting feature
## can be enabled when plot=TRUE is set. limits on the variance range
## can be set by providing a vector of length 2.

##' @export dispVariance
dispVariance = function(trackll, min = 7, plot = FALSE, limits = c(), 
    log = FALSE, output = FALSE) {
    if (min < 3) {
        stop("min value should be at least 3")
    }
    
    # filter by length
    filtered = filterTrack(trackll, filter = c(min = min, max = Inf))
    result = lapply(filtered, dispVariance.trackl)
    
    if (plot == TRUE) {
        melted = reshape2::melt(result)  #convert to dataframe
        # rename columns to more reasonable ones
        names(melted) = c("variance", "track.name", "trackList")
        if (log) {
            melted$variance = log10(melted$variance)
        }
        # different tracklists will show up as different colors with some
        # transparency
        plt = ggplot2::ggplot(melted, ggplot2::aes_string(x = "variance", 
            color = "trackList")) + ggplot2::geom_line(alpha = 0.5, 
            position = "identity", stat = "density")
        if (length(limits) != 2) {
            plot(plt)
        } else {
            plot(plt + ggplot2::xlim(limits))  #apply range limits
        }
    }
    
    if (output == TRUE) {
        for (i in seq_along(result)) {
            track.df = reshape2::melt(result[[i]])
            names(track.df) = c("dispVariance", "track.name")
            track.df = track.df[, c(2, 1)]
            fileName = paste("DispVariance Individual-", 
                .timeStamp(names(result)[i]), ".csv", sep = "")
            cat("\nOutput dispVariance for individual trajectories.\n")
            write.csv(file = fileName, track.df)
        }
    }
    return(result)
}
