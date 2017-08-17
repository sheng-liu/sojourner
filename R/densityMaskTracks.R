#### densityMaskTracks.R
#### Wu Lab, Johns Hopkins University
#### Author: Sun Jay Yoo
#### Date: July 21, 2017

## densityMaskTracks-methods
##
##
###############################################################################
##' @name densityMaskTracks
##' @aliases densityMaskTracks
##' @title densityMaskTracks
##' @rdname densityMaskTracks-methods
##' @docType methods
##'
##' @description mask track lists and lists of track lists using kernel density clusters

##' @usage 
##' .densityMaskTracks(trackll, scale = 128, removeEdge = F, automatic = F, separate = F)
##' 
##' densityMaskTracks(track.list, scale = 128, removeEdge = F, automatic = F, separate = F)
##' 
##' plotPoints(track.list, scale = 128)
##' 
##' plotLines(track.list, scale = 128)
##' 

##' @param trackll An uncensored/unfiltered list of track lists
##' @param scale X and Y scale (in pixels) of track video window
##' @param removeEdge Remove edge clusters with incomplete contour lines/ploygons
##' @param automatic Find p automatically using a model (not recommended)
##' @param separate Separate by cluster
##' @param track.list A single uncensored/filtered track list.

##' @details
##' 
##' When densityMaskTracks() is called by default with automatic = F, it will repeatedly ask the user for the kernel density probability (p)
##' and the number of smallest clusters to elimnate. The kernel density probability (p) is a factor that determines how dense the cluster contours are.
##' Low p creates smaller and/or fewer clusters and vice versa. Adjust p accordingly, but if there are still small extra clusters made in undesired
##' areas, raise the number of smallest clusters to eliminate accordingly (note: sometimes noise clusters are too small to see). 
##' Use .densityMaskTracksl() to apply this to only one track list.
##' 
##' The separate parameter allows users to separate each track list from one video into their cluster components, creating a list of track lists.
##' Applying this separate parameter to a list of track lists from multiple videos will simply append all the separated clusters together.
##' Each track list is named "c#" as the header. The # indicating the cluster number.
##' 
##' The removeEdge parameter allwos users to automatically remove any clusters that are on edges and/or have an incomplete contour line (discontinuous polgon)
##' 
##' Use plotTrackPoints and plotTrackLines to plot lists of track lists into separate scatter/line plots. 
##' Use .plotTrackPoints and .plotTrackLines for a single track list. These track lists can be plotted at any point in analysis.
##' 
##' EXTRA:
##' 
##' The general method for creating a masked track list from a single track list begins by 
##' first calculating its kernel density using kernelDensity(), 
##' creating the mask using this value createMask() (while adjusting the kernel density probability [p] as needed), 
##' then generating the masked track list using applyMask. The reason why these three steps are separated in the source code is
##' to allow for quick repeated adjustments to the kernel density probability (p), as the other two steps can take more time.
##' 
##' The value for the kernel density probability (p) is automatically calculated by default using 
##' a regression model estimating the approximate desired probability using the track list's average track length 
##' (derived from a specific data set).Thus, it is highly recommended to use uncensored/unfiltered track lists.
##' 
##' If one had to apply this function to a large number of typically unvariable track lists, a regression model can be made and integrated into the source code.
##' 

##' @examples
##' 
##' #Default call for masking a list of track lists with separation by cluster.
##' masked.trackll <- densityMaskTracks(trackll, removeEdge = T)
##' 
##' #Default call for masking a track list
##' masked.trackl <- .densityMaskTracks(trackl)


##' @export .densityMaskTracks
##' @export densityMaskTracks
##' @export plotPoints
##' @export plotLines
##' @export .plotPoints
##' @export .plotLines

##' @importFrom dplyr bind_rows
##' @importFrom MASS kde2d
##' @importFrom sp point.in.polygon

###############################################################################

#library(dplyr) #bind_rows, MASS::kde2d
#library(sp) #point.in.polygon

#### kernelDensity ####

#Returns kernel density

kernelDensity = function (track.list, scale = 128){
    
    #Merge track list into a single dataframe
    df <- mergeAllPoints(track.list)
    
    #Calculate kernel density from dataframe
    dens <- MASS::kde2d(df[[1]], df[[2]], n=200, lims=c(c(0, scale), c(0, scale)));
    
    return (dens);
}

#### createMask ####

#Returns binary mask and plots

createMask = function (track.list, scale = 128, kernel.density, p = NULL, eliminate = NULL, plot = T, separate = F, removeEdge = F){
    #Store all merged track coordinate points into a dataframe
    df <- mergeAllPoints(track.list)
    
    if (is.null(p)){
        p = 0.5
    }
    
    if (is.null(eliminate)){
        eliminate = 0
    }
    
    if (p <= 0 || p >= 1){
        cat("\np set to safe value of 0.5.\n")
        p = 0.5
    }
    
    # Calculate contours to plot
    prob <- c(p)
    dx <- diff(kernel.density$x[1:2])
    dy <- diff(kernel.density$y[1:2])
    sz <- sort(kernel.density$z)
    c1 <- cumsum(sz) * dx * dy 
    levels <- sapply(prob, function(x){ 
        approx(c1, sz, xout = 1 - x)$y
    })
    
    #Create the countour polygon with using coordinate points
    ls <- contourLines(kernel.density, level=levels)
    
    #Remove edge/border clusters by removing all contour lines that are not complete polygons
    if (removeEdge){
        i = 1
        repeat {
            
            #Compare the first and last points in the polygon vector for continuity
            if (tail(ls[[i]]$x, n = 1) != head(ls[[i]]$x, n = 1) || tail(ls[[i]]$y, n = 1) != head(ls[[i]]$y, n = 1)){
                ls[[i]] <- NULL
            } else {
                i = i + 1
            }
            
            #Break loop if complete
            if (is.null(ls) || i >= length(ls) + 1){
                break
            }
        }
    }
    
    #Keep only the largest user-specified number of clusters, if given
    if (eliminate > 0){
        num.clusters = length(ls) - eliminate;
        while (length(ls) >  num.clusters){
            noise = 0;
            min = Inf;
            for (i in 1:length(ls)){
                if(length(ls[[i]][[2]]) < min){
                    noise = i
                    min = length(ls[[i]][[2]])
                }
            }
            ls[[noise]] <- NULL
        }
    }
    
    #Use coordinate coordinate polygon to create the cluster shape
    cluster <- list()
    for (i in 1:length(ls)){
        cluster[[i]] <- point.in.polygon(df[[1]], df[[2]], ls[[i]]$x, ls[[i]]$y)
    }
    
    #Create binary mask of track coordinates
    df$region <- factor(Reduce("+", cluster))
    
    #Plot with mask and contour
    if(plot){
        title = paste(getTrackFileName(track.list),"Mask with Kernel Density Probability (p) of", round(p, digits = 3), sep = " ");
        plot(df[[2]] ~ df[[1]], col=region, data=df, xlim = c(0, scale), ylim = c(0, scale), xlab = "x (Pixels)", ylab = "y (Pixels)", main = title, cex = .1)
        contour(kernel.density, levels=levels, labels=prob, add=T)
    }
    cat("\n", length(ls), "clusters.", getTrackFileName(track.list), "masked at kernel density probability =", round(p, digits = 3), ", eliminate =", eliminate, "\n")
    
    if (!separate){
        return(df$region)
    } else {
        return(cluster)
    }
}

#### applymask ####

#Creates masked track list

applyMask = function(track.list, mask){
    
    #Instantiate a masked track list with indexing variables
    masked.track.list = list();
    masked.track.list.names = list();
    index.mask = 1;
    index = 1;
    
    #Loop through all tracks
    for(i in 1:length(track.list)){
        mask.bool = TRUE;
        
        #Remove any tracks outside mask
        for (j in 1:nrow(track.list[[i]])){
            if (mask[[index]] == 1){
                mask.bool = FALSE;
            }
            index = index + 1;
        }
        if (!mask.bool){
            masked.track.list[index.mask] <- track.list[i];
            index.mask = index.mask + 1;
            masked.track.list.names[1 + length(masked.track.list.names)] = names(track.list[i]);
        }		
    }
    names(masked.track.list) <- masked.track.list.names;
    #Return masked track list
    return (masked.track.list);
}

#### mergeAllPoints ####

mergeAllPoints = function(track.list){
    if (length(track.list[[1]]) == 3){
        df <- bind_rows(track.list, .id = "Trajectory")[, c("x", "y", "z")]
    } else {
        df <- bind_rows(track.list, .id = "Trajectory")[, c("x", "y", "z", "Frame")]
    }
    
    return (df);
}

#### .densityMaskTracks ###

.densityMaskTracks = function (track.list, scale = 128, removeEdge = F, separate = F, buildModel = F){
    
    track.name = getTrackFileName(track.list)
    model = NULL;
    model.file = list.files(path=getwd(),pattern="MODEL.csv",full.names=T)
    if (length(model.file) != 1 || basename(model.file) != "MODEL.csv"){
        cat("No model read. If desired, ensure there is one MODEL.csv.")
    } else {
        cat("MODEL.csv read.")
        model = read.csv(model.file)
    }

    cat("\n Masking", track.name, "...\n")
    
    #Calculate kernel density
    kd <- kernelDensity(track.list, scale = scale);
    
    
    #Set model variables
    points = nrow(mergeAllPoints(track.list))
    tracks = length(track.list)
    avg = points/tracks
    
    #Predict p using model
    if (is.null(model)){
        p = 0.5;
        #p = -0.1207484 + 0.3468734*avg
    } else {
        if (nrow(model) < 2){
            p = 0.5
            cat("Note: Model has less than 2 points.")
        } else {
            model.fit <- summary(fit <- lm(model$Probability ~ model$Average.Track.Length))
            p = model.fit$coefficients[[1]] + model.fit$coefficients[[2]] * avg
        }
    }
    
    #Option between automatic and manual
    if (!buildModel){
        
        #Automatically create mask using default model
        mask <- createMask(track.list, scale = scale, kd, p = p, separate = separate, removeEdge = removeEdge);
        
    } else {
        new.model <- NULL;

        #Instantiate to default parameters
        eliminate = 0;
        done = FALSE;
        
        #reate mask using default mask
        mask <- createMask(track.list, scale = scale, kd, p = p, separate = separate, removeEdge = removeEdge);
        
        #Repeatedly ask if satisfied with mask
        while (!done){
            
            #Done prompt
            cat("\n")
            done = as.integer(readline(prompt="Done (1 = YES; 0 = NO)?: "))
            
            #If invalid input / not a number, set to default 0
            if (is.na(done)){
                cat("\nInvalid input, set to default 0 = NO.\n")
                done = 0;
            }
            
            #Convert integer to logical
            done = as.logical(done)
            
            #Break if done
            if (done){
                new.model <- data.frame(track.name, points, tracks, points/tracks, p);
                colnames(new.model) <- c("File Name", "Points", "Tracks", "Average Track Length", "Probability");
                break;
            }
            
            #Request new p
            p = as.numeric(readline(prompt="New kernel density probability (p): "))
            
            #If invalid input / not a number, set to default
            if (is.na(p)){
                cat("\nInvalid input, set to default.\n")
                p = NULL;
            }
            
            #Request to eliminate
            eliminate = as.integer(readline(prompt="Number of smallest clusters to elimnate (recommended 0, unless last resort): "))
            
            #If invalid input / not a number, set to default
            if (is.na(eliminate)){
                cat("\nIncorrect input, set to default = 0.\n")
                eliminate = 0;
            }
            
            #Create mask using set parameter
            mask <- createMask(track.list, scale = scale, kd, p = p, eliminate = eliminate, separate = separate, removeEdge = removeEdge);
        }
        if (is.null(model)){
            cat("New MODEL.csv created.")
            write.table(new.model, file = "MODEL.csv", sep = ",");
         } else {
            cat("Data point added to MODEL.csv.")
            write.table(new.model, file = "MODEL.csv", sep = ",", append = T, col.names = F);
        }
    }
    
    #Apply mask depending on desire to separate clusters
    if (typeof(mask) == "integer"){
        
        #Apply full mask to track list and return
        masked.track.list <- applyMask(track.list, mask);
        cat("\n Masked track list for", getTrackFileName(track.list), "created.\n")
        return(masked.track.list)
    } else {
        
        #Create new list of track lists for separated clusters
        masked.trackll.cells <- list()
        masked.trackll.cells.names <- list()
        
        #Loop through separated masks and apply to track list
        for (i in 1:length(mask)){
            masked.trackll.cells[[length(masked.trackll.cells)+1]] <- applyMask(track.list, mask[[i]]);
            masked.trackll.cells.names[[length(masked.trackll.cells.names)+1]] <- i;
        }
        
        #Name separated list of track lists apppropriately and return
        names(masked.trackll.cells) <- paste("c", masked.trackll.cells.names, sep ="");
        cat("\n Masked list of track lists (separated by cluster) for", getTrackFileName(track.list), "created.\n")
        return(masked.trackll.cells)
    }
    
}

#### densityMaskTracks ####

densityMaskTracks = function (trackll, scale = 128, removeEdge = F, separate = F, buildModel = F){
    
    #Instantiate empty list
    masked.trackll <- list()
    
    #Option to separate
    if (separate){
        
        #Apply separation to each track list and append the resulting list of track lists to each other
        for (i in 1:length(trackll)){
            tracks <- .densityMaskTracks(trackll[[i]], scale = scale, removeEdge = removeEdge, separate = separate, buildModel = buildModel)
            masked.trackll <- append(masked.trackll, tracks)
        }
        names(masked.trackll) <- paste(names(trackll), names(masked.trackll), sep = "_")
    } else {
        
        #Apply algorithm to each track list
        for (i in 1:length(trackll)){
            masked.trackll[[i]] <- .densityMaskTracks(trackll[[i]], scale = scale, removeEdge = removeEdge, separate = separate, buildModel = buildModel)
        }
        names(masked.trackll) <- names(trackll)
    }
    
    #Confirmation text and return
    cat("\nAll tracks lists masked.\n")
    return(masked.trackll)
}

#### plotPoints ####

plotPoints = function(trackll, scale = 128){
    for (i in 1:length(trackll)){
        .plotPoints(trackll[[i]], scale = scale)
        title(sub= names(trackll)[[i]])
    }
}

.plotPoints = function(track.list, scale = 128){
    df <- mergeAllPoints(track.list)
    plot(df[[1]], df[[2]], xlim = c(0, scale), ylim = c(0, scale), xlab = "x (Pixels)", ylab = "y (Pixels)", main = paste("Tracks Plot for ", getTrackFileName(track.list), sep = ""), cex = .1);
}

#### plotLines ####

plotLines = function(trackll, scale = 128){
    for (i in 1:length(trackll)){
        .plotLines(trackll[[i]], scale = scale)
        title(sub= names(trackll)[[i]])
    }
}

.plotLines = function(track.list, scale = 128){
    plot(track.list[[1]][[1]], track.list[[1]][[2]], type = "l", xlim = c(0, scale), ylim = c(0, scale), xlab = "x (Pixels)", ylab = "y (Pixels)", main = paste("Tracks Plot for ", getTrackFileName(track.list), sep = ""));
    for(i in 2:length(track.list)){
        lines(track.list[[i]][[1]], track.list[[i]][[2]])
    }
}

