## densityMaskTracks-methods
##' @name densityMaskTracks
##' @aliases densityMaskTracks plotLines plotLinesTrackl plotPoints 
##' plotPointsTrackl
##' @title densityMaskTracks
##' @rdname densityMaskTracks-methods
##' @docType methods
##'
##' @description mask track lists and lists of track lists using kernel 
##' density clusters

##' @usage 
##' densityMaskTracks(trackll, scale = 128, removeEdge = FALSE, 
##' separate = FALSE, buildModel = FALSE)
##' plotLines(trackll, scale = 128)
##' plotLinesTrackl(track.list, scale = 128)
##' plotPoints(trackll, scale = 128)
##' plotPointsTrackl(track.list, scale = 128)
##' 
##' @param trackll An uncensored/unfiltered list of track lists
##' @param track.list Track list
##' @param scale X and Y scale (in pixels) of track video window
##' @param removeEdge Remove edge clusters with incomplete contour 
##' lines/ploygons
##' @param separate Separate by cluster
##' @param buildModel Manually configure the kernel density probability (p), 
##' while continuously building a model. T to create a model or improve an 
##' existing model. F to load in a MODEL.csv for automatic masking.
##' @return masked trackll
##' @details
##' 
##' This algorithm relies on one parameter, the kernel density probability (p),
##' to mask track lists. Following, describes a method to optimize a workflow 
##' to predict p actively.
##' 
##' When densityMaskTracks() is called with buildModel = TRUE, it will 
##' repeatedly ask the user for the kernel density probability (p) and the 
##' number of smallest clusters to elimnate. The kernel density probability 
##' (p) is a factor that determines the cluster contours density. Low p creates
##' smaller and/or fewer clusters and vice versa. Adjust p accordingly, but if
##' there are still small extra clusters made in undesired areas, raise the
##' number of smallest clusters to eliminate accordingly (note: sometimes noise
##' clusters are too small to see). Manual input will get progressively easier
##' after 3 data points as the model is continuously being improved and 
##' applied. Building the model will create a MODEL.csv in the working 
##' directory that can be used to mask automatically if buildModel = FALSE 
##' (this will look for one MODEL.csv in the working directory).
##' 
##' The separate parameter allows users to separate each track list from one
##' video into their cluster components, creating a list of track lists. 
##' Applying this separate parameter to a list of track lists from multiple
##' videos will simply append all the separated clusters together. Each track
##' list is named 'c#' as the header. The # indicating the cluster number.
##' 
##' The removeEdge parameter allwos users to automatically remove any clusters
##' that are on edges and/or have an incomplete contour line (discontinuous
##' polgon)
##' 
##' Use plotPoints and plotLines to plot lists of track lists into
##' separate scatter/line plots. Use plotPointsTrackl and plotLinesTrackl for a
##' single track list. These track lists can be plotted at any point in
##' analysis.
##' 
##' EXTRA:
##' 
##' The general method for creating a masked track list from a single 
##' track list begins by first calculating its kernel density using 
##' kernelDensity(), creating the mask using this value createMask() (while 
##' adjusting the kernel density probability [p] as needed), then generating 
##' the masked track list using applyMask. The reason why these three steps 
##' are separated in the source code is to allow for quick repeated 
##' adjustments to the kernel density probability (p), as the other two steps 
##' can take more time.
##' 
##' The value for the kernel density probability (p) is automatically 
##' calculated by default using a regression model estimating the approximate 
##' desired probability using the track list's average track length (derived 
##' from a specific data set).Thus, it is highly recommended to use
##' uncensored/unfiltered track lists.
##' 
##' If one had to apply this function to a large number of typically unvariable
##' track lists, a regression model can be made and integrated into the source
##' code.
##' 
##' 
##' when building initial masking model, MODEL.csv will be created in the
##' working directory (this can be renamed as long as it ends with MODEL.csv).
##' Each time the model is improved, the manual input will get progressively
##' easier. The output will be a masked trackll that is created exactly as set
##' during the program.

##' @examples
##' 


##' # create trackll 
##' track.folder=system.file('extdata','SWR1',package='sojourner')
##' trackll <- createTrackll(folder=track.folder, input = 3)
##'
##' # mask trackll using using default model (may not fit all data)
##' trackll.masked=densityMaskTracks(trackll)
##'
##' # plot to see the masking effect
##' plotNucTrackOverlay(folder=track.folder, trackll=trackll)
##' plotNucTrackOverlay(folder=track.folder, trackll=trackll.masked)
##' 
##' # One can also use plotLines and plotPoints for simple visualization
##' plotPoints(trackll.masked)
##' plotLines(trackll.masked)
##' 
##' # To visualize a single trackl (the trackll in the example only contains
##' # one trackl)
##' plotPointsTrackl(trackll.masked[[1]])
##' plotLinesTrackl(trackll.masked[[1]])
##' 
##' # create trackll by build model manually, useful when default model doesn't
##' # yield good masking either too strigent or too lose
##' ###trackll.masked.md <- densityMaskTracks(trackll, buildModel = TRUE)
##'
##' # model building is recommended for new dataset, currently one needs to
##' # manually evaluate the goodness of the masking, and strength /lose it by
##' # specifying p value when prompted by the command.
##'
##' # This is an example of masking model building process
##' 
##' # Masking mage6 ...
##' # No initial model read. If desired, ensure there is one file ending in
##' # MODEL.csv in working directory.
##' # 2 clusters. mage6 masked at kernel density probability = 0.5, 
##' # eliminate = 0 
##'
##' # Done (1 = YES; 0 = NO)?: 0 # assume too strigent answer no. 
##' # New kernel density probability (p): 0.6 # assume too strigent, 
##' # losen it to 0.6 
##' # Number of smallest clusters to elimnate (recommended 0, 
##' # unless last resort): 0 as recommended
##'
##' # 2 clusters. mage6 masked at kernel density probability = 0.6,
##' # eliminate = 0 
##' # Done (1 = YES; 0 = NO)?: 1 # yes. 
##'
##' # New MODEL.csv created.
##' # Masked track list for mage6 created.
##' # All tracks lists masked.
##'
##' # compare results
##' ###plotTrackOverlay(trackll) # original
##' ###plotTrackOverlay(trackll.masked) # masked using default masking model
##' ###plotTrackOverlay(trackll.masked.md) # masked using losened up masking 
##' # model 
##'
##' # now one can use this modified Model.csv to process data that was not
##' # masked well using default model by simply putting Model.csv in working
##' # directory and set buildModel=FALSE
##' ###trackll.masked2 <- densityMaskTracks(trackll, buildModel = FALSE)
##' ###plotTrackOverlay(trackll.masked.md) # masked by building new (lossen up)
##' # model
##' ###plotTrackOverlay(trackll.masked2) # masked using provided (lossen up) 
##' # model without building it

##' @export densityMaskTracks
##' @export plotPoints
##' @export plotLines
##' @export plotPointsTrackl
##' @export plotLinesTrackl


############################################################################### 

# library(dplyr) #bind_rows, MASS::kde2d library(sp) #point.in.polygon

#### kernelDensity ####

# Returns kernel density

kernelDensity = function(track.list, scale = 128) {
    
    # Merge track list into a single dataframe
    df <- mergeAllPoints(track.list)
    
    # Calculate kernel density from dataframe
    dens <- MASS::kde2d(df[[1]], df[[2]], n = 200, lims = c(c(0, scale), 
        c(0, scale)))
    
    return(dens)
}

#### createMask ####

# Returns binary mask and plots

createMask = function(track.list, scale = 128, kernel.density, p = NULL, 
    eliminate = NULL, plot = TRUE, separate = FALSE, removeEdge = FALSE) {
    # Store all merged track coordinate points into a dataframe
    df <- mergeAllPoints(track.list)
    
    if (is.null(p)) {
        p = 0.5
    }
    
    if (is.null(eliminate)) {
        eliminate = 0
    }
    
    if (p <= 0 || p >= 1) {
        cat("\np set to safe value of 0.5.\n")
        p = 0.5
    }
    
    # Calculate contours to plot
    prob <- c(p)
    dx <- diff(kernel.density$x[seq_len(2)])
    dy <- diff(kernel.density$y[seq_len(2)])
    sz <- sort(kernel.density$z)
    c1 <- cumsum(sz) * dx * dy
    levels <- sapply(prob, function(x) {
        approx(c1, sz, xout = 1 - x)$y
    })
    
    # Create the countour polygon with using coordinate points
    ls <- contourLines(kernel.density, levels = levels)
    
    # Remove edge/border clusters by removing all contour lines that are
    # not complete polygons
    if (removeEdge) {
        i = 1
        repeat {
            
            # Compare the first and last points in the polygon vector for
            # continuity
            if (tail(ls[[i]]$x, n = 1) != head(ls[[i]]$x, n = 1) || 
                tail(ls[[i]]$y, n = 1) != head(ls[[i]]$y, n = 1)) {
                ls[[i]] <- NULL
            } else {
                i = i + 1
            }
            
            # Break loop if complete
            if (is.null(ls) || i >= length(ls) + 1) {
                break
            }
        }
    }
    
    # Keep only the largest user-specified number of clusters, if given
    if (eliminate > 0) {
        num.clusters = length(ls) - eliminate
        while (length(ls) > num.clusters) {
            noise = 0
            min = Inf
            for (i in seq_along(ls)) {
                if (length(ls[[i]][[2]]) < min) {
                  noise = i
                  min = length(ls[[i]][[2]])
                }
            }
            ls[[noise]] <- NULL
        }
    }
    
    # Use coordinate coordinate polygon to create the cluster shape
    cluster <- list()
    for (i in seq_along(ls)) {
        cluster[[i]] <- point.in.polygon(df[[1]], df[[2]], ls[[i]]$x, ls[[i]]$y)
    }
    
    # Create binary mask of track coordinates
    df$region <- factor(Reduce("+", cluster))
    
    # Plot with mask and contour
    if (plot) {
        title = paste(getTrackFileName(track.list), 
            "Mask with Kernel Density Probability (p) of", 
            round(p, digits = 3), sep = " ")
        plot(df[[2]] ~ df[[1]], col = df$region, data = df, xlim = c(0, 
            scale), ylim = c(0, scale), xlab = "x (Pixels)", 
            ylab = "y (Pixels)", main = title, cex = 0.1)
        contour(kernel.density, levels = levels, labels = prob, add = TRUE)
    }
    cat("\n", length(ls), "clusters.", getTrackFileName(track.list), 
        "masked at kernel density probability =", 
        round(p, digits = 3), ", eliminate =", eliminate, "\n")
    
    if (!separate) {
        return(df$region)
    } else {
        return(cluster)
    }
}

#### applymask ####

# Creates masked track list

applyMask = function(track.list, mask) {
    
    # Instantiate a masked track list with indexing variables
    masked.track.list = list()
    masked.track.list.names = list()
    index.mask = 1
    index = 1
    
    # Loop through all tracks
    for (i in seq_along(track.list)) {
        mask.bool = TRUE
        
        # Remove any tracks outside mask
        for (j in seq_len(nrow(track.list[[i]]))) {
            if (mask[[index]] == 1) {
                mask.bool = FALSE
            }
            index = index + 1
        }
        if (!mask.bool) {
            masked.track.list[index.mask] <- track.list[i]
            index.mask = index.mask + 1
            masked.track.list.names[1 + length(masked.track.list.names)] = 
                names(track.list[i])
        }
    }
    names(masked.track.list) <- masked.track.list.names
    # Return masked track list
    return(masked.track.list)
}

#### mergeAllPoints ####

mergeAllPoints = function(track.list) {
    if (length(track.list[[1]]) == 3) {
        df <- bind_rows(track.list, .id = "Trajectory")[, c("x", "y", "z")]
    } else {
        df <- bind_rows(track.list, .id = "Trajectory")[, c("x", "y", "z", 
            "Frame")]
    }
    
    return(df)
}

#### .densityMaskTracks ###

.densityMaskTracks = function(track.list, scale = 128, removeEdge = FALSE, 
    separate = FALSE, buildModel = FALSE) {
    
    # Initial confirmation
    track.name = getTrackFileName(track.list)
    cat("\n Masking", track.name, "...\n")
    
    # Calculate kernel density
    kd <- kernelDensity(track.list, scale = scale)
    
    
    # Set model variables
    points = nrow(mergeAllPoints(track.list))
    tracks = length(track.list)
    avg = points/tracks
    
    # Find initial model
    model.fit = NULL
    model = NULL
    model.file = list.files(path = getwd(), pattern = "MODEL.csv", full.names = TRUE, 
        ignore.case = TRUE)
    if (length(model.file) != 1) {
        cat(paste("\nNo initial model read. If desired,", 
            "ensure there is one file ending in MODEL.csv", 
            "in working directory.\n"))
        p = 0.5
    } else {
        model = read.csv(model.file)
        if (nrow(model) < 3) {
            cat(paste("\n", basename(model.file), " read.\n", sep = ""))
            p = 0.5
            cat(paste("\nNote: Model has less than 3 points.", 
                "Need more data. p set to safe default 0.5.\n"))
        } else {
            model.fit <- summary(
                fit <- lm(model$Probability ~ model$Average.Track.Length))
            cat(paste("\n", basename(model.file), " read. R-squared = ", 
                round(model.fit$r.squared, digits = 3), "\n", sep = ""))
            cat("\nPredicting p...\n")
            p = model.fit$coefficients[[1]] + model.fit$coefficients[[2]] * 
                avg
            if (!(p > 0 || p < 1)) {
                cat("\nWarning: Model inaccuracy! p set to safe default 0.5\n")
                p = 0.5
            }
        }
    }
    
    # Option between automatic and manual
    if (!buildModel) {
        
        # Automatically initial mask
        mask <- createMask(track.list, scale = scale, kd, p = p, 
            separate = separate, removeEdge = removeEdge)
        
    } else {
        new.model <- NULL
        
        # Instantiate to default parameters
        eliminate = 0
        done = FALSE
        
        # reate mask using default mask
        mask <- createMask(track.list, scale = scale, kd, p = p, 
            separate = separate, removeEdge = removeEdge)
        
        # Repeatedly ask if satisfied with mask Infinite loop issue
        while (!done) {
            
            # Done prompt
            cat("\n")
            done = as.integer(readline(prompt = "Done (1 = YES; 0 = NO)?: "))
            
            # If invalid input / not a number, set to default 0
            if (is.na(done)) {
                cat("\nInvalid input, set to default 0 = NO.\n")
                done = 0
            }
            
            # Convert integer to logical
            done = as.logical(done)
            
            # If done, add to new model data frame and break
            if (done) {
                new.model <- data.frame(track.name, points, tracks, 
                    points/tracks, p)
                colnames(new.model) <- c("File Name", "Points", "Tracks", 
                  "Average Track Length", "Probability")
                break
            }
            
            # Request new p
            p = as.numeric(readline(
                prompt = "New kernel density probability (p): "))
            
            # If invalid input / not a number, set to default
            if (is.na(p)) {
                cat("\nInvalid input, set to default.\n")
                p = NULL
            }
            
            # Request to eliminate
            eliminate = as.integer(readline(prompt = paste(
                "Number of smallest clusters to elimnate", 
                "(recommended 0, unless last resort): ")))
            
            # If invalid input / not a number, set to default
            if (is.na(eliminate)) {
                cat("\nIncorrect input, set to default = 0.\n")
                eliminate = 0
            }
            
            # Create mask using set parameter
            mask <- createMask(track.list, scale = scale, kd, p = p, 
                eliminate = eliminate, separate = separate, 
                removeEdge = removeEdge)
        }
        
        # Create new MODEL.csv or append new model accordingly
        if (is.null(model)) {
            cat("\nNew MODEL.csv created.\n")
            write.table(new.model, file = "MODEL.csv", sep = ",", 
                        row.names = FALSE)
        } else {
            write.table(new.model, file = basename(model.file), sep = ",", 
                append = TRUE, col.names = FALSE, row.names = FALSE)
            cat(paste("\nData point added to ", basename(model.file), ".\n", 
                sep = ""))
        }
    }
    
    # Apply mask depending on desire to separate clusters
    if (typeof(mask) == "integer") {
        
        # Apply full mask to track list and return
        masked.track.list <- applyMask(track.list, mask)
        cat("\n Masked track list for", getTrackFileName(track.list), 
            "created.\n")
        return(masked.track.list)
    } else {
        
        # Create new list of track lists for separated clusters
        masked.trackll.cells <- list()
        masked.trackll.cells.names <- list()
        
        # Loop through separated masks and apply to track list
        for (i in seq_along(mask)) {
            masked.trackll.cells[[length(masked.trackll.cells) + 1]] <- 
                applyMask(track.list, mask[[i]])
            masked.trackll.cells.names[[length(masked.trackll.cells.names) + 
                1]] <- i
        }
        
        # Name separated list of track lists apppropriately and return
        names(masked.trackll.cells) <- paste("c", masked.trackll.cells.names, 
            sep = "")
        cat("\n Masked list of track lists (separated by cluster) for", 
            getTrackFileName(track.list), "created.\n")
        return(masked.trackll.cells)
    }
    
}

#### densityMaskTracks ####

densityMaskTracks = function(trackll, scale = 128, removeEdge = FALSE, 
    separate = FALSE, buildModel = FALSE) {
    
    # Instantiate empty list
    masked.trackll <- list()
    
    # Option to separate
    if (separate) {
        
        # Apply separation to each track list and append the resulting list of
        # track lists to each other
        for (i in seq_along(trackll)) {
            tracks <- .densityMaskTracks(trackll[[i]], scale = scale, 
                removeEdge = removeEdge, separate = separate, 
                buildModel = buildModel)
            masked.trackll <- append(masked.trackll, tracks)
        }
        names(masked.trackll) <- paste(names(trackll), names(masked.trackll), 
            sep = "_")
    } else {
        
        # Apply algorithm to each track list
        for (i in seq_along(trackll)) {
            masked.trackll[[i]] <- .densityMaskTracks(trackll[[i]], 
                scale = scale, removeEdge = removeEdge, separate = separate, 
                buildModel = buildModel)
        }
        names(masked.trackll) <- names(trackll)
    }
    
    # Confirmation text and return
    cat("\nAll tracks lists masked.\n")
    return(masked.trackll)
}

#### plotPoints ####

plotPoints = function(trackll, scale = 128) {
    for (i in seq_along(trackll)) {
        plotPointsTrackl(trackll[[i]], scale = scale)
        title(sub = names(trackll)[[i]])
    }
}

plotPointsTrackl = function(track.list, scale = 128) {
    df <- mergeAllPoints(track.list)
    plot(df[[1]], df[[2]], xlim = c(0, scale), ylim = c(0, scale), 
        xlab = "x (Pixels)", ylab = "y (Pixels)", 
        main = paste("Tracks Plot for ", getTrackFileName(track.list), 
        sep = ""), cex = 0.1)
}

#### plotLines ####

plotLines = function(trackll, scale = 128) {
    for (i in seq_along(trackll)) {
        plotLinesTrackl(trackll[[i]], scale = scale)
        title(sub = names(trackll)[[i]])
    }
}

plotLinesTrackl = function(track.list, scale = 128) {
    plot(track.list[[1]][[1]], track.list[[1]][[2]], type = "l", xlim = c(0, 
        scale), ylim = c(0, scale), xlab = "x (Pixels)", ylab = "y (Pixels)", 
        main = paste("Tracks Plot for ", getTrackFileName(track.list), 
            sep = ""))
    for (i in 2:length(track.list)) {
        lines(track.list[[i]][[1]], track.list[[i]][[2]])
    }
}

