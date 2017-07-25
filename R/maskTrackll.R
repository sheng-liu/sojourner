#### maskTrackll.R
#### Wu Lab, Johns Hopkins University
#### Author: Sun Jay Yoo
#### Date: July 21, 2017

## maskTrackll-methods
##
##
###############################################################################
##' @name maskTrackll
##' @aliases maskTrackll
##' @title maskTrackll
##' @rdname maskTrackll-methods
##' @docType methods
##'
##' @description mask track lists and lists of track lists using kernel density clusters

##' @usage 
##' maskTrackll(trackll, automatic = F)
##' 
##' maskTrackl(track.list, automatic = F)
##' 
##' plotPoints(track.list)
##' 
##' plotLines(track.list)
##' 

##' @param trackll An uncensored/unfiltered list of track lists.
##' @param automatic Find p automatically using a model (not recommended)
##' @param track.list A single uncensored/filtered track list.

##' @details
##' 
##' When maskTrackll() is called by default with automatic = F, it will repeatedly ask the user for the kernel density probability (p)
##' and the number of smallest clusters to elimnate. The kernel density probability (p) is a factor that determines how dense the cluster contours are.
##' Low p creates smaller and/or fewer clusters and vice versa. Adjust p accordingly, but if there are still small extra clusters made in undesired
##' areas, raise the number of smallest clusters to eliminate accordingly (note: sometimes noise clusters are too small to see). 
##' Use maskTrackl() to apply this to only one track list.
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
##' #Default call for masking a list of track lists.
##' masked.trackll <- maskTrackl(trackll)
##' 
##' #Default call for masking a track list
##' masked.trackl <- maskTrackl(trackl)


##' @export maskTrackl
##' @export maskTrackll
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

kernelDensity = function (track.list){
  
    #Merge track list into a single dataframe
    df <- mergeTracks(track.list)
    
    #Calculate kernel density from dataframe
    dens <- MASS::kde2d(df[[1]], df[[2]], n=200, lims=c(c(0, 128), c(0, 128)));

  	return (dens);
}

#### createMask ####

#Returns binary mask and plots

createMask = function (track.list, kernel.density, p = NULL, eliminate = NULL, plot = T, separate = F){
	#Store all merged track coordinate points into a dataframe
	df <- mergeTracks(track.list)
	
	if (is.null(p)){
	  p = -0.1207484 + 0.3468734*(nrow(df)/length(track.list))
	}
	
	if (is.null(eliminate)){
	  eliminate = 0
	}
	
	if (p <= 0 || p >= 1){
	  cat("p set to default 0.3.")
	  p = 0.3
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
		plot(df[[2]] ~ df[[1]], col=region, data=df, xlim = c(0, 128), ylim = c(0, 128), xlab = "x", ylab = "y", main = title, cex = .1)
		contour(kernel.density, levels=levels, labels=prob, add=T)
	}
	cat("\n", length(ls), "clusters extracted.", .getTrackFileName(track.list), "masked at a kernel density probability of", round(p, digits = 3), "\n")
  
	if (!separate){
	    return(df$region)
	} else {
	    return(cluster)
	}
}

#### applymask ####

#Creates masked track list

applyMask = function(track.list, mask){
    if (typeof(mask) == "integer"){
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
    } else {
        
    }
}

#### mergeTracks ####

mergeTracks = function(track.list){
  if (length(track.list[[1]]) == 3){
    df <- bind_rows(track.list, .id = "Trajectory")[, c("x", "y", "z")]
  } else {
    df <- bind_rows(track.list, .id = "Trajectory")[, c("x", "y", "z", "Frame")]
  }
  
  return (df);
}

#### maskTrackl ###

maskTrackl = function (track.list, automatic = F){
  cat("\n Mask for", getTrackFileName(track.list), "...\n")
  kd <- kernelDensity(track.list);
  if (automatic){
    mask <- createMask(track.list, kd, p = NULL);
  } else {
    eliminate = 0;
    p = NULL;
    done = FALSE;
    mask <- createMask(track.list, kd, p = p);
    while (!done){
      cat("\n")
      done = as.integer(readline(prompt="Done (1 = YES; 0 = NO)?: "))
      if (is.na(done)){
        cat("\nIncorrect input, set to default = 0.\n")
        done = 0;
      }
      done = as.logical(done)
      if (done){
        break;
      }
      p = as.numeric(readline(prompt="New kernel density probability (p): "))
      if (is.na(p)){
        cat("\nIncorrect input, set to default.\n")
        p = NULL;
      }
      eliminate = as.integer(readline(prompt="Number of smallest clusters to elimnate (recommended 0, unless last resort): "))
      if (is.na(eliminate)){
        cat("\nIncorrect input, set to default = 0.\n")
        eliminate = 0;
      }
      mask <- createMask(track.list, kd, p = p, eliminate = eliminate);
    }
  }
  masked.track.list <- applyMask(track.list, mask);
  cat("\n Mask created for", getTrackFileName(track.list), "\n")
  return(masked.track.list)
}

#### maskTrackll ####

maskTrackll = function (trackll, automatic = F){
  masked.trackll <- list()
  for (i in 1:length(trackll)){
    masked.trackll[[i]] <- maskTrackl(trackll[[i]], automatic = automatic)
  }
  names(masked.trackll) <- names(trackll)
  cat("\nAll tracks lists masked.\n")
  return(masked.trackll)
}

#### plotPoints ####

plotPoints = function(trackll){
  for (i in 1:length(trackll)){
    .plotPoints(trackll[[i]])
  }
}

.plotPoints = function(track.list){
  df <- mergeTracks(track.list)
  plot(df[[1]], df[[2]], xlim = c(0, 128), ylim = c(0, 128), xlab = "x", ylab = "y", main = getTrackFileName(track.list), cex = .1);
}

#### plotLines ####

plotLines = function(trackll){
  for (i in 1:length(trackll)){
    .plotLines(trackll[[i]])
  }
}

.plotLines = function(track.list){
  plot(track.list[[1]][[1]], track.list[[1]][[2]], type = "l", xlim = c(0, 128), ylim = c(0, 128), main = getTrackFileName(track.list))
  for(i in 2:length(track.list)){
    lines(track.list[[i]][[1]], track.list[[i]][[2]])
  }
}

