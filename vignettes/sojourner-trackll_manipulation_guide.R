## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- eval = FALSE-------------------------------------------------------
#  #Basic function call with interactive menu (optimzing 2 cores)
#  trackll <- createTrackll(interact = T, cores = 2)
#  
#  #Manual function call to process Diatrack session files (.mat)
#  trackll <- createTrackll("/DIRECTORYPATH/", input = 2, cores = 2)

## ---- eval = FALSE-------------------------------------------------------
#  #Basic function call to exportTrackll with 2 cores into current directory
#  exportTrackll(trackll, cores = 2)
#  
#  #Export one track list
#  .exportRowWise(trackl)
#  
#  #Get current working directory
#  getwd()
#  
#  #Import export save back into a trackll
#  trackll.2 <- createTrackll(folder = getwd(), input = 3, cores = 2)

## ---- eval = FALSE-------------------------------------------------------
#  #Basic function call of linkSkippedFrames
#  trackll.linked <- linkSkippedFrames(trackll, tolerance = 5, maxSkip = 10)
#  
#  #Export linked trackll into .csv files
#  exportTrackll(trackll.linked, cores = 2)

