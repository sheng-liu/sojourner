## diffusionMap
##
##
###############################################################################
##' @name diffusionMap
##' @aliases diffusionMap
##' @title diffusionMap
##' @rdname diffusionMap
##' @docType methods

##' @description Take in a .csv trackll file and perform diffusion map analysis.

##' @usage 
##' diffusionMap(filename)

##' @param filename trackll file path in .csv

##' @details
##' 

##' @examples
##' 

##' @export diffusionMap

##' @importFrom reticulate

##------------------------------------------------------------------------------
##
diffusionMap=function(filename){
    # Force python3 use... compatibility with Windows?
    # use_python('/usr/local/bin/python3', required = T)
    diffusionMap_py=system.file("python", "diffusionMap.py", package="sojourner")
    source_python(diffusionMap_py)
    diffusion_map(filename)
}