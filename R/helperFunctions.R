## helperFunctions
#
################################################################################

##------------------------------------------------------------------------------
## .timeStamp
# add time stamp and file name as a unique signature of the output file
##@export .timeStamp
.timeStamp=function(filename){

    basename=basename(filename)
    name=unlist(strsplit(basename,split="[.]"))
    fileName=paste(name[1],"-",format(Sys.time(),"%Y%m%_%H%M%S"),sep="")
    return(fileName)

}

##------------------------------------------------------------------------------
## .valid

# validity check for dt less than track length (-1)
.valid=function(dt,track){

    # get track length of all tracks
    tracklen=dim(track)[1]

    if (dt >(tracklen-1)){
        stop("\ntrack length:\t",dim(track)[1],
                "\ndt:\t\t",dt,
                "\nTime interval (dt) greater than track length-1\n")
    }
}


.valid=function(dt,track){

    # get track length of all tracks
    tracklen=dim(track)[1]

    if (dt >(tracklen-1)){
        stop("\ntrack length:\t",dim(track)[1],
                "\ndt:\t\t",dt,
                "\nTime interval (dt) greater than track length-1\n")
    }
}






##------------------------------------------------------------------------------
## tracks.msda2sjr

## @export tracks.msda2sjr
tracks.msda2sjr=function(file){

    tracks.file=readMat(file)
    # file.name=basename(file)
    tracks.mat=tracks.file$tracks

    trackl.sjr=lapply(tracks.mat, function(x){
        x=data.frame(x)
        x=x[,-1] # remove time column
        x=x/0.107  # change micrometer to pixel
        colnames(x)=c("x","y")
        x$z=rep(1,times=dim(x)[1])
        return(x)
    })

    trackll.sjr=list(trackl.sjr)
    names(trackll.sjr)=basename(file)
    return(trackll.sjr)

}

# convert msda to sojourner format
# not used in functions yet

##------------------------------------------------------------------------------
## same.scale

same.scale=function(mixmdl.lst){

    scale=list()
    length(scale)=length(mixmdl.lst)
    names(scale)=names(mixmdl.lst)

    for (i in 1:length(mixmdl.lst)){

        den=density(mixmdl.lst[[i]]$x)
        scale.x=c(min=min(den$x),max=max(den$x))
        scale.y=c(min=min(den$y),max=max(den$y))
        scale[[i]]=data.frame(scale.x,scale.y)
    }

    scale.df=do.call(rbind.data.frame,scale)
    scale.same=data.frame(apply(scale.df,2,function(x){
        c(min=min(x),max=max(x))}))

    return(scale.same)

}

##------------------------------------------------------------------------------
## getStartFrame
# returns starting frame of a track/trajectory (using its name) at a given index
# for a track list

#PARAMETERS:
#track.list = named track list output
#Note: Last five characters of the original file name without extension (cannot 
#contain ".")
#index = index of the track in the track list (track number)

##@export getStartFrame
getStartFrame = function(track.list, index){
    return(as.numeric(substr(names(track.list[index]),
                                gregexpr(pattern = '\\.', 
                                         names(track.list[index]))[[1]][1]+1,
                                gregexpr(pattern = '\\.', 
                                         names(track.list[index]))[[1]][2]-1)));
}

##------------------------------------------------------------------------------
## getTrackFileName
# returns the shortened file name of the track

#PARAMETERS:
#track.list = named track list output
#Note: Last five characters of the original file name without extension (cannot 
#contain ".")

##@export getTrackFileName
getTrackFileName = function(track.list){
    return(substr(names(track.list[1]), 1, 
                  gregexpr(pattern = '\\.', names(track.list[1]))[[1]][1]-1));
}

##------------------------------------------------------------------------------
## abTrack
# returns absolute corrdinates of a track (for plotting)

#PARAMETERS:
#track = track input to be transformed into absolute coordinates

##@export abTrack
abTrack = function(track){
    min.x = min(track$x);
    min.y = min(track$y);
    return(data.frame(x=track$x-min.x, y=track$y-min.y));
}
##------------------------------------------------------------------------------
## removeFrameRecord
# remove frame record for backwards compatibility

#PARAMETERS:
#track.list = track list with frame record in the fourth column

##@export removeFrameRecord
removeFrameRecord = function(track.list){
    for (i in 1:length(track.list)){
        track.list[[i]] <- track.list[[i]][-c(4)];
    }
    return (track.list);
}

