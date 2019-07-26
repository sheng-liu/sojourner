## hmmTrackll-methods
##
##
###############################################################################
##' @name hmmTrackll
##' @aliases hmmTrackll
##' @title hmmTrackll
##' @rdname hmmTrackll-methods
##' @docType methods
##'
##' @description Convert trackll to a (bivarte) displacement format for Hidden 
##' Markov Model fitting.
##' @usage
##' hmmTrackll(trackll, t.interval=0.01)

##' @param trackll list of trajectorise created by crreateTrackll().
##' @param t.interval time interval. 
##' @return Hidden Markov Model fitted trackll

##' @examples
##' folder=system.file("extdata","SWR1",package="sojourner")
##' trackll=createTrackll(folder=folder, input=3)
##' hTrackll=hmmTrackll(trackll)
##' 

##' @details
##' dx, displacement on x.
##' 
##' dy, displacement on y. 
##' 
##' time, lenght of time a particle has being traveling. 
##' 
##' trackStepIndex, steps index within a track. 
##' 
##' 

##' @export hmmTrackll

###############################################################################

##------------------------------------------------------------------------------

# conversion
# > trackll[[1]][[1]]
# x     y z Frame
# 1 40.32 49.36 1     1
# 2 41.03 47.33 1     2
# 3 41.05 48.43 1     3
# 4 40.90 49.09 1     4

# > dat[[1]][[1]]
# time trackStepIndex       dx       dy
# 2 0.01              1  0.07597 -0.21721
# 3 0.02              2  0.00214  0.11770
# 4 0.03              3 -0.01605  0.07062


# TODO: names are frame number, can be remvoed.

##'@export hmmTrackll
##'
##'
##

# data=displacement[["stepwise.displacement"]]

# > displacement[["stepwise.displacement"]]
# $hmmBaysExample
# stepwiseDisplacement trackIndex
# 1            0.21317064      track
# 2            0.20526420      track
# 3            0.09507908      track
# 4            0.16183674      track
# 5            0.15508866      track
# 6            0.11910026      track
# 7            0.03835341      track
# 8            0.07907094      track
# 9            0.19908813      track
# 10           0.29826794      track
# 11           0.05355240      track
# 12           0.13607008      track
# 13           0.05855978      track
# 14           0.05358287      track
# 15           0.31348425      track
# 16           0.03402509      track
# 17           0.04417534      track
# 18           0.08849055      track
# 19           0.09597544      track
# 20           0.09346633      track
# 21           0.08538280      track
# 22           0.08939897      track

hmmTrackll=function(trackll,t.interval=0.01){

    # numericIndex=FALSE
    # all the downstrain hmm analysis based on a bivariate model, so no need to 
    # have a bivar=TRUE paramter
    
    
    cat("calculating displacement...\n")

    # if (bivar == TRUE){
        data=displacement_trackll(trackll,dt=1,bivar=TRUE)
        
    # }else{
    #     dp=displacementCDF(trackll,dt=1,plot=FALSE,bivar=bivar)
    #     data=dp$stepwise.displacement
    # }
    
    #if (numericIndex == FALSE){

        dat=sapply(data,function(x){x["InidvidualDisplacement"]})

        # rownames(dat)
        # 
        # lapply(dat,function(x){
        #     
        #     for (i in seq_along(x)){
        #         y=rownames(x[[i]])
        #     }
        #     y
        # })
        
        rownames(dat[[1]][[1]])
        
        # > rownames(dat[[1]][[1]])
        # [1] "2" "3" "4"
        # > names(dat[[1]][[1]])
        # [1] "dx" "dy"
        # 
        
        names(dat[[1]][1])
        
        
        
        
        # > names(dat[[1]][1])
        # [1] "mage6.1.4.1.1"
        
        
        for (i in seq_along(dat)){

            ## add time
            #  track.len=table(dat[[i]]["trackIndex"])

        #             > head(track.len)
        #
        #             mage6.1.4.1.1 mage6.1003.11.108.108  mage6.1026.2.109.109
        #                         3                    10                     1
        #      mage6.1030.2.110.110  mage6.1042.2.111.111       mage6.11.12.5.5
        #                         1                     1                    11

            track.len=sapply(dat[[i]],function(x){dim(x)[1]})
            
            
            
            # from tb.index generate time
            time.lst=sapply(names(track.len),function(x){
                t=data.frame(seq(from=t.interval,to=track.len[x]*t.interval,
                                 by=t.interval))
                names(t)="time"
                return(t)
            },simplify = FALSE, USE.NAMES = TRUE)

            #         # lapply code is change from this for loop
            #         time=list()
            #         for (j in seq_along(names(track.len))){
            #             time[[j]]=seq(
            #                       from=t.interval,to=track.len[j]*t.interval,
            #                           by=t.interval)
            #         }
            #         names(time)=names(tb.index)
            #         identical(time,time.lst)
            #         [1] TRUE

            # put time into the list, then use do.call(rbind) to make them 
            # a data.frame
            # xx=as.data.frame(mapply(cbind.data.frame,time.lst,dat[[i]]))
            
            # mapply(cbind,time.lst,dat[[i]]) # returns a matrix
            # Map(cbind,time.lst,dat[[i]])
            
            # with Map, no need for the converstion
            # time.df=reshape2::melt(time.lst)
            # names(time.df)=c("time","index")
            # time.df=time.df[order(time.df$index),]

#             # notice the effect of subsetting using ["name"] vs $
#             > str(as.vector(dat[[i]]["trackIndex"]))
#            'data.frame':1249 obs. of  1 variable:
#                 $ trackIndex: chr  "mage6.1.4.1.1" "mage6.1.4.1.1" 
#                                                               "mage6.1.4.1.1"
#             > str(as.vector(dat[[i]]$trackIndex))
#             chr [1:1249] "mage6.1.4.1.1" "mage6.1.4.1.1" "mage6.1.4.1.1" ...

            # track unique id trackStepIndex
            id.lst=sapply(names(track.len),function(x,t.interval=1){
                id=data.frame(seq(from=t.interval,to=track.len[x]*t.interval,
                                  by=t.interval))
                names(id)="trackStepIndex"
                return(id)
            },simplify = FALSE, USE.NAMES = TRUE)

            
            # with Mapp no need for these
            # id.df=reshape2::melt(id.lst)
            # colnames(id.df)=c("stepNum","trackIndex")
            # id.df$trackStepIndex=paste(id.df$trackIndex,id.df$stepNum,sep=".")

            


            # order dat
            # dat[[i]]=dat[[i]][order(dat[[i]]$trackIndex),]
            # head(dat[[i]],20)
            # 
            # time.df=time.df[order(time.df$index),]
            # head(time.df,20)
            # 
            # identical(time.df$index,dat[[i]]$trackIndex)



            ## combine time and dat
            #         dat[[i]]=cbind(dat[[i]],time.df)
            #         dat[[i]]$index=NULL

            # use subsetting saves an operation
            # dat[[i]]$time=time.df$time
            # dat[[i]]$trackStepIndex=id.df$trackStepIndex
            #names(dat[[i]])=c("displacement","mlcnum","time")


            dat[[i]]=Map(cbind,time.lst,id.lst,dat[[i]])
        }

        
    #     }else{
    # 
    # 
    #     dat=list()
    #     length(dat)=length(data)
    #     names(dat)=names(data)
    # 
    #     for (i in seq_along(data)){
    #         ## subtract track index
    #         index.lst=strsplit(data[[i]]$trackIndex,split="[.]")
    #         index=sapply(index.lst,function(x){x[length(x)]})
    #         cat("converting",length(index),"trajectories to panel data\n")
    # 
    #         dat[[i]]=cbind(data[[i]]["StepwiseDisplacement"],index)
    #         colnames(dat[[i]])=c("displacement","mlcnum")
    # 
    # 
    #         ## add time
    #         track.len=table(index)
    #         # > tb.index
    #         # index
    #         # 10 11 12 14 20 22 27  3  4  5  6  8
    #         #  9 17  8  7  6  6  7  8  7 11  8  6
    # 
    # 
    #         # from tb.index generate time
    #         time.lst=sapply(names(track.len),function(x){
    #             seq(from=t.interval,to=track.len[x]*t.interval,by=t.interval)
    #         },simplify = FALSE, USE.NAMES = TRUE)
    # 
    #         #         # lapply code is change from this for loop
    #         #         time=list()
    #         #         for (j in seq_along(names(track.len))){
    #         #             time[[j]]=seq(from=t.interval,to=track.len[j]*
    #         #                           t.interval,
    #         #                           by=t.interval)
    #         #         }
    #         #         names(time)=names(tb.index)
    #         #         identical(time,time.lst)
    #         #         [1] TRUE
    # 
    #         time.df=reshape2::melt(time.lst)
    #         names(time.df)=c("time","index")
    #         time.df=time.df[order(as.integer(time.df$index)),]
    # 
    #         ## combine time and dat
    #         #         dat[[i]]=cbind(dat[[i]],time.df)
    #         #         dat[[i]]$index=NULL
    # 
    #         # use subsetting saves an operation
    #         dat[[i]]$time=time.df$time
    # 
    #     }
    # 
    # }

    return(dat)
}


## convert trackll back to data.frame, with time and mlcnum
## this is useful for performance enhancement next step








##------------------------------------------------------------------------------
## TODO
## in cat information add folder name


