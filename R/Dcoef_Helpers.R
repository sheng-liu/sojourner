## Dcoef helpers

## < r^2 > = 2 * d * D * △t = 4 D △t

## < r^2 >, mean square displacement
## d, dimensionality of the problem
## D, diffusion coefficient
## t, time

## D △t = < r^2 > / 4


# divide the values by 2 and 2
# first 2 comes from Einstein Brownian motion equation k = sqrt(2 * D * dT)
# second 2 comes from dimension

##-----------------------------------------------------------------------------
## Dcoef.static

## return a list of coefficients
Dcoef.static=function(MSD,lag.start=2,lag.end=5,t.interval=0.010){

    # linear fitting of the MSD curves between dt 2 and 5
    cat("lag.start  ",lag.start,"\t","lag.end  ",lag.end,"\n")

    # specifying x
    x=(lag.start:lag.end)*t.interval
    dimension=2

    # fitting y to x
    D.coef=lapply(MSD,function(msd){
            apply(msd[lag.start:lag.end,],MARGIN=2,function(y){
                fit=lm(y~x)
                MSDslope=coefficients(fit)[2]/(2*dimension)
                MSDcorr=summary(fit)$r.squared
                sc=c(MSDslope,MSDcorr)
                names(sc)=c("slope","corr")
                return(sc)
            })
        })
    # this is the lapply-for setup, see below for alternative for-lapply setup
    # and its comparison

    # change shape of the matrix
    D.coef=sapply(D.coef,function(x){t(x)},simplify = FALSE)

    return(D.coef)
}

# see below for alternative for-lapply setup; reverse lapply-for setup is
# better, list are named and can easily be paralleled; the cons being it can
# only pass in one parameter, if more than two is not as easily setup

# for-lapply setup
#     D.coef=list()
#     for (i in 1:length(MSD)){
#         D.coef[[i]]=apply(MSD[[i]][lag.start:lag.end,],MARGIN=2,function(y){
#             fit=lm(y~x)
#             MSDslope=coefficients(fit)[2]/2/dimension
#             MSDcorr=summary(fit)$r.squared
#             sc=c(MSDslope,MSDcorr)
#             names(sc)=c("slope","corr")
#             return(sc)
#         })
#
#     }
#     names(D.coef)=names(MSD)

# change shape of the matrix
#     D.coef=sapply(D.coef,function(x){
#         x=t(x)
#         # colnames(x)=c("slope","corr")
#         return(x)
#     },simplify = FALSE)



##------------------------------------------------------------------------------
## .Dcoef.roll

## cant use roll on MSD method = percentage, as its MSD is different length, MSD method = percentage is calculated differently, using msd_track_vecdt(), instead of msd().

Dcoef.roll=function(MSD,window.size=4,t.interval=0.010){

    D.coef=list()
    D.coef.roll=list()
    names.D.coef.roll=c()
    window=1:window.size
    dt=dim(MSD[[1]])[1]
    dimension=2


    for (i in 1:length(MSD)){

        # zero stands for the first window as subsetting using window+j
        for ( j in 0:(dt-window.size)){
            print(window+j)
            x=window+j
            D.coef.roll[[j+1]]=apply(MSD[[i]][window+j,],
                                MARGIN=2,
                                function(y){
                                    x=x*t.interval
                                    fit=lm(y~x)
                                    MSDslope=coefficients(fit)[2]/2/dimension
                                    MSDcorr=summary(fit)$r.squared
                                    sc=c(MSDslope,MSDcorr)
                                    names(sc)=c("slope","corr")
                                    return(sc)
                                    })
            names.D.coef.roll=c(names.D.coef.roll,
                                paste(as.character(window+j),collapse=" "))
        }

        names(D.coef.roll)=names.D.coef.roll

        D.coef[[i]]=D.coef.roll
        names.D.coef.roll=c()
    }
    names(D.coef)=names(MSD)


    # change shape of the matrix
    for (i in 1:length(D.coef)){
        for (j in 1:length(D.coef[[i]])){
            D.coef[[i]][[j]]=t(D.coef[[i]][[j]])
        }
    }

#     # alternative
#     D.coef=lapply(D.coef,function(x){
#         sapply(x,function(x){t(x)},simplify = FALSE)
#
#     })

    return(D.coef)
}

# alternative treat Dcoef.roll as special case of Dcoef.static, 
# except start and end is rolling

# D.coef.roll=function(MSD,window.size=4,t.interval=0.010){
# 
#     dt=dim(MSD[[1]])[1]
# 
#     start=1:(dt-window.size+1)
#     end=lag.start+window.size-1
# 
#     D.coef=list()
#     D.coef.vec=c()
#     for (i in start){
#         dd=Dcoef.static(MSD,lag.start=start[i],lag.end=end[i],t.interval=0.010)
#         D.coef.vec=c(D.coef.vec,list(dd))
#     }
# 
# 
#     # rearrange
#     D.coef.vec
# 
#     mapply(c,start,end)
# 
#     # generate names
#     name.vec=c()
#     for (i in start) {
#         name=paste(start[[i]]:end[[i]],collapse=" ")
#         name.vec=c(name.vec,name)
#     }
# 
#     names(D.coef.vec)=name.vec
# 
# 
# #n=c()
#     # take the first item of the list to form a new list
# #     for (i in start){
# #         for (file.number in 1:length(MSD)){
# #             n=c(n,list(D.coef.vec[[i]][file.number]))
# #         }
# #
# #     }
# #
# #     n=c()
# #     for (file.number in 1:length(MSD)){
# #
# #         for (i in start){
# #
# #             n=c(n,list(D.coef.vec[[i]][file.number]))
# #
# #         }
# #     }
# #
# 
# 
# 
#     name2=names(MSD)
#     l=list()
#     for (j in name2 ){
#         ind=which(name2 == j)
#         l[[ind]]=lapply(D.coef.vec,function(x,i){x[i]},i=j)
#     }
# 
#     names(l)=name2
# 
#     # remove last level of list
#     z=l
#     for (i in 1:length(l)){
#         for (j in 1:length(l[[i]])){
#             z[[i]][j]=l[[i]][[j]]
#         }
#     }
# 
#     #z[[1]][[1]]=l[[1]][1]
# 
# }





##------------------------------------------------------------------------------
## percentage
## To determine the diffusion constant from a trajectory, a line was fit to 
## MSD(n􏲄t) with n running from 1 to the largest integer less than or equal t
## o L/4 (Saxton, 1997).

## "The short-range diffusion coefficients D*(0: 4) and D*(O: 8) are well 
## determined; the longest-range diffusion coefficients D*(0: 512) and 
## D*(0: 1024) are so broadly dis- tributed as to be useless (Fig. 2 a)" 
## (Saxton, 1997)
## from this the maximum track length used for determining diffusion 
## coefficient should be restricted to 32, which yeilds 1/4*32=8. 
## We can then use the trimTrack() to realize this cut-off.

## 0~4 frames equals 4 steps, so the minimum frame taken into account in
## sojourner's numbering system (which start with frame 1 rather than 0) 
## should be 1~5, D(1:5) and D(1:9) to allow sampling from 0~4. Anything 
## below 8, should be directly using percentage =1, what tracks 
## that has length 9

# "D*(2: 4), the short-range diffusion coefficient used by Kusumi et al.(1993).
# A short-range D* has the advantages that it is accurately obtained and the
# influence of directed and confined motion is minimized. D*(2: 4) advantageous
# for analyzing experi- mental data. The range of D is wide enough that it is
# convenient to plot the distribution of log D. (Saxton, 1997)

# any track length >32, take 32
# 22~32, percentage 1/4, use first 5~8 points, excluding initial point (ie.2~6-2~8);
# 15~21, percentage 0.4, use first 5~8 points, excluding initial point (ie.2~6-2~8);
# 10~14, percentage 0.6, use first 5~8 points, excluding initial point (i.e.2~6-2~8);(1)
# 7~9, percentage 1, use all 5~9 frames, excluding initital point (i.e dt 2~6-2~8).(2)
# 5~6, percentage 1, use all points (i.e. dt 2~4-2~5)(3)

# (1) 10*0.6 = 6 frames, 6 time lags, exclude 1st, 5 points for fitting
# (2) 7 frames, 6 time lags, exclude 1 st, 5 points for fitting.
# (3) 5 frames, 4 time step, remove 2, left 2, use all points if less than 7

# this tired percentage setup makes it always use first 2~4-2~8 points 
# for fitting and deriving diffusion coefficient.



Dcoef.perc=function(trackll,percentage=0.25,weighted=FALSE,filter=c(min=7,max=Inf),
                    trimmer=c(min=1,max=31),resolution=0.107,t.interval=0.010){

    # calculate msd using msd_perc()
    msd.lst=msd_perc(trackll,percentage=percentage,filter=filter,trimmer=trimmer,
                    resolution=resolution,output=FALSE)

    # exclude the first time lag for fitting for all category
    msd.remove1st=lapply(msd.lst,function(x){
        for (i in 1:length(x)){x[[i]]=x[[i]][-1]}
        return(x)
    })


#     # the reverse setup is not as convenient
#     msd.remove2st=list()
#     for (i in 1:length(msd.list)){
#         msd.remove2st[[i]]=lapply(msd.list[[i]],function(x){
#             x=x[-1]
#         })
#     }
#     names(msd.remove2st)=names(msd.lst)


    # copy trackll's structure
    D.coef=list()
    length(D.coef)=length(msd.lst)
    names(D.coef)=names(msd.lst)


    dimension=2


    ## this works for trajecotries length >=7 frames
    for (i in 1:length(msd.remove1st)){

        for (j in 1:length(msd.remove1st[[i]])){

            y=msd.remove1st[[i]][[j]]
            len=length(msd.remove1st[[i]][[j]])

            # exclude 1st x value for fitting
            x=seq(from=t.interval*2,to=(len+1)*t.interval,by=t.interval)

            if (weighted == TRUE){
                w=1:len; fit=lm(y~x,weights =w )
            }else{
                fit=lm(y~x)
            }

            MSDslope=coefficients(fit)[2]/2/dimension
            MSDcorr=summary(fit)$r.squared
            sc=c(MSDslope,MSDcorr)
            names(sc)=c("slope","corr")
            D.coef[[i]][[j]]=sc
        }

    }


    # lapply can't be used as need two variable function
    # however you can use lapply and for loop inside

#     for (i in 1:length(y)){
#         lapply(y[[i]],function(x,t.interval){
#             len=length(x)
#
#         })
#    }


    # this changes shape/format into a matrix
    D.coef=sapply(D.coef,function(x){
        do.call(rbind,x)},simplify=FALSE)

    return(D.coef)

}


# # an attempt to incorperate trajecotry has length less than 6
#     # remove first element of msd.lst
#     msd.remove1st=lapply(msd.lst,function(x){
#         for (i in 1:length(x)){
#
#             # if 3=<frames<=6 are selected, all points needs to be used
#
#             if (length(x[[i]])<=3) {stop(
#                 "tracks must at least have 3 time steps (i.e. 4 frames) for confident coefficient fitting, please filter track first.")
#
#             }else if (length(x[[i]])>=5 & length(x[[i]])<=6){
#                 x[[i]]=x[[i]]
#             }else{
#                 x[[i]]=x[[i]][-1]
#             }
#
#         }
#         return(x)
#     })
#






# x=lst()
# for (i in 1:length(D.coef)){
#     x[[i]]=do.call(rbind,D.coef[i])
# }
#
# dx=sapply(D.coef,function(x){
#     do.call(rbind,x)})
# # this returns a matrix
#
# }
# names(D.coef)=names(MSD)


#     Weights are set to be the number of points (length of trajectory?) averaged to generate the mean square displacement value at the given delay (in this case, it is the 25%). Thus, we give more weight to MSD curves with greater certainty (larger number of elements averaged).

# weights are essentially the lenght of the msd list

#     % - M the weighted mean of MSD for each delay
#     % - STD the weighted standard deviation
#     % - N the number of degrees of freedom in the weighted mean
#     % (see http://en.wikipedia.org/wiki/Weighted_mean)

# plot those coef and get the mean of all

# The only requirement for weights is that the vector supplied must be the same length as the data.

# simplest weights  index of the msd (as it shows how many points is used to generate the msd, steps)

# more sophistacted  1/theta^2 (variance)



# would be nice to have a subsetting method for a sojourner class, there are so many levels of subsetting, each time it needs a lapply

## instead of trackll, it maybe better to store tracks in data.table, then folders
## instead of (second) list of data.frame, it may worth the effort simply making it a data.frame (data.table) with fourth column as trajectory numbers.

## it makes program (maybe) easier, computation faster

##------------------------------------------------------------------------------
## rsquare.filter
## r.squared >= rsquaae as quality control

##export rsquare.filter
rsquare.filter=function(D.coef,rsquare=0.8){

    cat("\nApplying r square filter...",rsquare,"\n")

    slope=lapply(D.coef,function(x){x[,"slope"]})   # x[colnames(x) == "slope"]
    corr=lapply(D.coef,function(x){x[,"corr"]})    # x[colnames(x) == "corr"]

    # the "still" molecule wil generate a NA in correlation, thus is.na(x) == FALSE
    corr.filter=lapply(corr,function(x){x>=rsquare & is.na(x) == FALSE})

    # add corr and slope in the output

    # mapply("[",D.coef,corr.filter,SIMPLIFY=FALSE)
    # directly mapply to D.coef, lost the matrix structure

    D.coef.slope.subset=mapply("[",slope,corr.filter,SIMPLIFY=FALSE)
    D.coef.corr.subset=mapply("[",corr,corr.filter,SIMPLIFY=FALSE)
    D.coef.subset=mapply(cbind,D.coef.slope.subset,D.coef.corr.subset,SIMPLIFY=FALSE)

    # add colnames
    D.coef.subset=lapply(D.coef.subset,function(x){
        colnames(x)=c("slope","corr")
        return(x)
    })

    return(D.coef.subset)
}

rsquare.filter.roll=function(D.coef,rsquare=0.8){
    D.coef.subset=list()
    length(D.coef.subset)=length(D.coef)
    names(D.coef.subset)=names(D.coef)
    for (i in 1:length(D.coef)) {
        for(j in 1:length(D.coef[[i]])){
            D.coef.subset[[i]][j]=rsquare.filter(D.coef[[i]][j],
                                                rsquare=rsquare)
            #names(D.coef.subset[[i]][j])=names(D.coef[[i]][j])

        }
        # increase a level to name the list
        names(D.coef.subset[[i]])=names(D.coef[[i]])
    }
    return(D.coef.subset)
}




# to varify the fit
# fit=lm(MSD[[1]][2:5,][,1]~x); plot(fit)

# DONE: output goodness of fit


## the next two filterTrack blocks maybe combined to increase efficiency,
## however for now the efficiency is secondary, let the logic stand
## clear, then improve the efficiency. as many times efficiency is at
## the expense of sacrifice clearness of the code, hard to read or
## interpretate later.


## alternative  works
#             for (m in 1:length(corr)){
#                 for (n in 1:length(corr[[m]])){
#                     if (corr[[m]][n]<rsquare)
#                         slope[[m]][n]=NaN
#                 }
#             }

# D.coef.subset[[i]]=slope



## filter without losing location information
## or filter in the last step, has to be replaced before log
## if corr <rsquare, replace slope with NaN

## alternative
#     for (i in 1: length(D.coef)){
#
#         for (j in 1: length(D.coef[[i]])){
#
#             # dim(D.coef[[i]][[j]])[2] is the length of the matrix
#             for (k in 1:dim(D.coef[[i]][[j]])[2]){
#
#
#                 if (D.coef[[i]][[j]][,k]["corr"]<rsquare)
#                     D.coef[[i]][[j]][,k]["slope"]=NaN
#
#             }
#         }
#     }







##------------------------------------------------------------------------------
## Dcoef.log
##@export Dcoef.log
# Dcoef.log=function(D.coef.subset,static=TRUE){
#     if (static){
#
#         #Log.D.coef=suppressWarnings(lapply(D.coef,log))
#         # worth noting "log computes logarithms, by default natural logarithms"
#         Log.D.coef=lapply(D.coef.subset,log10)
#
#         # remove NaN if wanted
#         # Log.D.coef=lapply(Log.D.coef, function(x){
#         #    x[!is.nan(x)]
#         #})
#     }else{
#         ## logorithm
#         Log.D.coef=list()
#         for (i in 1:length(D.coef.subset)){
#             #Log.D.coef=suppressWarnings(lapply(D.coef,log))
#             Log.D.coef[[i]]=lapply(D.coef.subset[[i]],log)
#         }
#         names(Log.D.coef)=names(D.coef.subset)
#         return(Log.D.coef)
#
#     }
#     return(Log.D.coef)
#
# }

# not used but keep
Dcoef.log=function(D.coef.subset){


        #Log.D.coef=suppressWarnings(lapply(D.coef,log))
        # worth noting "log computes logarithms, by default natural logarithms"
        Log.D.coef=lapply(D.coef.subset,log10)

        # remove NaN if wanted
        # Log.D.coef=lapply(Log.D.coef, function(x){
        #    x[!is.nan(x)]
        #})

    return(Log.D.coef)
}
