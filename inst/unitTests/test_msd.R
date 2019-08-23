test_msd <- function() {
    folder=system.file('extdata','SWR1',package='sojourner')
    trackll = createTrackll(folder=folder, input=3)
    csv = read.csv(file = paste(folder,"/SWR1_WT_140mW_image6.csv", sep=""))
    checkEquals(trackll[[1]][[1]][1:3], csv[1:4,][4:6])
    checkEquals(max(csv[[2]]), length(trackll[[1]]))
    checkEquals(csv[csv$Trajectory == 5,4:6]$x, trackll[[1]][[5]][1:3]$x)
    checkTrue(is.list(trackll))
    checkTrue(is.list(trackll[[1]]))
    checkTrue(is.list(trackll[[1]][[1]]))
    checkTrue(is.double(trackll[[1]][[1]][[1]]))
}
