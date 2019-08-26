test_createTrackll <- function() {
    folder=system.file('extdata','SWR1',package='sojourner')
    trackll = createTrackll(folder=folder, input=3)
    csv = read.csv(file = paste(folder,"/SWR1_WT_140mW_image6.csv", sep=""))
    RUnit::checkEquals(trackll[[1]][[1]][1:3], csv[1:4,][4:6])
    RUnit::checkEquals(max(csv[[2]]), length(trackll[[1]]))
    RUnit::checkEquals(csv[csv$Trajectory == 5,4:6]$x, trackll[[1]][[5]][1:3]$x)
    RUnit::checkTrue(is.list(trackll))
    RUnit::checkTrue(is.list(trackll[[1]]))
    RUnit::checkTrue(is.list(trackll[[1]][[1]]))
    RUnit::checkTrue(is.double(trackll[[1]][[1]][[1]]))
}
