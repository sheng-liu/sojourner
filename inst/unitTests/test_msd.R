test_msd <- function() {
    folder=system.file("extdata","SWR1",package="sojourner")
    trackll=createTrackll(folder=folder, input=3)
    trackll.flt=filterTrack(trackll,filter=c(min=5,max=Inf))
    msd=msd(trackll.flt,dt=6,summarize=TRUE,plot=TRUE)
    RUnit::checkEquals(msd[[1]][1,][[3]], 45)
    RUnit::checkEqualsNumeric(msd[[1]][1,][[2]], 0.003288, tolerance = 1e-4)
    RUnit::checkEqualsNumeric(msd[[1]][1,][[1]], 0.028727, tolerance = 1e-4)
    RUnit::checkEquals(nrow(msd[[1]]), 6)
    RUnit::checkEquals(ncol(msd[[1]]), 3)
}
