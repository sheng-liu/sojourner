test_filterTrack <- function() {
    folder=system.file("extdata","SWR1",package="sojourner")
    trackll=createTrackll(folder=folder, input=3)
    trackll.flt=filterTrack(trackll,filter=c(min=5,max=Inf))
    RUnit::checkTrue(lapply(trackLength(trackll),min)==2)
    RUnit::checkTrue(lapply(trackLength(trackll.flt),min)==5)
    RUnit::checkTrue(lapply(trackLength(trackll.flt),max)==30)
    RUnit::checkTrue(length(trackll.flt[[1]])==70)
}
