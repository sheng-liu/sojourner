test_trimTrack <- function() {
    folder=system.file("extdata","SWR1",package="sojourner")
    trackll=createTrackll(folder=folder, input=3)
    trackll.trm=trimTrack(trackll,trimmer=c(min=1,max=20))
    RUnit::checkTrue(lapply(trackLength(trackll.trm),max)==20)
    RUnit::checkTrue(lapply(trackLength(trackll.trm),min)==2)
    RUnit::checkTrue(lapply(trackLength(trackll),max)==30)
    RUnit::checkTrue(length(trackll.trm[[1]])==207)
}