#### server.R

shinyServer(function(input, output, session){
    
    #Instantiate reactive values
    trackll <- reactiveValues(data= NULL);
    trackll.save <- reactiveValues(data= NULL);
    msd.trackll <- reactiveValues(data= NULL);
    cdf <- reactiveValues(data= NULL);
    fitcdf <- reactiveValues(data= NULL);
    dt <- reactiveValues(data= NULL);
    folder <- reactiveValues(data = NULL);
    nuclei.list <- reactiveValues(data = NULL);
    
    observeEvent(input$folder, {
        folder$data <- dirname(file.choose());
        setwd(folder$data);
        cat("# NEW SESSION CODE: ", Sys.time(), "\n",file="command_history.R", sep = "", append=TRUE);
        cat("folder <- dirname(file.choose());",file="command_history.R", sep = "\n", append=TRUE);
        cat("setwd(folder);",file="command_history.R", sep = "\n", append=TRUE);
        nuclei.list$data = list.files(path=folder$data,pattern="_Nuclei.tif",full.names = TRUE)
        output$folderConfirm <- renderText({
            paste("New working directory: ", folder$data, sep = "")
        })
    })
    
    #Read
    observeEvent(input$read, {
        withBusyIndicatorServer("read", {    
            ab.track = FALSE;
            frameRecord = FALSE;
            for(i in 1:length(input$parameters)){
                if (input$parameters[[i]] == 1){
                    ab.track = TRUE;
                } else if (input$parameters[[i]] == 2) {
                    frameRecord = TRUE;
                }
            }
            trackll$data <- createTrackll(folder = folder$data,
                input = input$input,
                ab.track = ab.track,
                cores = input$cores,
                frameRecord = frameRecord)
            cat("trackll <- createTrackll(folder = folder, input = ", input$input, ", ab.track = ", ab.track, ", cores = ", input$cores, ", frameRecord = ", frameRecord, "); #READ\n", file = "command_history.R", sep = "", append = TRUE);
            trackll.save$data <- trackll$data;
            cat("trackll.save <- trackll",file="command_history.R", sep = "\n", append=TRUE);
            output$readConfirm <- renderText({
                print("Files read.")
            })
        })
    })
    
    output$readNote <- renderText({
        print("Note: reading may take time.")
    })
    
    #Reset
    observeEvent(input$reset, {
        trackll$data <- trackll.save$data;
        cat("trackll <- trackll.save; #RESET",file="command_history.R", sep = "\n", append=TRUE);
        output$resetConfirm <- renderText({
            print("Trackll reset.")
        })
    })
    
    #Link
    observeEvent(input$link, {
        withBusyIndicatorServer("link",{    
            trackll$data <- linkSkippedFrames(trackll = trackll$data, 
                tolerance = input$tolerance, 
                maxSkip = input$maxSkip, 
                cores = input$cores);
            cat("trackll <- linkSkippedFrames(trackll = trackll, tolerance = ", input$tolerance, ", maxSkip = ", input$maxSkip, ", cores = ", input$cores, "); #LINK\n", file = "command_history.R", sep = "", append = TRUE);
            output$linkConfirm <- renderText({
                print("Linking completed.")
            })
        })
    })
    
    #Filter
    observeEvent(input$filter, {
        withBusyIndicatorServer("filter",{
            if (input$maxFilter == 0){
                trackll$data <- filterTrack(trackll$data, 
                    filter = c(min = input$minFilter, max = Inf))
                cat("trackll <- filterTrack(trackll = trackll, filter = c(min = ", input$minFilter, ", max = Inf)); #FILTER\n", file = "command_history.R", sep = "", append = TRUE);
            } else {
                trackll$data <- filterTrack(trackll$data, 
                    filter = c(min = input$minFilter, max = input$maxFilter))
                cat("trackll <- filterTrack(trackll = trackll, filter = c(min = ", input$minFilter, ", max = ", input$maxFilter, ")); #FILTER\n", file = "command_history.R", sep = "", append = TRUE);
            }
            output$filterConfirm <- renderText({
                print("Filtering completed.")
            })
        })
    })
    
    #Trim
    observeEvent(input$trim, {
        withBusyIndicatorServer("trim",{
            trackll$data <- trimTrack(trackll$data, 
                    trimmer = c(min = input$trimRange[[1]], max = input$trimRange[[2]]));
            cat("trackll <- trimTrack(trackll = trackll, trimmer = c(min = ", input$trimRange[[1]], ", max = ", input$trimRange[[2]], ")); #TRIM\n", file = "command_history.R", sep = "", append = TRUE);
            output$trimConfirm <- renderText({
                print("Trimming completed.")
            })
        })
    })
    
    #Mask
    observeEvent(input$mask, {
        withBusyIndicatorServer("mask",{
            trackll$data <- maskTracks(folder$data, trackll$data);
            cat("trackll <- maskTracks(folder = folder, trackll = trackll); #MASK",file="command_history.R", sep = "\n", append=TRUE);
            output$maskConfirm <- renderText({
                print("Masking completed.")
            })
        })
    })
    
    #Merge
    observeEvent(input$merge, {
        withBusyIndicatorServer("merge",{
            trackll$data <- mergeTracks(folder$data, trackll$data)
            cat("trackll <- mergeTracks(folder = folder, trackll = trackll); #MERGE",file="command_history.R", sep = "\n", append=TRUE);
            output$mergeConfirm <- renderText({
                print("Merging completed.")
            })
        })
    })
    
    #Update trackl number slider to match data
    observe({
        updateSliderInput(session, 
            inputId = "tracklNum",
            max = length(trackll$data))
    })
    
    #Plot trackl
    output$plotPoints <- renderPlot({
        if (!is.null(trackll$data)){
            if (input$plotType == 1){
                .plotPoints(trackll$data[[input$tracklNum]])
            } else {
                if (length(nuclei.list$data) > 0){
                    .plotNucTrackOverlay(trackl = trackll$data[input$tracklNum], image.file = nuclei.list$data[input$tracklNum])
                } else {
                    output$noNucOverlay <- renderText({
                        paste("(No nuclei images present for overlay.)")
                    })
                    .plotLines(trackll$data[[input$tracklNum]])
                }
            }
        }
    }, width = 600, height = 600)
    
    #Print trackll info
    output$trackllInfo <- renderText({
        paste("Total number of videos: ", length(trackll$data), sep = " ")
    })
    
    #Print trackl info
    output$tracklInfo <- renderText({
        paste("Number of tracks in video ",  input$tracklNum, ":  ", length(trackll$data[[input$tracklNum]]), sep ="")
    })
    
    #Export current state trackll
    observeEvent(input$export, {
        withBusyIndicatorServer("export",{ 
            exportTrackll(trackll$data, cores = input$cores);
            cat("exportTrackll(trackll = trackll, cores = ", input$cores,"); #EXPORT\n",file="command_history.R", sep = "", append=TRUE);
            output$exportConfirm <- renderText({
                paste("Exported to: ", getwd(), sep = "")
            })
        })
    })
    
    #MSD
    observeEvent(input$calculateMSD, {
        withBusyIndicatorServer("calculateMSD",{ 
            if (input$plotMSD){
                output$plotMSD <- renderPlot({
                    msd.trackll$data <- isolate(msd(trackll$data, 
                        dt = input$dtMSD, 
                        resolution = input$resolutionMSD,
                        summarize = input$summarizeMSD, 
                        cores = input$cores,
                        plot = TRUE,
                        output = input$outputMSD))
                }, width = 900, height = 400)
                cat("msd.trackll <- msd(trackll = trackll, dt = ", input$dtMSD, ", resolution = ", input$resolutionMSD, ", summarize = ", input$summarizeMSD, ", cores = ", input$cores, ", plot = TRUE, output = ", input$outputMSD, "); #MSD\n", file = "command_history.R", sep = "", append = TRUE);
            } else {
                msd.trackll$data <- isolate(msd(trackll$data, 
                    dt = input$dtMSD, 
                    resolution = input$resolutionMSD,
                    summarize = input$summarizeMSD, 
                    cores = input$cores,
                    plot = FALSE,
                    output = input$outputMSD))
                cat("msd.trackll <- msd(trackll = trackll, dt = ", input$dtMSD, ", resolution = ", input$resolutionMSD, ", summarize = ", input$summarizeMSD, ", cores = ", input$cores, ", plot = FALSE, output = ", input$outputMSD, "); #MSD\n", file = "command_history.R", sep = "", append = TRUE);
            }
            if (input$outputMSD){
                output$MSDConfirm <- renderText({
                    paste("MSD calculated. Files exported to: ", getwd(), sep = "")
                })
            } else {
                output$MSDConfirm <- renderText({
                    print("MSD calculated.")
                })
            }
        })
    })
    
    #Dcoef
    observeEvent(input$calculateDcoef, {
        withBusyIndicatorServer("calculateDcoef",{ 
            if (input$methodDcoef == 1){
                method <- "static"
            } else if (input$methodDcoef == 2){
                method <- "percentage"
            } else if (input$methodDcoef == 3){
                method <- "rolling.window"
            }
            
            if (input$binwidthDcoef == 0){
                binwidth <- NULL
            } else {
                binwidth <- input$binwidthDcoef
            }
            
            if (input$plotDcoef){
                output$plotDcoef <- renderPlot({
                    isolate(Dcoef(MSD = msd.trackll$data,
                        trackll = trackll$data, 
                        dt = input$dtDcoef, 
                        rsquare = input$rsquareDcoef, 
                        resolution = input$resolutionDcoef, 
                        binwidth = binwidth, 
                        method = method, 
                        plot = TRUE, 
                        output = input$outputDcoef, 
                        t.interval = input$t.intervalDcoef))
                }, width = 900, height = 400)
                cat("Dcoef(MSD = msd.trackll, trackll = trackll, dt = ", input$dtDcoef, ", rsquare = ", input$rsquareDcoef, ", resolution = ", input$resolutionDcoef, ", binwidth = ", binwidth, ", method = ", method, ", plot = TRUE, output =", input$outputDcoef, ", t.interval = ", input$t.intervalDcoef, "); #Dcoef\n", file = "command_history.R", sep = "", append = TRUE);
                
                updateTabsetPanel(session, "mainTabsetPanel",
                    selected = "Analysis Plots")
                
            } else {
                isolate(Dcoef(MSD = msd.trackll$data,
                    trackll = trackll$data, 
                    dt = input$dtDcoef, 
                    rsquare = input$rsquareDcoef, 
                    resolution = input$resolutionDcoef, 
                    binwidth = binwidth, 
                    method = method, 
                    plot = FALSE, 
                    output = input$outputDcoef, 
                    t.interval = input$t.intervalDcoef))
                cat("Dcoef(MSD = msd.trackll, trackll = trackll, dt = ", input$dtDcoef, ", rsquare = ", input$rsquareDcoef, ", resolution = ", input$resolutionDcoef, ", binwidth = ", binwidth, ", method = ", method, ", plot = FALSE, output =", input$outputDcoef, ", t.interval = ", input$t.intervalDcoef, "); #Dcoef\n", file = "command_history.R", sep = "", append = TRUE);
                
            }
            if (input$outputDcoef){
                output$DcoefConfirm <- renderText({
                    paste("Dcoef calculated. Files exported to: ", getwd(), sep = "")
                })
            } else {
                output$DcoefConfirm <- renderText({
                    print("Dcoef calculated.")
                })
            }
        })
    })
    
    #MSD present notification
    output$MSDpresent <- renderText({
        if (is.null(msd.trackll$data)){
            paste("WARNING: MSD has not been calculated.")
        } else {
            paste("MSD already calculated, ready for diffusion coefficient.")
        }
    })
    
    #Displacement CDF
    observeEvent(input$calculateDCDF, {
        withBusyIndicatorServer("calculateDCDF",{ 
            if (input$plotDCDF){
                output$plotDCDF <- renderPlot({
                    cdf$data <- isolate(displacementCDF(trackll = trackll$data,
                        dt = input$dtDCDF,
                        resolution = input$resolutionDCDF, 
                        plot = TRUE,
                        output = input$outputDCDF))
                }, width = 900, height = 600)
                cat("cdf <- displacementCDF(trackll = trackll, dt = ", input$dtDCDF, ", resolution = ", input$resolutionDCDF, ", plot = TRUE, output = ", input$outputDCDF, "); #DCDF\n", file = "command_history.R", sep = "", append = TRUE);
                
                updateTabsetPanel(session, "mainTabsetPanel",
                                  selected = "Analysis Plots")
                
            } else {
                cdf$data <- isolate(displacementCDF(trackll = trackll$data,
                    dt = input$dtDCDF,
                    resolution = input$resolutionDCDF, 
                    plot = FALSE,
                    output = input$outputDCDF))
                cat("cdf <- displacementCDF(trackll = trackll, dt = ", input$dtDCDF, ", resolution = ", input$resolutionDCDF, ", plot = FALSE, output = ", input$outputDCDF, "); #DCDF\n", file = "command_history.R", sep = "", append = TRUE);
                
                
            }
            if (input$outputDCDF){
                output$DCDFConfirm <- renderText({
                    paste("Displacement CDF calculted. Files exported to: ", getwd(), sep = "")
                })
            } else {
                output$DCDFConfirm <- renderText({
                    print("Displacement CDF calculated.")
                })
            }
        })
    })
    
    #Fit CDF 
    observeEvent(input$calculateFCDF, {
        withBusyIndicatorServer("calculateFCDF",{ 
            if (isolate(input$componentsFCDF) == 1){
                fitcdf$data <- isolate(fitCDF(cdf = cdf$data, 
                    components="one",
                    start=list(
                        oneCompFit=list(D=c(input$D_1,input$D_2))
                    ),
                    t.interval=input$t.intervalFCDF,
                    maxiter.search=input$maxiter.searchFCDF,
                    maxiter.optim=input$maxiter.optimFCDF,
                    output = input$outputFCDF,
                    seed=NULL))
                cat("fitdf <- fitCDF(cdf = cdf, components=\"one\", start=list(oneCompFit=list(D=c(",input$D_1,",",input$D_2,"))), t.interval=",input$t.intervalFCDF,", maxiter.search=",input$maxiter.searchFCDF,", maxiter.optim=",input$maxiter.optimFCDF,", output = ",input$outputFCDF,", seed=NULL) #fitCDF\n", file = "command_history.R", sep = "", append = TRUE);
            } else if (isolate(input$componentsFCDF) == 2){
                fitcdf$data <- isolate(fitCDF(cdf = cdf$data, 
                    components="two",
                    start=list(
                        twoCompFit=list(D1=c(input$D1_1,input$D1_2),
                            D2=c(input$D2_1,input$D2_2),
                            alpha=c(input$alpha_1,input$alpha_2))
                    ),
                    t.interval=input$t.intervalFCDF,
                    maxiter.search=input$maxiter.searchFCDF,
                    maxiter.optim=input$maxiter.optimFCDF,
                    output = input$outputFCDF,
                    seed=NULL))
                cat("fitdf <- fitCDF(cdf = cdf, components=\"two\", start=list(twoCompFit=list(D1=c(",input$D1_1,",",input$D1_2,"),D2=c(",input$D2_1,",",input$D2_2,"),alpha=c(",input$alpha_1,",",input$alpha_2,"))), t.interval=",input$t.intervalFCDF,", maxiter.search=",input$maxiter.searchFCDF,", maxiter.optim=",input$maxiter.optimFCDF,", output = ",input$outputFCDF,", seed=NULL) #fitCDF\n", file = "command_history.R", sep = "", append = TRUE);
            } else if (isolate(input$componentsFCDF) == 3){
                fitcdf$data <- isolate(fitCDF(cdf = cdf$data, 
                    components="three",
                    start=list(
                        threeCompFit=list(D1=c(input$D1_1,input$D1_2), 
                            D2=c(input$D2_1,input$D2_2), 
                            D3=c(input$D3_1,input$D3_2), 
                            alpha=c(input$alpha_1,input$alpha_2),
                            beta=c(input$beta_1,input$beta_2))
                    ),
                    t.interval=input$t.intervalFCDF,
                    maxiter.search=input$maxiter.searchFCDF,
                    maxiter.optim=input$maxiter.optimFCDF,
                    output = input$outputFCDF,
                    seed=NULL))
                cat("fitdf <- fitCDF(cdf = cdf, components=\"three\", start=list(threeCompFit=list(D1=c(",input$D1_1,",",input$D1_2,"),D2=c(",input$D2_1,",",input$D2_2,"), D3=c(",input$D3_1,",",input$D3_2,"), alpha=c(",input$alpha_1,",",input$alpha_2,", beta=c(",input$beta_1,",",input$beta_2,"))))), t.interval=",input$t.intervalFCDF,", maxiter.search=",input$maxiter.searchFCDF,", maxiter.optim=",input$maxiter.optimFCDF,", output = ",input$outputFCDF,", seed=NULL) #fitCDF\n", file = "command_history.R", sep = "", append = TRUE);
            }
            if (input$outputFCDF){
                output$FCDFConfirm <- renderText({
                    paste("Fit CDF calculated. Files/plots exported to: ", getwd(), sep = "")
                })
            } else {
                output$FCDFConfirm <- renderText({
                    paste("Fit CDF calculated. Plots exported to: ", getwd(), sep = "")
                })
            }
        })
    })
    
    #Dwell Time
    observeEvent(input$calculateDT, {
        withBusyIndicatorServer("calculateDT",{ 
            if (input$plotDT){
                output$plotDT <- renderPlot({
                    dt$data <- isolate(dwellTime(trackll$data,
                        t.interval = input$t.intervalDT, 
                        x.scale = c(min = input$x.scale.minDT, max = input$x.scale.maxDT), 
                        plot = TRUE, 
                        output = input$outputDT))
                }, width = 900, height = 600)
                cat("dt <- dwellTime(trackll = trackll, t.interval = ", input$t.intervalDT,", x.scale = c(min = ",input$x.scale.minDT,", max = ",input$x.scale.maxDT,"), plot = TRUE, output = ",input$outputDT,"); #DWELLTIME\n", file = "command_history.R", sep = "", append = TRUE);
                updateTabsetPanel(session, "mainTabsetPanel",
                    selected = "Analysis Plots")
                
            } else {
                dt$data <- isolate(dwellTime(trackll$data,
                    t.interval = input$t.intervalDT, 
                    x.scale = c(min = input$x.scale.minDT, max = input$x.scale.maxDT), 
                    plot = FALSE, 
                    output = input$outputDT))
                cat("dt <- dwellTime(trackll = trackll, t.interval = ", input$t.intervalDT,", x.scale = c(min = ",input$x.scale.minDT,", max = ",input$x.scale.maxDT,"), plot = FALSE, output = ",input$outputDT,"); #DWELLTIME\n", file = "command_history.R", sep = "", append = TRUE);
            }
            if (input$outputDT){
                output$DTConfirm <- renderText({
                    paste("DT calculated. Files exported to: ", getwd(), sep = "")
                })
            } else {
                output$DTConfirm <- renderText({
                    print("DT calculated.")
                })
            }
        })
    })

})
