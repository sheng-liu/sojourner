#### ui.R

# library(smt)
# library(shiny)
# library(shinyjs) 
# this has to be sourced as well to load into memory, even if the package has
# been imported to smt

source("helpers.R")

shinyUI(fluidPage(#theme = "bootstrap.css",
    
    shinyjs::useShinyjs(),
    tags$style(appCSS),
    
    titlePanel(title = "Single Molecule Tracking",
               #title = div("Single Molecule Tracking", img(src = 'university.logo.small.horizontal.blue.png', height = 40, align = "right")), 
               windowTitle = "smt"),
    
    tabsetPanel(id = "section",
                tabPanel("1. Read Tracks",
                         
                         h3("Read Tracks"),
                         
                         sidebarPanel(
                             
                             h3("Select Input Folder"),
                             actionButton(inputId = "folder", 
                                          label = "Select any file in the folder...",
                                          icon = icon("folder-open")),
                             textOutput("folderConfirm"),
                             
                             h3("Set Parameters"),
                             radioButtons(inputId = "input", 
                                          label = h5("Input type:"), 
                                          c("Diatrack .txt file" = 1, 
                                            "Diatrack .mat session file" = 2, 
                                            "ImageJ .csv file" = 3, 
                                            "SlimFast .txt file" = 4), 
                                          selected = 2),
                             
                             checkboxGroupInput(inputId = "parameters", 
                                                label = h5("More options:"), 
                                                choices =  c("Use absolute coordinates" = 1,
                                                             "Keep frame record" = 2),
                                                selected =  2),
                             
                             sliderInput(inputId = "cores", 
                                         label = h5("Cores for parallel computation:"),
                                         min=1, max=parallel::detectCores(logical = F), 
                                         value = 1, step = 1),
                             
                             textOutput("readNote"),
                             
                             withBusyIndicatorUI((actionButton(inputId = "read", 
                                          label = "Read folder",
                                          icon = icon("import", lib = "glyphicon")))),
                             
                             textOutput("readConfirm")
                         )
                         
                         
                ),
                tabPanel("2. Process Tracks",
                         
                         h3("Process Tracks"),
                         
                         sidebarPanel(
                             
                             actionButton(inputId = "reset", 
                                          label = "Reset",
                                          icon = icon("undo")),
                             textOutput("resetConfirm"),
                             
                             h3("Link"),
                             numericInput(inputId <- "tolerance", 
                                          label = h5("Tolerance level (pixels):"),
                                          value = 5),
                             numericInput(inputId <- "maxSkip", 
                                          label = h5("Frame skips:"), 
                                          value = 10,
                                          step = 1),
                             withBusyIndicatorUI(actionButton(inputId = "link", 
                                          label = "Link",
                                          icon = icon("link"))),
                             textOutput("linkConfirm"),
                             
                             h3("Filter"),
                             numericInput(inputId <- "minFilter", 
                                          label = h5("Minimum track length:"),
                                          value = 7),
                             
                             numericInput(inputId <- "maxFilter", 
                                          label = h5("Maximum track length (enter 0 for infinity):"), 
                                          value = 0),
                             withBusyIndicatorUI(actionButton(inputId = "filter", 
                                          label = "Filter",
                                          icon = icon("filter"))),
                             textOutput("filterConfirm"),
                             
                             h3("Trim"),
                             sliderInput(inputId = "trimRange", 
                                         label = h5("Track length range:"),
                                         min = 1, max = 50, 
                                         value = c(1, 10),
                                         step = 1),
                             withBusyIndicatorUI(actionButton(inputId = "trim", 
                                          label = "Trim",
                                          icon = icon("cut"))),
                             textOutput("trimConfirm"),
                             
                             h3("Mask"),
                             withBusyIndicatorUI(actionButton(inputId = "mask", 
                                          label = "Apply image masks",
                                          icon = icon("braille"))),
                             textOutput("maskConfirm"),
                             
                             h3("Merge"),
                             withBusyIndicatorUI(actionButton(inputId = "merge", 
                                          label = "Merge",
                                          icon = icon("compress"))),
                             textOutput("mergeConfirm")
                         )
                         
                ),
                
                navbarMenu("3. Analysis",
                           tabPanel("Mean Squared Displacement",
                                    
                                    h3("Mean Squared Displacement"),
                                    
                                    sidebarPanel(
                                        
                                        h3("Set Parameters"),
                                        
                                        numericInput(inputId <- "dtMSD", 
                                                     label = h5("Time interval (dt): "),
                                                     min = 1,
                                                     value = 6,
                                                     step = 1),
                                        
                                        numericInput(inputId <- "resolutionMSD", 
                                                     label = h5("Resolution (pixels to μM): "),
                                                     min = 0,
                                                     value = 0.107),
                                        
                                        checkboxInput(inputId = "summarizeMSD", 
                                                      label = "Summarize", 
                                                      value = FALSE),
                                        
                                        h3("Output"),
                                        
                                        checkboxInput(inputId = "plotMSD", 
                                                      label = "Plot", 
                                                      value = TRUE),
                                        
                                        checkboxInput(inputId = "outputMSD", 
                                                      label = "Export .csv", 
                                                      value = FALSE),
                                        
                                        withBusyIndicatorUI(actionButton(inputId = "calculateMSD", 
                                                     label = "Calculate MSD",
                                                     icon = icon("stats", lib = "glyphicon"))),
                                        
                                        textOutput("MSDConfirm")
                                    )
                                    
                           ),
                           
                           tabPanel("Diffusion Coefficient",
                                    
                                    h3("Diffusion Coefficient"),
                                    
                                    sidebarPanel(
                                        
                                        h3("Set Parameters"),
                                        
                                        numericInput(inputId <- "dtDcoef", 
                                                     label = h5("Time interval (dt): "),
                                                     min = 1,
                                                     value = 6,
                                                     step = 1),
                                        
                                        numericInput(inputId <- "rsquareDcoef", 
                                                     label = h5("R-squared filter (set to 0 if undesired): "),
                                                     min = 0,
                                                     value = 0.8, 
                                                     step = 0.1),
                                        
                                        numericInput(inputId <- "resolutionDcoef", 
                                                     label = h5("Resolution (pixels to μM): "),
                                                     min = 0,
                                                     value = 0.107, 
                                                     step = 0.001),
                                        
                                        numericInput(inputId <- "binwidthDcoef", 
                                                     label = h5("Bin width (enter 0 for automatic): "),
                                                     min = 0,
                                                     value = 0),
                                        
                                        numericInput(inputId <- "t.intervalDcoef", 
                                                     label = h5("Time interval between frames (s): "),
                                                     min = 0,
                                                     value = 0.01, 
                                                     step = 0.001),
                                        
                                        radioButtons(inputId = "methodDcoef", 
                                                     label = h5("Method:"), 
                                                     c("Static" = 1, 
                                                       "Percentage" = 2, 
                                                       "Rolling window" = 3), 
                                                     selected = 1),
                                        
                                        h3("Output"),
                                        
                                        checkboxInput(inputId = "plotDcoef", 
                                                      label = "Plot", 
                                                      value = TRUE),
                                        
                                        checkboxInput(inputId = "outputDcoef", 
                                                      label = "Export .csv", 
                                                      value = FALSE),
                                        
                                        textOutput("MSDpresent"),
                                        
                                        withBusyIndicatorUI(actionButton(inputId = "calculateDcoef", 
                                                     label = "Calculate Dcoef",
                                                     icon = icon("stats", lib = "glyphicon"))),
                                        
                                        textOutput("DcoefConfirm")
                                    )
                                    
                           ),
                           
                           #tabPanel("Normal Distribution",
                            #        
                            #        h3("Normal Distibution"),
                            #        
                            #        sidebarPanel(
                            #            
                            #            h3("Set Parameters"),
                            #            
                            #            checkboxInput(inputId = "log.transformND", 
                            #                          label = "Logarithmic transformation", 
                            #                          value = FALSE),
                            #            
                            #            numericInput(inputId <- "binwidthND", 
                            #                         label = h5("Bin width (enter 0 for automatic): "),
                            #                         min = 0,
                            #                         value = 0),
                            #            
                            #            h3("Output"),
                            #            
                            #            checkboxInput(inputId = "combine.plotND", 
                            #                          label = "Combine plot", 
                            #                          value = FALSE),
                            #            
                            #            checkboxInput(inputId = "outputND", 
                            #                          label = "Export .csv", 
                            #                          value = FALSE),
                            #            
                            #            textOutput("Dcoefpresent"),
                            #            
                            #            withBusyIndicatorUI(actionButton(inputId = "calculateND", 
                            #                                             label = "Calculate ND",
                            #                                             icon = icon("stats", lib = "glyphicon"))),
                            #            
                            #            textOutput("NDConfirm")
                            #        )
                            #        
                           #),
                           
                           tabPanel("Cummulative Distribution Function",
                                    
                                    h3("Cummulative Distribution Function"),
                                    
                                    sidebarPanel(
                                        
                                        h3("Displacement CDF"),
                                        
                                        numericInput(inputId <- "dtDCDF", 
                                                     label = h5("Time interval (dt): "),
                                                     min = 1,
                                                     value = 6,
                                                     step = 1),
                                        
                                        numericInput(inputId <- "resolutionDCDF", 
                                                     label = h5("Resolution (pixels to μM): "),
                                                     min = 0,
                                                     value = 0.107),
                                        
                                        checkboxInput(inputId = "plotDCDF", 
                                                      label = "Plot", 
                                                      value = TRUE),
                                        
                                        checkboxInput(inputId = "outputDCDF", 
                                                      label = "Export .csv", 
                                                      value = FALSE),
                                        
                                        withBusyIndicatorUI(actionButton(inputId = "calculateDCDF", 
                                                     label = "Calculate displacement CDF",
                                                     icon = icon("stats", lib = "glyphicon"))),
                                        
                                        textOutput("DCDFConfirm"),
                                        
                                        h3("Fit CDF"),
                                        
                                        numericInput(inputId <- "t.intervalFCDF", 
                                                     label = h5("Time interval between frames (s): "),
                                                     min = 0,
                                                     value = 0.01, 
                                                     step = 0.001),
                                        
                                        numericInput(inputId <- "maxiter.searchFCDF", 
                                                     label = h5("Maximum iteration in random search start value process: "),
                                                     min = 0,
                                                     value = 1e3, 
                                                     step = 10),
                                        
                                        numericInput(inputId <- "maxiter.optimFCDF", 
                                                     label = h5("Maximum iteration in local optimization process: "),
                                                     min = 0,
                                                     value = 1e3, 
                                                     step = 10),
                                        
                                        radioButtons(inputId = "componentsFCDF", 
                                                     label = h5("Components:"), 
                                                     c("One" = 1, 
                                                       "Two" = 2, 
                                                       "Three" = 3), 
                                                     selected = 2),
                                        
                                        conditionalPanel(
                                            condition = "input.componentsFCDF == 1",
                                            numericInput(inputId <- "D_1", 
                                                         label = h5("D min value: "),
                                                         value = 1e-3),
                                            
                                            numericInput(inputId <- "D_2", 
                                                         label = h5("D max value: "), 
                                                         value = 2)
                                        ),
                                        
                                        conditionalPanel(
                                            condition = "input.componentsFCDF == 2",
                                            numericInput(inputId <- "D1_1", 
                                                         label = h5("D1 min value: "),
                                                         value = 1e-3),
                                            
                                            numericInput(inputId <- "D1_2", 
                                                         label = h5("D1 max value: "), 
                                                         value = 2),
                                            
                                            numericInput(inputId <- "D2_1", 
                                                         label = h5("D2 min value: "),
                                                         value = 1e-3),
                                            
                                            numericInput(inputId <- "D2_2", 
                                                         label = h5("D2 max value: "), 
                                                         value = 2), 
                                            
                                            numericInput(inputId <- "alpha_1", 
                                                         label = h5("Alpha min value: "),
                                                         value = 1e-3),
                                            
                                            numericInput(inputId <- "alpha_2", 
                                                         label = h5("Alpha max value: "),
                                                         value = 1)
                                        ),
                                        
                                        conditionalPanel(
                                            condition = "input.componentsFCDF == 3",
                                            numericInput(inputId <- "D1_1", 
                                                         label = h5("D1 min value: "),
                                                         value = 1e-3),
                                            
                                            numericInput(inputId <- "D1_2", 
                                                         label = h5("D1 max value: "), 
                                                         value = 2),
                                            
                                            numericInput(inputId <- "D2_1", 
                                                         label = h5("D2 min value: "),
                                                         value = 1e-3),
                                            
                                            numericInput(inputId <- "D2_2", 
                                                         label = h5("D2 max value: "), 
                                                         value = 2), 
                                            
                                            numericInput(inputId <- "D3_1", 
                                                         label = h5("D3 min value: "),
                                                         value = 1e-3),
                                            
                                            numericInput(inputId <- "D3_2", 
                                                         label = h5("D3 max value: "), 
                                                         value = 2),
                                            
                                            numericInput(inputId <- "alpha_1", 
                                                         label = h5("Alpha min value: "),
                                                         value = 1e-3,
                                                         step = 1e-3),
                                            
                                            numericInput(inputId <- "alpha_2", 
                                                         label = h5("Alpha max value: "),
                                                         value = 1), 
                                            
                                            numericInput(inputId <- "beta_1", 
                                                         label = h5("Beta min value: "),
                                                         value = 1e-3,
                                                         step = 1e-3),
                                            
                                            numericInput(inputId <- "beta_2", 
                                                         label = h5("Beta max value: "),
                                                         value = 1)
                                            
                                        ),
                                        
                                        checkboxInput(inputId = "outputFCDF", 
                                                      label = "Export .csv", 
                                                      value = FALSE),
                                        
                                        withBusyIndicatorUI(actionButton(inputId = "calculateFCDF", 
                                                     label = "Fit CDF and export plot .pdf",
                                                     icon = icon("stats", lib = "glyphicon"))),
                                        
                                        textOutput("FCDFConfirm")
                                        
                                    )
                           ),
                           
                           tabPanel("Dwell Time",
                                    
                                    h3("Dwel Time"),
                                    
                                    sidebarPanel(
                                        
                                        h3("Set Parameters"),
                                        
                                        numericInput(inputId <- "t.intervalDT", 
                                                     label = h5("Time interval between frames (ms):"),
                                                     min = 0,
                                                     value = 10,
                                                     step = 1),
                                        
                                        numericInput(inputId <- "x.scale.minDT", 
                                                     label = h5("x-axis scale minimum:"),
                                                     value = 0,
                                                     step = 1),
                                        
                                        numericInput(inputId <- "x.scale.maxDT", 
                                                     label = h5("x-axis scale maximum"),
                                                     value = 250,
                                                     step = 1),
                                        
                                        h3("Output"),
                                        
                                        checkboxInput(inputId = "plotDT", 
                                                      label = "Plot", 
                                                      value = TRUE),
                                        
                                        checkboxInput(inputId = "outputDT", 
                                                      label = "Export .csv and plot .pdf", 
                                                      value = FALSE),
                                        
                                        withBusyIndicatorUI(actionButton(inputId = "calculateDT", 
                                                     label = "Calculate DT",
                                                     icon = icon("stats", lib = "glyphicon"))),
                                        
                                        textOutput("DTConfirm")
                                    )
                           )
                           
                )
                
    ),
    conditionalPanel(condition = "input.section == '1. Read Tracks' || input.section == '2. Process Tracks'", 
                     mainPanel(
                         
                         h3("Track Info"),
                         textOutput("trackllInfo"),
                         textOutput("tracklInfo"),
                         sliderInput(inputId = "tracklNum", 
                                     label = h5("Select video number:"),
                                     min = 1, max = 1,
                                     value = 1, 
                                     step = 1),
                         
                         h3("Export All Tracks"),
                         withBusyIndicatorUI(actionButton(inputId = "export", 
                                      label = "Export to working directory",
                                      icon = icon("download"))),
                         textOutput("exportConfirm"),
                         
                         h3("Track Plot"),
                         radioButtons(inputId = "plotType", 
                                      label = h5("Plot type:"), 
                                      c("Coordinate Points" = 1, 
                                        "Trajectory Lines" = 2),
                                      selected = 2),
                         
                         plotOutput(outputId = "plotPoints", inline = T)
                     )
    ),
    
    conditionalPanel(condition = "input.section == 'Mean Squared Displacement'", 
                     mainPanel(
                         plotOutput(outputId = "plotMSD", inline = T)
                     )
    ),
    
    conditionalPanel(condition = "input.section == 'Diffusion Coefficient'", 
                     mainPanel(
                         plotOutput(outputId = "plotDcoef", inline = T)
                     )
    ),  
    
    conditionalPanel(condition = "input.section == 'Cummulative Distribution Function'", 
                     mainPanel(
                         plotOutput(outputId = "plotDCDF", inline = T)
                     )
    ),  
    
    conditionalPanel(condition = "input.section == 'Dwell Time'", 
                     mainPanel(
                         plotOutput(outputId = "plotDT", inline = T)
                     )
    )
))