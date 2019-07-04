# Description ------------------------------------

# creates peak image matrices for a classifier and plots for visual review

# User Input ------------------------------------------------------------------


# for Windows operating system use \\ or / instead of \ to separate directories
# example:
# outPutParentDir <- 'M:\\folder1\\foldwerWithFile'

# path to folder where output directory will be created
outPutParentDir <- "M:\\folder1\\foldwerWithFile"

# csv filename with peak windows in parent directory
windowCSV <- "filename1.csv"

# Path to folder with ms1 tables
ms1TableDir <- "M:\\folder2\\foldwerWithFile"

# Path to max intensity table
maxIntensityDir <- "M:\\folder3\\foldwerWithFile"

# csv filename of max intensity table
csvMaxIntensity <- "filename2.csv"

# retention time duration to include in each peak image
# must be longer than the widest peak window
# in units of retention time
RtDur <- 1

# max number of injection files to include in matrix
# if fewer than this number exist, remaining rows will be filled with 0s
maxSignalFiles <- 30

# minimum scan separation in units of retention time
minRtSep <- 0.005

# proportion of processors to use 0-1
propProcessors <- 0.75

# Restart previous run?
restart <- FALSE

# plot output folder path
existingOutputPath <- NA


# End user input ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

# advanced user input (probably don't need to change) ------------------------


# column numbers for ID, min m/z, max m/z, min RT, max RT
# adjust if using non-standard window table
colNumWindowID <- 1
colNumWindMzMin <- 2
colNumWindMzMax <- 3
colNumWindRtMin <- 4
colNumWindRtMax <- 5

# proportion of plot lines to take from highest peaks in window
propHighPeaks <- 0.7


# Library (packages)--------------------------------------------------------
library(data.table) # faster version of data frames

library(ggplot2) # for plotting
library(gridExtra) # for organizing plots

library(parallel) # for parallel processing

library(lcmsMetab) # LCMS data processing tools

# Timing --------------------------------------------------------------
startTime <- Sys.time()
startTimeStamp <- format(startTime, "(%Y%m%dh%H%M)")

# set or create the output folders -------------------------------------------

OutputFolderName <- "peakWindowMatricesAndPlots"

if (restart == TRUE) {

  outputFolderPath <- existingOutputPath

} else {

  outputFolderPath <- CreateOutputDir(outPutParentDir, OutputFolderName)

  outputMatrixDir <- CreateOutputDir(outputFolderPath, "Matrices")
  outputPlotDir <- CreateOutputDir(outputFolderPath, "Plots")

}

# Initialize summary file -------------------------------------------

logFilePath <- paste0(outputFolderPath,"/",
                      OutputFolderName,"LogFile",
                      startTimeStamp, ".txt")
logText <- paste0("start time ", startTimeStamp)
logText <- UpdateLogText(logText,"Initialize Summary File")
UpdateLogFile(logFilePath, logText)

# import window table and max intensity table ---------------------------------

# create vector of column numbers to load
vect.loadCols <- c(colNumWindowID,
                   colNumWindMzMin,
                   colNumWindMzMax,
                   colNumWindRtMin,
                   colNumWindRtMax)

dt.windows <- fread(paste0(outPutParentDir,"\\",windowCSV),
                    select = vect.loadCols)

colnames(dt.windows)[1] <- "windowID"
colnames(dt.windows)[2] <- "windMzMin"
colnames(dt.windows)[3] <- "windMzMax"
colnames(dt.windows)[4] <- "windRtMin"
colnames(dt.windows)[5] <- "windRtMax"


# import max intensity table first row to get column number
dt.maxIntensity <- fread(paste0(maxIntensityDir,"\\",csvMaxIntensity),
                      nrows = 1)

# column number to import
importCols <- ncol(dt.maxIntensity)

dropCols = NULL

# if window midpoints are included in the table, don't import them
if (colnames(dt.maxIntensity)[2] == "windMzMid") {

  # column number to import - don't impor windMzMid and WindRtMid
  importCols <- ncol(dt.maxIntensity) - 2

  dropCols = c(2,3)

}

# import table assign classes
dt.maxIntensity <- fread(paste0(maxIntensityDir,"\\",csvMaxIntensity),
                         drop = dropCols,
                         colClasses = list(double = 2:importCols),
                         integer64 = "double")

dt.maxIntensity <- NaZeroEmptyRmDup(dt.maxIntensity)

# for restart, remove completed windows ------------------------------------

# if restarting, remove windows that already have a matrix

# list windows that already have a plot
dt.plottedWindows <- data.table(windowID = list.files(path = outputMatrixDir,
                                           pattern = ".csv",
                                           full.names = FALSE,
                                           ignore.case = TRUE))
# remove extension
dt.plottedWindows[, windowID := gsub("\\.csv","",windowID)]

# remove windows with plots from dt.windows and dt.windowDirs
dt.windows <- dt.windows[!dt.plottedWindows, on = .(windowID)]

rm(dt.plottedWindows)

# sequence of integers for each file
windowSeq <- 1:nrow(dt.windows)

print(paste0("Windows to plot: ", length(windowSeq)))

# determine max lines to plot -----------------------------------------------

maxLinesHigh <- ceiling(maxSignalFiles*propHighPeaks)
maxLinesOther <- maxSignalFiles - maxLinesHigh

# defined functions ---------------------------------------------------------

# function for creating a table for a single peak window and context
# from a single MS1 table
GenWindowTable <- function(fileName, fndt.windows, ms1TableDir) {

  # pull values for subsetting
  lp.mzMin <- fndt.windows$windMzMin[1]
  lp.mzMax <- fndt.windows$windMzMax[1]
  lp.windRtMin <- fndt.windows$windRtMin[1]
  lp.windRtMax <- fndt.windows$windRtMax[1]
  lp.rtMin <- fndt.windows$rtStart[1]
  lp.rtMax <- fndt.windows$rtEnd[1]

  # read in the ms1 table
  dt.ms1Table <- fread(paste0(ms1TableDir, "/", fileName,".csv"),
                       colClasses = list(double = 1:3))

  # replace NA with 0
  dt.ms1Table[is.na(dt.ms1Table)] <- 0

  # subset table to within context scan window
  lpdt.ms1Table.context <- dt.ms1Table[scanRT %between% c(lp.rtMin,lp.rtMax)]

  if (nrow(lpdt.ms1Table.context) > 0) {

    # create a table of unique retention times in that range
    dt.contextWindScans <- unique(lpdt.ms1Table.context, by = "scanRT")[, list(scanRT)]

    # subset to mz range
    lpdt.ms1Table.context <-
      lpdt.ms1Table.context[mz %between% c(lp.mzMin,lp.mzMax)]

    # subset to RT and mz of window
    lpdt.ms1Table.window <-
      lpdt.ms1Table.context[scanRT %between% c(lp.windRtMin,lp.windRtMax)]

    # check if there are rows in the window
    if (nrow(lpdt.ms1Table.window) > 0) {

      # get maximum intensity for each scan
      lpdt.scanMaxInt <- lpdt.ms1Table.context[, max(intensity),
                                                     keyby = scanRT]

      # add max scan window intensities to scan table
      setkey(dt.contextWindScans, scanRT)
      dt.contextWindScans[lpdt.scanMaxInt, intensity := V1]

      # replace NA with 0 for scans with nothing in mz window
      dt.contextWindScans[is.na(dt.contextWindScans)] <- 0


    } else {
      # create empty scan table
      dt.contextWindScans <- data.table(scanRT = numeric(length = 0),
                               intensity = numeric(length = 0))

    } # end rows in the window if-else statement
  } else {
    # if no rows in context RT
    # create empty scan table
    dt.contextWindScans <- data.table(scanRT = numeric(length = 0),
                             intensity = numeric(length = 0))
  } # end rows in context RT window if statement

  return(dt.contextWindScans)

} # end WriteWindowTable function

# function to generate row plots
PlotPeakWindow <- function(list.peakWindowTables, outputPlotDir,
                           fndt.windows, maxPlotIntensity) {

  windowName <- fndt.windows$windowID[1]
  windMzMin <- fndt.windows$windMzMin[1]
  windMzMax <- fndt.windows$windMzMax[1]
  windRtMin <- fndt.windows$windRtMin[1]
  windRtMax <- fndt.windows$windRtMax[1]

  # combine the data tables
  dt.plotData <- rbindlist(list.peakWindowTables, use.names = TRUE,
                        fill = FALSE, idcol = TRUE)

  # make id column a factor
  dt.plotData[, .id := as.factor(.id)]

  # set NA intensity to zero
  dt.plotData[is.na(intensity), intensity := 0]

  # get max value in window
  maxWindowVal <- maxPlotIntensity

  if (nrow(dt.plotData) > 0) {

    # if there are scans in the window, generate plot

    titleText <- paste0("window: ", windowName, "\n",
                        "m/z: ", windMzMin, " - ", windMzMax, "\n",
                        "RT: ", windRtMin, " - ", windRtMax)

    plotTheme <-
      theme(title = element_text(size = 10),
            axis.title = element_text(size = 9),
            axis.text.x = element_text(size = 8),
            axis.text.y = element_text(size = 8),
            panel.background = element_blank(),
            panel.border = element_rect(fill = NA, linetype = 1, size = .5),
            panel.grid.major.y = element_line(colour = "grey50", size = 0.1),
            panel.grid.major.x = element_blank(),
            plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
            panel.spacing = unit(0, "pt"),
            strip.switch.pad.grid = unit(0, "pt"),
            strip.switch.pad.wrap = unit(0, "pt"))

    plot <-
      ggplot() +
      geom_line(aes(x = scanRT, y = intensity, color = .id), data = dt.plotData,
                show.legend = FALSE,
                size = 0.3) +
      scale_y_continuous(expand = c(0,0),
                         labels = function(x) formatC(x, digits = 1,
                                                      format = "E")) +
      scale_x_continuous(expand = c(0,0),
                         breaks = function(x) seq(from = round(2*min(x),1)/2,
                                                  to = max(x), by = 0.05)) +
      geom_vline(xintercept = windRtMin, linetype = "longdash", size = .6) +
      geom_vline(xintercept = windRtMax, linetype = "longdash", size = .6) +
      coord_cartesian(ylim = c(0, maxWindowVal*1.05)) +
      ggtitle(titleText) +
      xlab("Retention Time") +
      plotTheme

    ggsave(filename = paste0(windowName,".jpg"),
           plot = plot,
           device = "jpeg",
           path = outputPlotDir,
           width = 9,
           height = 4.5,
           units = "in",
           dpi = 400)

  } else {

    # no signal in window generate text file
    textFile <- file(paste0(outputPlotDir,"/",windowName,".txt"))
    textMessage <- "no signal in this window from specified tables"
    writeLines(textMessage, textFile )
    close(textFile)

  }

}

# function to get interpolated points from each file
InterpToSlices <- function(dt.ms1.window, fndt.windows,
                           vect.RtSlicePoints) {

  if (nrow(dt.ms1.window) > 2) {

    # get interpolated vector
    vect.interpIntensity <- approx(x = dt.ms1.window$scanRT,
                                   y = dt.ms1.window$intensity,
                                   xout = vect.RtSlicePoints,
                                   method = "linear",
                                   n = NULL, rule = 2, ties = mean)

    dt.interpIntensity <- data.table(x = vect.interpIntensity$x,
                                     y = vect.interpIntensity$y)

    # replace negative time values with 0
    dt.interpIntensity[, interpIntensity := ifelse(x <= 0, 0, y)]

    vect.filePeak <- dt.interpIntensity$interpIntensity

  } else {

    # if too few rows to intperpolate create zero vector
    vect.filePeak <- rep(0, length.out = length(vect.RtSlicePoints))

  }

  return(vect.filePeak)

}

# function to create peak group image matrices
PeakGroupMatrix <- function(windNum, dt.windows, dt.maxIntensity, ms1TableDir,
                            outputFolderPath, outputMatrixDir,
                            outputPlotDir, maxSignalFiles,
                           maxLinesHigh, maxLinesOther, RtDur, minRtSep) {

  # get dt.windows row corresponding to windNum
  fndt.windows <- dt.windows[windNum]

  # get window boundaries
  fn.windID <- fndt.windows$windowID[1]
  windRtMin <- fndt.windows$windRtMin[1]
  windRtMax <- fndt.windows$windRtMax[1]
  windMzMin <- fndt.windows$windMzMin[1]
  windMzMax <- fndt.windows$windMzMax[1]

  # print windowName to add to log
  print(fn.windID)

  # get corresponidng row of max intensity table
  dt.maxIntensity.window <- dt.maxIntensity[windowID == fn.windID]

  # melt to long table of files and intensities
  dt.maxIntensity.window <- melt(dt.maxIntensity.window,id.vars = "windowID",
                   measure.vars = c(2:(ncol(dt.maxIntensity.window))),
                   variable.name = "fileName",
                   value.name = "intensity",
                   na.rm = TRUE,
                   variable.factor = FALSE)

  # remove unnecessary windowID column
  dt.maxIntensity.window$windowID <- NULL

  # find RT mid point
  windRtMid <- (windRtMax + windRtMin) / 2

  # number of slices
  numRtSlices <- floor(RtDur / minRtSep)

  # find RT start point for matrix and end point
  rtSliceStart <- windRtMid - RtDur/2

  rtSliceEnd <- windRtMid + RtDur/2

  # add to table
  fndt.windows[, rtStart := rtSliceStart]
  fndt.windows[, rtEnd := rtSliceEnd]

  # create RT slice point vector
  vect.RtSlicePoints <- round(seq(from = rtSliceStart, to = rtSliceEnd,
                            length.out =  numRtSlices),4)

  # remove files with max signal of 0 in window
  dt.maxIntensity.window <- dt.maxIntensity.window[intensity > 0]

  # check that some files have signal in the window
  if (nrow(dt.maxIntensity.window) > 0) {

    # get max value in window
    maxWindowVal <- max(dt.maxIntensity.window$intensity,
                        na.rm = TRUE)

    # if more files remain than maximum plot lines select based on input
    if (nrow(dt.maxIntensity.window) > maxLinesHigh + maxLinesOther) {
      # sort by intensity
      setorder(dt.maxIntensity.window, -intensity)

      # create list based on inputs for number of lines to plot
      fileListPlot <- c(dt.maxIntensity.window$fileName[1:maxLinesHigh],
                    sample(dt.maxIntensity.window$fileName[maxLinesHigh + 1:maxLinesOther]))

    } else {
      # plot all the files with signal in window
      fileListPlot <- dt.maxIntensity.window$fileName
    }

    # create list of peak window tables
    list.peakWindowTables <- lapply(fileListPlot, GenWindowTable,
                                    fndt.windows, ms1TableDir)

    peakWindSeq <- 1:length(list.peakWindowTables)

    # feed list to plot function
    # get max plot y-value
    maxPlotIntensity <- max(dt.maxIntensity.window$intensity, na.rm = TRUE)
    # apply function
    PlotPeakWindow(list.peakWindowTables, outputPlotDir,
                   fndt.windows, maxPlotIntensity)

    # feed list create matrix
    dataset <- lapply(list.peakWindowTables, InterpToSlices,
                      fndt.windows, vect.RtSlicePoints)

    dataset <- do.call(rbind, dataset)

    # if not enough files with signal to complete the matrix,
    # fill the remaining rows with 0s
    if (nrow(dataset) < maxSignalFiles) {

      # number of rows to fill
      numFillRows <- maxSignalFiles - nrow(dataset)

      # create 0 matrix for remaining rows
      mat.fill <- matrix(0L, nrow = numFillRows, ncol = numRtSlices)

      dataset <- rbind(dataset, mat.fill)

    }

    # get matrix maximum
    maxIntMatrix <- max(dataset)

    # divide values by maximum
    dataset <- dataset / maxIntMatrix

    # add vector for with 1s for window

    # create table
    dt.binaryWindow <- data.table(rt = vect.RtSlicePoints)

    # set intensities within window to 1
    dt.binaryWindow[, intensity :=
                      ifelse(rt >= windRtMin & rt <= windRtMax, 1, 0)]

    # add to matrix
    dataset <- rbind(dataset, matrix(data = dt.binaryWindow$intensity, nrow = 1))

    # visual check
    # image(dataset, useRaster = TRUE, axes = FALSE,
    # col = grey(seq(0, 1, length = 256)))

    # set 0 to NA to reduce file size
    dataset[dataset == 0] <- NA

    # write the matrix
    writePath <- paste0(outputMatrixDir, "/", fn.windID, ".csv")
    fwrite(dataset, file = writePath, row.names = FALSE, col.names = FALSE)



  } else {

    # no signal in window generate text file
    textFile <- file(paste0(outputMatrixDir,"/",fn.windID,".txt"))
    textMessage <- "no signal in this window from specified tables"
    writeLines(textMessage, textFile )
    close(textFile)

  }

}


# setup cluster ------------------------------------------

# create path for cluster outfile
clusterLogPath <- paste0(outputFolderPath,"/","clusterLog",
                         startTimeStamp,".txt")

# use 75% of cores
numCores <- max(ceiling(detectCores() * propProcessors),2)

# make a cluster
clustr <- makePSOCKcluster(names = numCores,
                      outfile = clusterLogPath)

# load packages into cluster
clusterEvalQ(clustr, {
  library(data.table)
  library(ggplot2)
  library(lcmsMetab)
})

# pass needed objects to cluster
clusterExport(clustr, c("PeakGroupMatrix",
                        "InterpToSlices",
                        "PlotPeakWindow",
                        "GenWindowTable",
                        "dt.windows", "dt.maxIntensity", "ms1TableDir",
                        "outputFolderPath", "outputMatrixDir",
                        "outputPlotDir", "maxSignalFiles",
                        "maxLinesHigh", "maxLinesOther",
                        "RtDur", "minRtSep"))

logText <-
  UpdateLogText(logText,"cluster setup",runTime(startTime))
UpdateLogFile(logFilePath, logText)

# write tables - parallel processing - use functions ---------------------------------------

# non parallel application for development
# lapply(windowSeq, PeakGroupMatrix,
#           dt.windows, dt.maxIntensity, ms1TableDir,
#           outputFolderPath, outputMatrixDir,
#           outputPlotDir, maxSignalFiles,
#           maxLinesHigh, maxLinesOther, RtDur, minRtSep)

# apply in parallel with cluster
parLapply(clustr, windowSeq, PeakGroupMatrix,
          dt.windows, dt.maxIntensity, ms1TableDir,
          outputFolderPath, outputMatrixDir,
          outputPlotDir, maxSignalFiles,
          maxLinesHigh, maxLinesOther, RtDur, minRtSep)

logText <-
  UpdateLogText(logText, "peak windows processed", runTime(startTime))
UpdateLogFile(logFilePath, logText)

# close the cluster ---------------------------------------------

stopCluster(clustr)
rm(clustr)

# Endscript (record keeping) -----------------------------------------

# get script file name
scriptName <- basename(sys.frame(1)$ofile)

# remove the file extension ".R"
scriptName <- gsub("\\.R","",scriptName)

logText <-
  UpdateLogText(logText,"end script",runTime(startTime))
UpdateLogFile(logFilePath, logText)

if (1 == 1) {
  print("script complete")
  print(runTime(startTime))
}

