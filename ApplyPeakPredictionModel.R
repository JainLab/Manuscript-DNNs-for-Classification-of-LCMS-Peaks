# Description -----------------------------------------------------------------

# prepare image matrices for model. Use model for prediction

# User Input ------------------------------------------------------------------

# Path to folder with individual peak table directories
# for Windows operating system use \\ or / instead of \ to separate directories
# example:
# peakTableParentDir <- 'M:\\folder1\\folderWithFile'
peakMatrixParentDir <- "E:\\Example1\\Example"

# path to folder where output directory will be created
outPutParentDir <- "E:\\Example2\\Example"

# model file name
# model must be located in the parent folder
modelFile <- "modelFileName.h5"

# peak window plot directory
peakWindowPlotDir <- "E:\\Example3\\Example"


# End user input ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

# advanced user input (probably don't need to change) ------------------------

# create peak review table? This is used for reviewing the peak window plots
# using the generated jpegs and the reivew app
genPlotReviewTable <- TRUE

# restart not fully implmented TODO
restart <- FALSE
existingOutputPath <- NA

# Library (packages)--------------------------------------------------------
library(data.table) # faster version of data frames

library(ggplot2) # for plotting
library(gridExtra) # for organizing plots

library(keras) # for neural networks

library(abind) # for binding multi-dimension arrays

library(parallel) # for parallel processing

library(lcmsMetab) # LCMS data processing tools

# Timing --------------------------------------------------------------
startTime <- Sys.time()
startTimeStamp <- format(startTime, "(%Y%m%dh%H%M)")

# set or create the output folder -------------------------------------------

OutputFolderName <- "ApplyPeakPredictionModel"

if (restart == TRUE) {

  outputFolderPath <- existingOutputPath

} else {

  # create output folder if it doesn't exist
  outputFolderPath <- paste0(outPutParentDir,"/",
                             OutputFolderName,startTimeStamp)

  if (!dir.exists(outputFolderPath)) {
    dir.create(outputFolderPath)
  }
}

# Initialize summary file -------------------------------------------

logFilePath <- paste0(outputFolderPath,"/",
                      OutputFolderName,"LogFile",
                      startTimeStamp, ".txt")
logText <- paste0("start time ", startTimeStamp)
logText <- UpdateLogText(logText,"Initialize Summary File")
UpdateLogFile(logFilePath, logText)

# determine files to load ---------------------------------------------------

# create a list of peak window tables
fileListPath <- list.files(peakMatrixParentDir, full.names = TRUE,
                           pattern = ".csv")

# make table
dt.peakFiles <- data.table(path = fileListPath)

# create column with just window names
dt.peakFiles[, windName := gsub("^.*/", "", path)]
dt.peakFiles[, windName := gsub("\\.csv", "", windName)]

# import function ----------------------------------------------------------

# function for reading and converting to matrix
ReadToMatrix <- function(path) {

  # load
  data <- fread(file = path, header = FALSE)

  # replace zeros removed in export
  data[is.na(data)] <- 0

  # convert to matrix
  data <- as.matrix(data)

}

# import and re-format data ---------------------------------------------------

# read in tables to list
array.data <- lapply(dt.peakFiles$path, ReadToMatrix)

# convert to array
array.data <- abind(array.data, along = 0, force.array = TRUE)

# import the model ----------------------------------------------------------

model <- load_model_hdf5(filepath = paste0(outPutParentDir,"/",modelFile))

# apply model ----------------------------------------------------------------

linesPerImage <- dim(array.data)[2]
rtBins <- dim(array.data)[3]

# reshape the data
array.data <- array_reshape(array.data, c(nrow(array.data),
                                          linesPerImage, rtBins, 1))

# create data table of predicted probabilities
dt.prediction  <- predict_proba(model, array.data)
dt.prediction <- data.table(dt.prediction)
colnames(dt.prediction)[1] <- "predProb"

# write prediction table -----------------------------------------------------

dt.windowPred <- data.table(windName = dt.peakFiles$windName,
                            predProb = dt.prediction$predProb)

tableName <- paste0("windowPredictionTable", startTimeStamp)

WrtTable(dt.windowPred, outputFolderPath, tableName)

# write peak review format table ---------------------------------------------

if (genPlotReviewTable == TRUE) {

  # create a list of peak window plots
  fileListPath <- list.files(peakWindowPlotDir, full.names = TRUE,
                             pattern = ".jpg")

  # make table
  dt.plotFiles <- data.table(path = fileListPath)

  # create column with just window names
  dt.plotFiles[, windName := gsub("^.*/", "", path)]
  dt.plotFiles[, windName := gsub("\\.jpg", "", windName)]

  # add path to plot image to table
  setkey(dt.plotFiles, windName)
  setkey(dt.windowPred, windName)
  dt.windowPred[dt.plotFiles, path := path]

  # add columns for peak review app table
  dt.windowPred[, rowNum := 1:nrow(dt.windowPred)]
  dt.windowPred[, viewed := rep.int(0, times = nrow(dt.windowPred))]
  dt.windowPred[, good1bad0 := rep.int(1, times = nrow(dt.windowPred))]
  dt.windowPred[, windowBad := rep.int(0, times = nrow(dt.windowPred))]
  dt.windowPred[, borderline := rep.int(0, times = nrow(dt.windowPred))]
  dt.windowPred[, image := paste0(windName,".jpg")]

  # reorder columns
  setcolorder(dt.windowPred, c("rowNum",
                               "image",
                               "viewed",
                               "good1bad0",
                               "windowBad",
                               "borderline",
                               "path",
                               "predProb"))

  dt.windowPred$windName <- NULL

  # use 0.5 threshold to as boundary between good and bad
  dt.windowPred[, good1bad0 := ifelse(predProb >= 0.5,1,0)]

  # reset row numbering - shuffle
  dt.windowPred[, rowNum := sample(nrow(dt.windowPred))]

  # sort by the new row order
  setorder(dt.windowPred, rowNum)

  tableName <- paste0("windowReviewAppTable", startTimeStamp)
  WrtTable(dt.windowPred, outputFolderPath, tableName)

} # end if statement for peak review table generation

# plot histogram -------------------------------------------------------------

# plot title
titleText <- paste0("Prediction Score Histogram")

# plot object
plotHist <-
  ggplot(dt.windowPred, aes(predProb)) +
  geom_histogram(binwidth = 0.025, center = 0.0125) +
  ggtitle(titleText) +
  ylab('training peak count') +
  xlab('Probability Good') +
  ggThemeLcmsMetab()

# save plot
ggsave("plotHist.pdf", plot = plotHist, device = "pdf",
       path = outputFolderPath,
       width = 11,
       height = 8.5,
       units = "in")

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

