# Description ------------------------------------

# uses tables from manual review and model peak scoring to remove likely
# bad peak windows from the window list

# clear existing variables-----------------------------------------------------
# remove variables from workspace
rm(list = ls())
gc()

# User Input ------------------------------------------------------------------

# for Windows operating system use \\ or / instead of \ to separate directories
# example:
# outPutParentDir <- 'M:\\folder1\\foldwerWithFile'

# path to folder where output directory will be created
# put copies of the CSV files in this folder
outPutParentDir <- "E:\\2019.06_MicrobiomeStrainMediaMultiTime_Lipidomics_Pos\\2019.06.18_Extraction\\8.FilterPeaks"

# CSV peak window table - contains window bounds - standard format
csvWindowTable <- "MZminePeakWindowViewerTable(20190620h0041).csv"

# use manual review table?
useManualReview <- TRUE

# CSV manual review tracking table - standard format
# if not used, set to NA
csvManualReviewTable <- "reviewTracking(2019.06.21.1059).csv"

# use model predictions score?
useModelScore <- TRUE

# CSV peak score table - standard format
# Variable not used if useModelScore set to FALSE
csvPredictionScoreTable <- "windowPredictionTable(20190626h1345).csv"

# peak model score threshold - scores at or above this number will be retained
# Variable not used if useModelScore set to FALSE
scoreThresh <- .4


# End user input ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

# advanced user input (probably don't need to change) ------------------------


# Library (packages)--------------------------------------------------------
library(data.table) # faster version of data frames
library(lcmsMetab) # LCMS data processing tools

# Timing --------------------------------------------------------------
startTime <- Sys.time()
startTimeStamp <- format(startTime, "(%Y%m%dh%H%M)")

# set or create the output folders -------------------------------------------

OutputFolderName <- paste0("peakWindowFiltering")

# create output folder if it doesn't exist
outputFolderPath <- paste0(outPutParentDir,"/",OutputFolderName,startTimeStamp)

if (!dir.exists(outputFolderPath)) {
  dir.create(outputFolderPath)
}

# Initialize summary file -------------------------------------------

logFilePath <- paste0(outputFolderPath,"/",
                      OutputFolderName,"_LogFile",
                      startTimeStamp, ".txt")
logText <- paste0("start time ", startTimeStamp)
logText <- UpdateLogText(logText,"Initialize Summary File")
UpdateLogFile(logFilePath, logText)

# import tables --------------------------------------------------------------

dt.windows <- fread(paste0(outPutParentDir,"/", csvWindowTable))

if (useManualReview == TRUE) {

  dt.manualReview <- fread(paste0(outPutParentDir,"/", csvManualReviewTable))

}

if (useModelScore == TRUE) {

  dt.predScore <- fread(paste0(outPutParentDir,"/", csvPredictionScoreTable))

}

# create combined table -------------------------------------------------------

setkey(dt.windows, windowID)

if (useManualReview == TRUE) {

  dt.manualReview[, windowID := gsub(".jpg$", "", image)]

  setkey(dt.manualReview, windowID)

  dt.windows[dt.manualReview, viewed := viewed]
  dt.windows[dt.manualReview, good1bad0 := good1bad0]

}

if (useModelScore == TRUE) {

  colnames(dt.predScore)[1] <- "windowID"
  setkey(dt.predScore, windowID)
  dt.windows[dt.predScore, predProb := predProb]

}

# filter ----------------------------------------------------------------------

if (useManualReview == TRUE & useModelScore == FALSE) {

  # remove viewed windows marked as bad
  # peaks are set as 1 (good) by default so zeros were manually labeled bad
  dt.windows <- dt.windows[good1bad0 == 1]

  # remove filtering columns
  dt.windows$viewed <- NULL
  dt.windows$good1bad0 <- NULL

}

if (useModelScore == TRUE & useManualReview == FALSE) {

  # remove unviewed windows with low score
  dt.windows <- dt.windows[predProb > scoreThresh]

  dt.windows$predProb <- NULL

}

if (useModelScore == TRUE & useManualReview == TRUE) {

  # remove unviewed windows with low score
  dt.windows <- dt.windows[(viewed == 1 & good1bad0 == 1) | (viewed == 0 & predProb > scoreThresh)]

  # remove filtering columns
  dt.windows$predProb <- NULL
  dt.windows$viewed <- NULL
  dt.windows$good1bad0 <- NULL

}

# write filtered table -------------------------------------------------------

# write table
WrtTable(dt.windows, outputFolderPath,
         paste0("filteredWindowTable",startTimeStamp))

# Endscript (record keeping) -----------------------------------------

# get script file name
scriptName <- basename(sys.frame(1)$ofile)

# remove the file extension ".R"
scriptName <- gsub("\\.R","",scriptName)

# copy script to folder to save as record
file.copy(sys.frame(1)$ofile,
          to = file.path(outputFolderPath,
                         paste0(scriptName, startTimeStamp, ".R")))

logText <-
  UpdateLogText(logText,"end script",runTime(startTime))
UpdateLogFile(logFilePath, logText)

if (1 == 1) {
  print("script complete")
  print(runTime(startTime))
}