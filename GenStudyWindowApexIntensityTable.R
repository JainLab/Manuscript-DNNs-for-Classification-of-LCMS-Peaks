# Description ------------------------------------

# Create a table that includes the max intensity for each MS1 table for
# each window

# clear existing variables-----------------------------------------------------
# remove variables from workspace
rm(list = ls())
gc()

# User Input ------------------------------------------------------------------

# Path to folder with MS1 table csv files
# use \\ instead of \ to separate directories
# example:
# mzXMLDir <- 'M:\\folder1\\foldwerWithFile'
ms1TableDir <- "M:\\folder1\\foldwerWithFile"

# path to folder where output directory will be created
outPutParentDir <- "M:\\folder2\\foldwerWithFile"

# csv filename with peak windows in parent directory
windowCSV <- "filename.csv"

# proportion of processors to use 0-1
propProcessors <- 0.2

# include window midpoints?
windowMidpoints <- TRUE

# End user input ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

# advanced user input (probably don't need to change) ------------------------

# column numbers - change if not using standard window table
colNumWindowID <- 1
colNumWindMzMin <- 2
colNumWindMzMax <- 3
colNumWindRtMin <- 4
colNumWindRtMax <- 5

# Library (packages)--------------------------------------------------------
library(data.table) # faster version of data frames

library(parallel) # for parallel processing

library(lcmsMetab) # LCMS data processing tools

# Timing --------------------------------------------------------------
startTime <- Sys.time()
startTimeStamp <- format(startTime, "(%Y%m%dh%H%M)")

# create the output folder -------------------------------------------

OutputFolderName <- paste0("WindowMaxIntensityTable",startTimeStamp)

# create output folder if it doesn't exist
outputFolderPath <- paste0(outPutParentDir,"/",OutputFolderName)

if (!dir.exists(outputFolderPath)) {
  dir.create(outputFolderPath)
}

# Initialize summary file -------------------------------------------

logFilePath <- paste0(outputFolderPath,"/",
                      OutputFolderName,"LogFile",
                      startTimeStamp, ".txt")
logText <- paste0("start time ", startTimeStamp)
logText <- UpdateLogText(logText,"Initialize Summary File")
UpdateLogFile(logFilePath, logText)

# import window table -------------------------------------

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

# List MS1 files ---------------------------------------

filePathList <- list.files(path = ms1TableDir, pattern = "\\.csv$",
                           full.names = TRUE, recursive = FALSE,
                           ignore.case = TRUE)

fileNameList <- list.files(path = ms1TableDir, pattern = "\\.csv$",
                           full.names = FALSE, recursive = FALSE,
                           ignore.case = TRUE)

# sequence of integers for each file
fileSeq <- 1:length(filePathList)

# create table
dt.ms1TableFiles <- data.table(path = filePathList,
                               name = fileNameList)

# functions: create window max table for single ms1 table -------------------

GetWindowMaxPoint <- function(windowNum, dt.ms1Table, tableName,
                              dt.windows) {

  # for a single window and single dt.ms1Table extract the max value in the
  # window

  # pull values for subsetting
  windowID <- dt.windows$windowID[windowNum]
  windMzMin <- dt.windows$windMzMin[windowNum]
  windMzMax <- dt.windows$windMzMax[windowNum]
  windRtMin <- dt.windows$windRtMin[windowNum]
  windRtMax <- dt.windows$windRtMax[windowNum]

  # subset table to rt range
  dt.ms1Table.window <-
    dt.ms1Table[scanRT %between% c(windRtMin,windRtMax)]
  # subset to mz range
  dt.ms1Table.window <-
    dt.ms1Table.window[mz %between% c(windMzMin,windMzMax)]

  if (nrow(dt.ms1Table.window) > 0) {

    # if there is at least one row in the window take the max intensity
    maxIntensity <- max(dt.ms1Table.window$intensity, na.rm = TRUE)

  } else {

    maxIntensity <- NA

  }

  # create single row data table
  dt.windowMaxIntensity <- data.table(file = tableName,
                                      windowID = windowID,
                                      maxIntensity = as.numeric(maxIntensity))

  return(dt.windowMaxIntensity)

}

WindowMaxPointTable <- function(fileNum, dt.ms1TableFiles, dt.windows) {

  # function pulls the highest intensity value for each window in a file
  # and creates a long form table

  # dt.windows is the data table of the peak windows ID, MzMin, MzMax, RtMin,
  # RtMax, extractRtMin, extractRtMax
  # dt.ms1TableFiles is the table of ms1 tables with the path in column 1
  # and the name in column 2
  # fileNum is the row number corresponding to the MS1 file in dt.ms1TableFiles
  # outputFolderPath is the folder where the individual peak directories are
  # created

  # print for progress reporting
  print(dt.ms1TableFiles$name[fileNum])

  # Read the MS1 table
  dt.ms1Table <- fread(dt.ms1TableFiles$path[fileNum])

  # get the table name (from the original LCMS file)
  tableName <- gsub(".csv","",dt.ms1TableFiles$name[fileNum])

  # create a sequence for iterating the windows
  windowSeq <- 1:nrow(dt.windows)

  # lapply function to get list of single row window max tables
  dt.FileWindowMax <- lapply(windowSeq, GetWindowMaxPoint,
                             dt.ms1Table, tableName, dt.windows)

  # combine list into single table
  dt.FileWindowMax <- rbindlist(dt.FileWindowMax, use.names = FALSE,
                                fill = FALSE, idcol = NULL)

  return(dt.FileWindowMax)


}

# setup cluster ------------------------------------------

# create path for cluster outfile
clusterLogPath <- paste0(outputFolderPath,"/","clusterLog",".txt")

# use 75% of cores
numCores <- max(floor(detectCores() * propProcessors),2)

# make a cluster
clustr <- makePSOCKcluster(names = numCores,
                           outfile = clusterLogPath)

# load packages into cluster
clusterEvalQ(clustr, {
  library(data.table)
  library(lcmsMetab)
})

# pass needed objects to cluster
clusterExport(clustr, c("WindowMaxPointTable",
                        "GetWindowMaxPoint",
                        "dt.windows",
                        "dt.ms1TableFiles"))

logText <-
  UpdateLogText(logText,"cluster setup",runTime(startTime))
UpdateLogFile(logFilePath, logText)

# write tables - parallel processing - use functions ---------------------------------------

# apply in parallel with cluster
dt.fileWindowTable <- parLapply(clustr, fileSeq, WindowMaxPointTable,
          dt.ms1TableFiles, dt.windows)

logText <-
  UpdateLogText(logText,"peak max table complete",runTime(startTime))
UpdateLogFile(logFilePath, logText)

# close the cluster ---------------------------------------------

stopCluster(clustr)
rm(clustr)

# combine list -------------------------------------------------------------

dt.fileWindowTable <- rbindlist(dt.fileWindowTable)

# create wide-format export table -----------------------------------------

dt.fileWindowIntWide <- dcast(dt.fileWindowTable, windowID ~ file)

# add window midpoint values -----------------------------------------------------

if (windowMidpoints == TRUE) {

  # create columns
  dt.windows[, windMzMid := round((windMzMin + windMzMax)/2, 4)]
  dt.windows[, windRtMid := round((windRtMin + windRtMax)/2,3)]

  # add columns to height table
  setkey(dt.windows, windowID)
  setkey(dt.fileWindowIntWide, windowID)
  dt.fileWindowIntWide[dt.windows, windMzMid := windMzMid]
  dt.fileWindowIntWide[dt.windows, windRtMid := windRtMid]

  # reorder
  setcolorder(dt.fileWindowIntWide, c("windowID", "windMzMid", "windRtMid"))

}

# write table --------------------------------------------------------------

WrtTable(dt.fileWindowIntWide, outputFolderPath,
         paste0("WindowMaxIntensityTable", startTimeStamp))

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