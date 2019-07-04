# Description ------------------------------------

# writes MS1 tables for each mzXML in a given directory

# User Input ------------------------------------------------------------------

# Path to folder with mzXML files
# use \\ instead of \ to separate directories
# example:
# mzXMLDir <- 'M:\\folder1\\foldwerWithFile'
mzXMLDir <- "M:\\folder1\\foldwerWithFile"

# path to folder where output directory will be created
outPutParentDir <- "M:\\folder2\\foldwerWithFiles"

# End user input ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

# Library (packages)--------------------------------------------------------
library(data.table) # faster version of data frames
library(mzR) # for reading mzXML files
library(parallel) # for parallel processing

library(lcmsMetab) # LCMS data processing tools

# Timing --------------------------------------------------------------
startTime <- Sys.time()
startTimeStamp <- format(startTime, "(%Y.%m.%d.%H%M)")

# create the output folder -------------------------------------------

OutputFolderName <- "lcmsTables"

outputFolderPath <- CreateOutputDir(outPutParentDir, OutputFolderName)

# Initialize summary file -------------------------------------------

logFilePath <- paste0(outputFolderPath,"/",
                      OutputFolderName,"LogFile",
                      startTimeStamp, ".txt")
logText <- paste0("start time ", startTimeStamp)
logText <- UpdateLogText(logText,"Initialize Summary File")
UpdateLogFile(logFilePath, logText)

# defined functions ---------------------------------------------------------

# given a scan export the mz-intensity table and add scan RT
scanExtract <- function(iScan, mzRobj, dt.scanRT) {

  # iScan is the seqNum of the scan
  # mzRobj is an MzR object
  # dt.header is a data.table of scan header information from
  # mzR header function

  # extract the table for the first relevant scan
  dt.rtMzIntensity <- data.table(peaks(mzRobj, iScan))
  colnames(dt.rtMzIntensity)[1] <- "mz"
  colnames(dt.rtMzIntensity)[2] <- "intensity"

  # get retention time for relevant scan
  dt.scanRT.select <- subset(dt.scanRT, seqNum == iScan)

  scanRT <- dt.scanRT.select$retentionTime[1]

  # add retention time as column scanRT
  dt.rtMzIntensity[, scanRT := scanRT]

  return(dt.rtMzIntensity)

}

# write an MS1 table for an mzXML file
writeMS1Table <- function(fileNum, dt.mzXMLFiles, outputFolderPath) {

  # mzXMLpath is the path to the mzXML file
  # outputFolderPath is the path the output folder
  # fileNum is the number that corresponds to the row in dt.mzXMLFiles
  # for the file to extract from
  # dt.window is the table of windows, passed to featExtract function
  # dt.mzXMLFiles is a 2 column table of the mzXML files to extract
  # from. First column is path, second is name


  # ifile is the file path of an individual mzXML
  mzXMLpath <- dt.mzXMLFiles$path[fileNum]

  # pull file name
  mzXMLname <- dt.mzXMLFiles$name[fileNum]

  # open the file
  mzRobj <- openMSfile(mzXMLpath)

  # create data table of header info for all scans
  dt.header <- data.table(header(mzRobj))

  # MS1 scan vector
  vect.MS1scanNums <- subset(dt.header,
                             msLevel == 1,
                             select = c("seqNum"))
  vect.MS1scanNums <- vect.MS1scanNums$seqNum

  dt.scanRT <- subset(dt.header,
                      select = c("seqNum","retentionTime"))

  # use defined function scanExtract to
  # create a list of rt-Mz-Intensities tables
  dt.rtMzIntensityMS1 <- lapply(vect.MS1scanNums, scanExtract,
                                mzRobj = mzRobj,
                                dt.scanRT = dt.scanRT)

  # combine the list of tables into a single table
  dt.rtMzIntensityMS1 <- rbindlist(dt.rtMzIntensityMS1, use.names = TRUE)

  # reduce precision to legitimate levels and convert RT to minutes
  dt.rtMzIntensityMS1[,mz := round(mz,4)]
  dt.rtMzIntensityMS1[,intensity := round(intensity,0)]
  dt.rtMzIntensityMS1[,scanRT := round(scanRT/60,3)]

  # write table
  WrtTable(dt.rtMzIntensityMS1,outputFolderPath,gsub(".mzXML","",mzXMLname))
}

# List mzXML files ---------------------------------------

filePathList <- list.files(path = mzXMLDir, pattern = "\\.mzXML$",
                           full.names = TRUE, recursive = FALSE,
                           ignore.case = TRUE)

fileNameList <- list.files(path = mzXMLDir, pattern = "\\.mzXML$",
                           full.names = FALSE, recursive = FALSE,
                           ignore.case = TRUE)

# sequence of integers for each file
fileSeq <- 1:length(filePathList)

# create table
dt.mzXMLFiles <- data.table(path = filePathList,
                            name = fileNameList)

# setup cluster ------------------------------------------

# use 75% of cores
numCores <- ceiling(detectCores() * .75)

# make a cluster
clustr <- makeCluster(numCores)

# load packages into cluster
clusterEvalQ(clustr, {
  library(data.table)
  library(mzR)
  library(lcmsMetab)
})

# pass needed objects to cluster
clusterExport(clustr, c("outputFolderPath",
                        "filePathList",
                        "fileSeq",
                        "dt.mzXMLFiles",
                        "writeMS1Table",
                        "scanExtract"))

logText <-
  UpdateLogText(logText,"cluster setup",runTime(startTime))
UpdateLogFile(logFilePath, logText)

# write tables - parallel processing - use functions ---------------------------------------

# apply in parallel with cluster
dt.peakExtract <- parLapply(clustr, fileSeq, writeMS1Table,
                            dt.mzXMLFiles, outputFolderPath)

logText <-
  UpdateLogText(logText,"peak height extraction complete",runTime(startTime))
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






