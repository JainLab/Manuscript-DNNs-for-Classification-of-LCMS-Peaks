# Description ------------------------------------

# This script generates a peak window table from
# a peak list CSV file exported from MZmine as
# described in the User Input section

# User Input ------------------------------------------------------------------

# for Windows operating system use \\ instead of \ to separate directories
# example:

# Path to folder with input csv
# parentDir <- 'M:\\folder1\\folderWithFile'
parentDir <- "M:\\folder1\\folderWithFile"

# generate the CSV file in the folder from MZmine
# by creating a peak list of the peaks for which
# you want windows generated
# Peak List Methods > Export/Import > Export to CSV file
# Export common elements: check the following boxes
# Export row m/z
# Export row retention time
# Export data file elements: check the following boxes
# Peak m/z
# Peak RT
# Peak RT start
# Peak RT end
# Peak height
# Peak m/z min
# Peak m/z max
# Ensure that no other boxes are checked

# name of csv file
inputCSV <- "MZmine2Export.csv"


# set bounds on window size

# set initial m/z bounds on ppm?
# If set to TRUE, uses maxWindowPpm variable below to set the inital m/z bounds
# If set to FALSE, uses min(Peak m/z min) and max(Peak m/z max)
# Initial bounds may be altered if addressOverlap is set to TRUE below
mzBoundOnPpm <- FALSE

# maximum peak window ppm
# m/z window will not be wider than twice the value
# computed from this tolerance, and will only be
# narrower when necessary to split from overlapping
# window that could not be merged
# This should be as high as or higher than tolerances
# used in previous processing steps
maxWindowPpm <- 15

# minimum total peaks detected in window
# set to 0 to prevent subsetting by this metric
minTotalDetected <- 0

# minimum proportion of files with peaks detected in window
# value between 0 and 1
# set to 0 to prevent subsetting by this metric
minPropDetected <- 0

# consider observed peak tails when setting windows?
# When setting windows for viewing in MZmine, it is important to see more
# than just the apexes when determining if a peak is true. However, this
# increases the chance that parts of neighboring peaks will be included in
# windows. If set to TRUE, use rtQuantile input below to indicate what quantile
# of peak tails to use to set the window.
considerTails <- TRUE

# quantile of RT start and RT end to use to set initial retention time window
# value between 0 and 1. Larger leads to wider windows
# If considerTails is set to FALSE, this variable is not used
rtQuantile <- 1

##### Limit max RT span initial windows #
# maximum window RT duration in the same units as retention time
# RT window will not be wider than this value
# This should be as high as or higher than tolerances
# used in previous processing steps.
# enter NULL to avoid using this criteria
maxWindowRtDur <- 1
# maximum standard deviations from the median apex
# enter NULL to avoid using this criteria
maxWindowRtSd <- NULL

##### Limit min RT span initial windows #
# minimum window RT duration in the same units as retention time
# RT window will not be narrower than this value
# set greater than the maximum distance
# between two consecutive scans or window may be left
# empty.
minWindowRtDur <- 0.04
# minimum standard deviations from the median apex
# Set to 0 if you don't want this parameter used
# enter NULL to avoid using this criteria
minWindowRtSd <- 0

# Note that standard deviation thresholds are applied
# first and overwritten by absolute time thresholds
# where absolute time minimum is above the standard
# deviation maximum or time maximum is below the
# standard deviation minimum. Only the maximum absolute
# time threshold is used in the case of window merges.

# additional chromatogram to show before and after peak window
# exported as additional windows in a separate table
# units are the same as input retention time
rtContextWindowBuffer <- .5

##### address overlap
# address overlap?
# If set to TRUE, the script will address overlaps by splitting or combining
# windows such that no window overlaps with any other window.
# final windows may not meet criteria entered above.
# This option uses a clustering algorithm which can greatly increase the
# length of time required to complete the script and may be problematic
# with certain data sets.
addressOverlap <- FALSE

# End User Input Section ###########################################
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

# Advanced user inputs (probably don't need to change) ------------

# overlapping windows are resolved by clustering and creating new windows
# from clusters.
# What is the minimum number of apexes for a cluster to be considered valid?
minApexPerCluster <- 4

# plot windows?
# probably only useful for development or small test tables
plotWindows <- FALSE

# Timing --------------------------------------------------------------
startTime <- Sys.time()
startTimeStamp <- format(startTime, "(%Y%m%dh%H%M)")

# define function for current runtime
runTime <- function(){
  runTimeText <- format(round(Sys.time() -  startTime,3))
  return(runTimeText)
}

# timer <- Sys.time()
# print(format(as.numeric(Sys.time() - timer), digits = 10))

# create the output folder -------------------------------------------

OutputFolderName <- paste0("PeakWindowTable",startTimeStamp)

# create output folder if it doesn't exist
outputFolderPath <- paste0(parentDir,"/",OutputFolderName)

if (!dir.exists(outputFolderPath)) {
  dir.create(outputFolderPath)
}

# Initialize summary file -------------------------------------------

SummaryFileName <- paste0(outputFolderPath,"/",
                          "PeakWindowTableSummaryFile",
                          startTimeStamp, ".txt")
summaryText <- c("start time",startTimeStamp)

# define function for updating the summary file
updateSummary <- function(x,...){
  summaryText <- c(summaryText,"",x,..., runTime())
  SummaryFile <- file(SummaryFileName)
  writeLines(summaryText, SummaryFile)
  close(SummaryFile)
  return(summaryText)
}

summaryText <-
  updateSummary("Initialize Summary File")

# Library (packages)--------------------------------------------------------
library(data.table) # faster version of data frames
library(ggplot2) # for plotting - included in tidyverse
library(gridExtra) # for arranging multiple plot outputs

library(mclust)

# defined functions ----------------------------------------------------

# function for replacing the 0 and empty cells in imported
# data tables with NA and removing any blank columns
# also removes duplicated rows
repZeroEmpty <- function(x) {
  
  # input x is data table
  
  # change 0 cells to NA
  x[x == 0] <- NA
  # change empty cells to NA
  x[x == ""] <- NA
  
  # remove any columns that are all NA - resulting from blank column
  x <- x[,which(unlist(lapply(x, function(x)!all(is.na(x))))), with = FALSE]
  # x <- x[,which(unlist(lapply(x, function(x)!all(is.na(x)))))]
  
  # remove duplicate rows
  x <- unique(x)
  
  return(x)
  
}

# function for creating a featID column from the mz and RT
# of each feature and placing it before the first injection
# column
addFeatID <- function(x, mzColNum, rtColNum, firstInjColNum) {
  
  # x is the data table
  # mzColNum is the column number with the feature mz
  # rtColNum is the column number with the feature RT
  # firstInjColNum is the first column number that contains
  # injection measurements
  
  # Last label column
  LastLblColNum <- firstInjColNum - 1
  
  # save mz and RT column names
  nameMz <- colnames(x)[mzColNum]
  nameRt <- colnames(x)[rtColNum]
  
  # rename mz and RT columns
  colnames(x)[mzColNum] <- "mz"
  colnames(x)[rtColNum] <- "RT"
  
  # convert mz and RT columns to character vectors
  mzColVect <- sprintf("%09.4f",round(x[,mz],4))
  rtColVect <- sprintf("%4.3f",round(x[,RT],3))
  
  # create column in data set for feature ID combining mz and RT as text
  x[, featID :=  paste0("rt",rtColVect,"mz",mzColVect)]
  
  # revert names
  colnames(x)[mzColNum] <- nameMz
  colnames(x)[rtColNum] <- nameRt
  
  rm(mzColVect)
  rm(rtColVect)
  gc()
  
  # Put feature ID columns before injection
  newFirstInjColNum <- firstInjColNum + 1
  setcolorder(x, c((1:LastLblColNum),ncol(x), (newFirstInjColNum:ncol(x) - 1)))
  
  return(x)
  
}

# function for importing MZmine exported table of multiple
# peak characteristics into a table of just the specified
# characteristic
singleCharTable <- function(fullCharTable, idColNum,
                            firstInjColNum,
                            charNum, numChars) {
  
  # column number of first character column
  firstCharCol <- firstInjColNum + charNum - 1
  
  # total columns
  numTotCharCols <- ncol(fullCharTable)
  
  # character column sequence for specified characteristic
  charColNums <- seq(from = firstCharCol,
                     to = numTotCharCols,
                     by = numChars)
  
  # create subset of imported full characteristic table
  # with mz, rt, and specified characteristic
  # columns for specific characteristic
  dt.charTable <- subset(fullCharTable,
                         select = c(idColNum, charColNums))
  
  return(dt.charTable)
  
}

# function creates a long format data table with columns:
# featID, injection, value.name
featInjLongFormat <- function(x, featIDColNum, value.name, na.rm) {
  
  # x is the data table
  # featIDColNum is the column number of the featID column
  # which must immediately precede the injection columns
  # data in preceding columns will not be included in the
  # output data table
  # value.name is the name of the measured variable e.g. "pkHt"
  # for peak height measurements
  # na.rm as defined in melt - TRUE for false for removing NA rows
  
  # drop unneeded columns
  x <- x[,featIDColNum:ncol(x)]
  
  # convert to long form with melt
  x <- melt(x, id.vars = "featID",
            measure.vars = c(2:ncol(x)),
            variable.name = "injection",
            value.name = value.name,
            na.rm = na.rm)
  
  return(x)
  
}

# define function for writing tables to output folder
wrtTable <- function(table, folderPath, name) {
  # table = table object to write
  # folderPath = path of folder to write to
  # name = base name of output file, time stamp
  # and .csv are added
  
  # export table
  # add data and file extension to the file name
  OutputFileName <- paste0(folderPath,"/",
                           name,
                           startTimeStamp,".csv")
  
  #Write the adjusted dataframe to csv
  fwrite(table,
         file = OutputFileName,
         row.names = FALSE)
  
}


pdfPlotsSingleFile <- function(list.plotPage, orientation, folderPath, name) {
  
  # this function creates a pdf of the plots in list.plotPage
  # and exports them to the file path/name PdfFilePath
  
  # list.plotPage is the list of lists of plot objects
  # orientation is the orientation of the paper
  # "portrait" or "landscape"
  # folderPath = path of folder to write to
  # name = base name of output file, time stamp
  # and .csv are added
  
  PdfFilePath <- paste0(folderPath, "/", name,
                        startTimeStamp, ".pdf")
  
  # paper orientation
  if (orientation == "portrait") {
    paper <- "letter"
    width <- 7.5
    height <- 10
  } else if (orientation == "landscape") {
    paper <- "USr"
    width <- 10
    height <- 7.5
  } else {
    stop("please enter portrait or landscape orientation")
  }
  
  pdf(PdfFilePath, width = width, height = height, paper = paper, onefile = TRUE)
  for (i in seq(length(list.plotPage))) {
    grid.arrange(grobs = list.plotPage[[i]], ncol = 1)
    
    percComplete <- sprintf("%02.1f",round(i/length(list.plotPage) * 100, 1))
    
    flush.console()
    cat("\r", paste0(name, " pdf completion ", percComplete,"%"))
    flush.console()
  }
  dev.off()
}

# plot theme --------------------------------------------------

# create theme object for use in plots

theme.std <-
  theme(title = element_text(size = 10),
        axis.title = element_text(size = 9),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA, linetype = 1, size = .5),
        panel.grid.major.y = element_line(colour = "grey50", size = 0.1),
        panel.grid.major.x = element_blank())

# import table -----------------------------------------------------------

# table path
inputCSV <- paste0(parentDir,"/",inputCSV)

# import table
dt.inj <- fread(inputCSV, sep = ",", header = TRUE,
                check.names = FALSE)

# change column names
colnames(dt.inj)[1] <- "MZmineMz"
colnames(dt.inj)[2] <- "MZmineRt"

dt.inj <- repZeroEmpty(dt.inj)

# create featID column -------------------------------------

dt.inj <- addFeatID(dt.inj, 1, 2, 3)

# add "MZmine" to featID column header


# count features -----------------------------------------------

numFeatures <- nrow(dt.inj)

# create feature summary table -------------------------------

dt.featSummary <- subset(dt.inj, select = c("featID", "MZmineMz", "MZmineRt"))

# create individual characteristic tables -----------------------------

# m/z, RT, Ht
numChars <- 7

# use defined function to create individual
# characteristic tables
dt.PkMz <- singleCharTable(dt.inj, 3, 4,
                           1, numChars)

dt.PkRT <- singleCharTable(dt.inj, 3, 4,
                           2, numChars)

dt.PkRtStart <- singleCharTable(dt.inj, 3, 4,
                                3, numChars)

dt.PkRtEnd <- singleCharTable(dt.inj, 3, 4,
                              4, numChars)

dt.PkHt <- singleCharTable(dt.inj, 3, 4,
                           5, numChars)

dt.PkMzMin <- singleCharTable(dt.inj, 3, 4,
                              6, numChars)

dt.PkMzMax <- singleCharTable(dt.inj, 3, 4,
                              7, numChars)

# count injection files --------------------------------------------

numInjFiles <- ncol(dt.PkMz) - 1

# create long form tables ---------------------------------------------

# apply to tables
dt.PkMz.long <- featInjLongFormat(dt.PkMz, 1,
                                  "pkMz", na.rm = FALSE)

dt.PkRT.long <- featInjLongFormat(dt.PkRT, 1,
                                  "pkRT", na.rm = FALSE)

dt.PkRtStart.long <- featInjLongFormat(dt.PkRtStart, 1,
                                       "pkRtStart", na.rm = FALSE)

dt.PkRtEnd.long <- featInjLongFormat(dt.PkRtEnd, 1,
                                     "pkRtEnd", na.rm = FALSE)

dt.PkMzMin.long <- featInjLongFormat(dt.PkMzMin, 1,
                                     "pkMzMin", na.rm = FALSE)

dt.PkMzMax.long <- featInjLongFormat(dt.PkMzMax, 1,
                                     "pkMzMax", na.rm = FALSE)

dt.PkHt.long <- featInjLongFormat(dt.PkHt, 1,
                                  "pkHt", na.rm = FALSE)

# write dt.PkHt.long -----------------------------------------

# export this table for later comparison with the values
# extracted by the peakWindowExtractor script

# strip information added after file name
dt.PkHt.long[,injection := gsub("mzXML.*$","mzXML",injection)]


# write using defined function
wrtTable(dt.PkHt.long, outputFolderPath, "Table_featureInjectionHeight")

# create long peak table --------------------------------------

# start with dt.PkMz.long
dt.peak.long <- copy(dt.PkMz.long)

# remove missing peaks
dt.peak.long <- subset(dt.peak.long,
                       !is.na(pkMz))

# remove characters after file name in injection column
dt.peak.long[,injection := gsub("mzXML.*$","mzXML",injection)]
dt.PkRT.long[,injection := gsub("mzXML.*$","mzXML",injection)]

# add Rt
setkeyv(dt.peak.long, c("featID","injection"))
setkeyv(dt.PkRT.long, c("featID","injection"))
dt.peak.long[dt.PkRT.long, pkRT := pkRT]

# create a copy to use as a reference after values are changed
dt.peak.long_og <- copy(dt.peak.long)

# Peak characteristic statistics -------------------------------

# generate a table with the medians of each feature characteristic
# for each feature as well as the max and min of the RT and m/z
dt.PkMz.long$pkMz<-as.numeric(dt.PkMz.long$pkMz)
dt.medMz <- dt.PkMz.long[,median(pkMz, na.rm = TRUE),
                         keyby = .(featID)]
colnames(dt.medMz)[2] <- "medMz"
dt.PkRT.long$pkRT<-as.numeric(dt.PkRT.long$pkRT)
dt.medRT <- dt.PkRT.long[,median(pkRT, na.rm = TRUE),
                         keyby = .(featID)]
colnames(dt.medRT)[2] <- "medRT"

dt.sdRT <- dt.PkRT.long[,sd(pkRT, na.rm = TRUE),
                        keyby = .(featID)]
colnames(dt.sdRT)[2] <- "sdRT"
dt.PkRtStart.long$pkRtStart<-as.numeric(dt.PkRtStart.long$pkRtStart)
dt.medRtStart <- dt.PkRtStart.long[,median(pkRtStart, na.rm = TRUE),
                                   keyby = .(featID)]
colnames(dt.medRtStart)[2] <- "medRtStart"

# quantile RT start based on user input
dt.quantileRtStart <- dt.PkRtStart.long[,quantile(pkRtStart, probs = c(1 - rtQuantile), na.rm = TRUE),
                                        keyby = .(featID)]
colnames(dt.quantileRtStart)[2] <- "quantRtStart"

dt.PkRtEnd.long$pkRtEnd<-as.numeric(dt.PkRtEnd.long$pkRtEnd)
dt.medRtEnd <- dt.PkRtEnd.long[,median(pkRtEnd, na.rm = TRUE),
                               keyby = .(featID)]
colnames(dt.medRtEnd)[2] <- "medRtEnd"

# quantile RT end based on user input
dt.quantileRtEnd <- dt.PkRtEnd.long[,quantile(pkRtEnd, probs = c(rtQuantile), na.rm = TRUE),
                                    keyby = .(featID)]
colnames(dt.quantileRtEnd)[2] <- "quantRtEnd"

dt.minApexRt <- dt.PkRT.long[,min(pkRT, na.rm = TRUE),
                             keyby = .(featID)]
colnames(dt.minApexRt)[2] <- "minApexRt"

dt.maxApexRt <- dt.PkRT.long[,max(pkRT, na.rm = TRUE),
                             keyby = .(featID)]
colnames(dt.maxApexRt)[2] <- "maxApexRt"

dt.PkMzMin.long$pkMzMin<-as.numeric(dt.PkMzMin.long$pkMzMin)
dt.medMzMin <- dt.PkMzMin.long[,median(pkMzMin, na.rm = TRUE),
                               keyby = .(featID)]
colnames(dt.medMzMin)[2] <- "medMzMin"

dt.minMzMin <- dt.PkMzMin.long[,min(pkMzMin, na.rm = TRUE),
                               keyby = .(featID)]
colnames(dt.minMzMin)[2] <- "minMzMin"

dt.PkMzMax.long$pkMzMax<-as.numeric(dt.PkMzMax.long$pkMzMax)
dt.medMzMax <- dt.PkMzMax.long[,median(pkMzMax, na.rm = TRUE),
                               keyby = .(featID)]
colnames(dt.medMzMax)[2] <- "medMzMax"

dt.maxMzMax <- dt.PkMzMax.long[,max(pkMzMax, na.rm = TRUE),
                               keyby = .(featID)]
colnames(dt.maxMzMax)[2] <- "maxMzMax"

dt.minMz <- dt.PkMz.long[,min(pkMz, na.rm = TRUE),
                         keyby = .(featID)]
colnames(dt.minMz)[2] <- "minMz"

dt.maxMz <- dt.PkMz.long[,max(pkMz, na.rm = TRUE),
                         keyby = .(featID)]
colnames(dt.maxMz)[2] <- "maxMz"

dt.PkHt.long$pkHt<-as.numeric(dt.PkHt.long$pkHt)
dt.medHt <- dt.PkHt.long[,median(pkHt, na.rm = TRUE),
                         keyby = .(featID)]
colnames(dt.medHt)[2] <- "medHt"

dt.sumHt <- dt.PkHt.long[, sum(pkHt, na.rm = TRUE),
                         keyby = .(featID)]
colnames(dt.sumHt)[2] <- "sumHt"

dt.sumDetected <- dt.PkHt.long[,sum(!is.na(pkHt)),
                               keyby = .(featID)]
colnames(dt.sumDetected)[2] <- "sumDetected"

dt.sumDetected[, propDetect := sumDetected / numInjFiles]

# add to dt.featSummary
setkey(dt.featSummary, featID)
dt.featSummary[dt.medMz, medMz := medMz]
dt.featSummary[dt.medRT, medRT := medRT]
dt.featSummary[dt.sdRT, sdRT := sdRT]
dt.featSummary[dt.medRtStart, medRtStart := medRtStart]
dt.featSummary[dt.medRtEnd, medRtEnd := medRtEnd]
dt.featSummary[dt.quantileRtStart, quantRtStart := quantRtStart]
dt.featSummary[dt.quantileRtEnd, quantRtEnd := quantRtEnd]
dt.featSummary[dt.minApexRt, minApexRt := minApexRt]
dt.featSummary[dt.maxApexRt, maxApexRt := maxApexRt]
dt.featSummary[dt.medMzMax, medRtEnd := medRtEnd]
dt.featSummary[dt.medMzMin, medMzMin := medMzMin]
dt.featSummary[dt.minMzMin, minMzMin := minMzMin]
dt.featSummary[dt.medMzMax, medMzMax := medMzMax]
dt.featSummary[dt.maxMzMax, maxMzMax := maxMzMax]
dt.featSummary[dt.minMz, minMz := minMz]
dt.featSummary[dt.maxMz, maxMz := maxMz]
dt.featSummary[dt.medHt, medHt := medHt]
dt.featSummary[dt.sumHt, sumHt := sumHt]
dt.featSummary[dt.sumDetected, sumDetected := sumDetected]
dt.featSummary[dt.sumDetected, propDetect := propDetect]

rm(dt.medMz)
rm(dt.medRT)
rm(dt.sdRT)
rm(dt.medRtStart)
rm(dt.medRtEnd)
rm(dt.minApexRt)
rm(dt.maxApexRt)
rm(dt.medMzMin)
rm(dt.minMzMin)
rm(dt.medMzMax)
rm(dt.maxMzMax)
rm(dt.minMz)
rm(dt.maxMz)
rm(dt.medHt)
rm(dt.sumHt)
rm(dt.quantileRtEnd)
rm(dt.quantileRtStart)

# subset by sum detected and proportion detected ------------------------------

dt.featSummary <- dt.featSummary[sumDetected >=
                                   minTotalDetected][propDetect >=
                                                       minPropDetected]

# set feature hiearchy criteria ---------------------------

# add column for medHt * propDetect
dt.featSummary[, medHtPrpDet := medHt * propDetect]

# mz Window (initial) --------------------------------------

if (mzBoundOnPpm == TRUE) {
  
  #  calculate mz window size extremes based on medMz and user input
  dt.featSummary[,maxMzDiff := maxWindowPpm / 1000000 *  medMz]
  # dt.featSummary[,minMzDiff := minWindowPpm / 1000000 * medMz]
  
  # set window based on maxMzDiff
  dt.featSummary[, windMzMin := medMz - maxMzDiff]
  dt.featSummary[, windMzMax := medMz + maxMzDiff]
  
}

if (mzBoundOnPpm == FALSE) {
  
  # set mz window bounds based on exported m/z min and m/z max
  dt.featSummary[, windMzMin := minMzMin]
  dt.featSummary[, windMzMax := maxMzMax]
  
}


# RT Window (initial) ------------------------------------------

# Start using the minimum and maximum apexes to set the windows
# Provide a small buffer to ensure extreme observations are inside window
# rather than forming the boundary
dt.featSummary[, windRtMin := minApexRt - .001]
dt.featSummary[, windRtMax := maxApexRt + .001]


if (considerTails == TRUE) {
  
  # set RT window based on the input quantiles of the min and max reported
  # by MZmine or by the min or max of the RT of the apex,
  # whichever is more extreme.
  dt.featSummary[, windRtMin :=
                   ifelse(quantRtStart < windRtMin,
                          quantRtStart, windRtMin)]
  
  dt.featSummary[, windRtMax :=
                   ifelse(quantRtEnd > windRtMax,
                          quantRtEnd, windRtMax)]
  
}


# if RT window doesn't meet standard deviation requirements,
# use standard deviation to set.

if (!is.null(minWindowRtSd)) {
  
  dt.featSummary[, windRtMin :=
                   ifelse(windRtMin < medRT - minWindowRtSd * sdRT,
                          windRtMin,
                          medRT - minWindowRtSd * sdRT)]
  
  dt.featSummary[, windRtMax :=
                   ifelse(windRtMax > medRT + minWindowRtSd * sdRT,
                          windRtMax,
                          medRT + minWindowRtSd * sdRT)]
}

if (!is.null(maxWindowRtSd)) {
  
  dt.featSummary[, windRtMin :=
                   ifelse(windRtMin > medRT - maxWindowRtSd * sdRT,
                          windRtMin,
                          medRT - maxWindowRtSd * sdRT)]
  
  dt.featSummary[, windRtMax :=
                   ifelse(windRtMax < medRT + maxWindowRtSd * sdRT,
                          windRtMax,
                          medRT + maxWindowRtSd * sdRT)]
}


# Apply RT span requirements

if (!is.null(minWindowRtDur)) {
  
  minWindowRtDiff <- minWindowRtDur/2
  
  dt.featSummary[, windRtMin :=
                   ifelse(windRtMin < medRT - minWindowRtDiff,
                          windRtMin, medRT - minWindowRtDiff)]
  
  dt.featSummary[, windRtMax :=
                   ifelse(windRtMax > medRT + minWindowRtDiff,
                          windRtMax, medRT + minWindowRtDiff)]
}

if (!is.null(maxWindowRtDur)) {
  
  maxWindowRtDiff <- maxWindowRtDur/2
  
  dt.featSummary[, windRtMin :=
                   ifelse(windRtMin > medRT - maxWindowRtDiff,
                          windRtMin, medRT - maxWindowRtDiff)]
  
  dt.featSummary[, windRtMax :=
                   ifelse(windRtMax < medRT + maxWindowRtDiff,
                          windRtMax, medRT + maxWindowRtDiff)]
  
}


# rank by duration -------------------------------------------------

# add duration column
dt.featSummary[, windRtDur := windRtMax - windRtMin]

# sort by this column
setorder(dt.featSummary, -windRtDur)

# create column for rank
dt.featSummary[,rank := seq(from = 1, to = nrow(dt.featSummary))]

# write full feature summary ----------------------------------------

# write using defined function
wrtTable(dt.featSummary, outputFolderPath, "Table_featureSummary")

# function write window table ----------------------------------------------

WrtWindowTable <- function(dt.windowFinal, outputFolderPath, tableName) {
  
  # add window ID ------------------------------------------------------------
  
  # formats window table for writing and writes
  
  # replace featID with identifier based on RT start
  setorder(dt.windowFinal, windRtMin, na.last = TRUE)
  
  # create column for windowID - only give IDs to remaining windows
  dt.windowFinal[,windowID := paste0("window",
                                     sprintf("%09i",
                                             seq_len(nrow(dt.windowFinal))))]
  
  setcolorder(dt.windowFinal,
              c("windowID",
                "windMzMin",
                "windMzMax",
                "windRtMin",
                "windRtMax"))
  
  # round before writing ------------------------------------------------------
  
  dt.windowFinal[,windMzMin := round(windMzMin,5)]
  dt.windowFinal[,windMzMax := round(windMzMax,5)]
  dt.windowFinal[,windRtMin := round(windRtMin,4)]
  dt.windowFinal[,windRtMax := round(windRtMax,4)]
  
  # write peakWindowViewer tables ------------------------------------
  
  # write using defined function
  wrtTable(dt.windowFinal, outputFolderPath, tableName)
  
  return(dt.windowFinal)
  
}

dt.windowOne <- subset(dt.featSummary,
                       select = c("windMzMin",
                                  "windMzMax",
                                  "windRtMin",
                                  "windRtMax"))

# if no address overlap ---------------------------------------------------

if (addressOverlap == FALSE) {
  
  dt.windowFinal <- copy(dt.windowOne)
  
  tableName <- "Table_MzRtWindows"
  
  dt.windowFinal <- WrtWindowTable(dt.windowFinal, outputFolderPath, tableName)
}

# if address overlap --------------------------------------------

if (addressOverlap == TRUE) {
  
  # _write window table pre-adjustment ----------------------------------------
  
  tableName <- "Table_MzRtWindows_PreAdj"
  
  WrtWindowTable(dt.windowOne, outputFolderPath, tableName)
  
  # _WindowPlot function -----------------------------------------------------
  
  WindowPlot <- function(plot.review, numPlot, dt.windows,
                         lpdt.peak.long, dom.featID) {
    
    if (is.null(lpdt.peak.long)) { # null points table entered
      
      points <- NULL
      
    } else { # points table entered
      
      points <- geom_point(aes(x = pkRT, y = pkMz, shape = featID),
                           data = lpdt.peak.long)
      
      
    }
    
    # plot the apex points
    p <- ggplot() +
      points
    
    # add a layer for each window
    for (i_feat in 1:nrow(dt.windows)) {
      
      dt.windows.feat <- dt.windows[i_feat]
      
      lpdt.Window <- data.table(windMz = c(dt.windows.feat$windMzMax,
                                           dt.windows.feat$windMzMax,
                                           dt.windows.feat$windMzMin,
                                           dt.windows.feat$windMzMin,
                                           dt.windows.feat$windMzMax),
                                windRt = c(dt.windows.feat$windRtMin,
                                           dt.windows.feat$windRtMax,
                                           dt.windows.feat$windRtMax,
                                           dt.windows.feat$windRtMin,
                                           dt.windows.feat$windRtMin))
      
      p <- p +
        geom_path(aes(x = windRt, y = windMz),
                  data = lpdt.Window, linetype = 1)
      
    }
    
    
    
    titleText <- paste0(" dom:", dom.featID)
    
    p <- p +
      ggtitle(titleText) +
      theme.std
    
    
    plot.review[[numPlot]]  <- list()
    
    plot.review[[numPlot]][[1]] <- p
    
    return(plot.review)
    
  }
  
  # _while loop prep ------------------------------------------------------
  
  # sort by rank
  setorder(dt.featSummary, rank)
  
  # add overlap, RtSplit, and mergeToDom columns
  dt.featSummary[, overlap := 0]
  dt.featSummary[, RtSplit := 0]
  dt.featSummary[, mergeToDom := NA]
  
  # track number of merges and splits
  numWindowMerges <- 0
  numWindowRtSplits <- 0
  numWindowMzSplits <- 0
  
  # counter for plots
  numPlot <- 0
  
  # initialize plot list
  plot.review <- list()
  
  # initialize the window table to shrink
  dt.loopShrink <- copy(dt.featSummary)
  setorder(dt.loopShrink, rank)
  
  # initialize overlap counter
  numWindowsInitially <- nrow(dt.loopShrink)
  numWindowsRemain <- nrow(dt.loopShrink)
  
  # create empty table to add new windows to
  dt.windowGrow <- data.table(windMzMin = numeric(),
                              windMzMax = numeric(),
                              windRtMin = numeric(),
                              windRtMax = numeric())
  
  # loop iteration counter
  loopIteration <- 0
  
  # _while loop overlap ----------------------------------------------------
  
  while (numWindowsRemain > 0) { # run until dt.loopShrink is gone
    
    # update loop iteration counter and percent completion tracking
    loopIteration <- loopIteration + 1
    
    # pull window values from table
    dom.windMzMin <- dt.loopShrink$windMzMin[1]
    dom.windMzMax <- dt.loopShrink$windMzMax[1]
    dom.windRtMin <- dt.loopShrink$windRtMin[1]
    dom.windRtMax <- dt.loopShrink$windRtMax[1]
    dom.mzMin <- dt.loopShrink$minMz[1]
    dom.mzMax <- dt.loopShrink$maxMz[1]
    dom.RtMin <- dt.loopShrink$minApexRt[1]
    dom.RtMax <- dt.loopShrink$maxApexRt[1]
    dom.featID <- dt.loopShrink$featID[1]
    dom.medRT <- dt.loopShrink$medRT[1]
    dom.sdRT <- dt.loopShrink$sdRT[1]
    dom.rank <- dt.loopShrink$rank[1]
    dom.medMz <- dt.loopShrink$medMz[1]
    
    lpvar.mzDiff <- dt.loopShrink$medMz[1] * maxWindowPpm / 1000000
    
    # subset dt.loopShrink to overlapping mz windows
    lpdt.overlap <- subset(dt.loopShrink,
                           windMzMin <= dom.windMzMax)
    
    lpdt.overlap <- subset(lpdt.overlap,
                           windMzMax >= dom.windMzMin)
    
    # Expand window in RT direction to incorporate the RT bounds of any
    # overlapping window
    
    # initialize the while loop variables to start loop
    prev.windRtMin <- dom.windRtMin + 1
    prev.windRtMax <- dom.windRtMax - 1
    
    # __check for overlap ------------------------------------------------
    
    # continue checking for RT window overlaps until window stops expanding
    while (dom.windRtMin < prev.windRtMin | dom.windRtMax > prev.windRtMax) {
      
      # update loop variables
      prev.windRtMin <- dom.windRtMin
      prev.windRtMax <- dom.windRtMax
      
      # subset to overlap
      lpdt.overlap <- lpdt.overlap[windRtMax >= dom.windRtMin]
      lpdt.overlap <- lpdt.overlap[windRtMin <= dom.windRtMax]
      
      # set the new window RT limits to the RT edges of the group of
      # overlapping windows
      dom.windRtMin <- min(lpdt.overlap$windRtMin)
      dom.windRtMax <- max(lpdt.overlap$windRtMax)
      
    }
    
    numGroupOverlap <- nrow(lpdt.overlap)
    
    # remove features in lpdt.overlap from dt.loopShrink
    dt.loopShrink <- dt.loopShrink[!dt.loopShrink$featID %in% lpdt.overlap$featID]
    
    # __if overlap -------------------------------------------------
    
    # if there are overlapping windows, cluster for new windows
    if (numGroupOverlap > 1) {
      
      # subset dt.peak.long to features in lpdt.overlap
      lpdt.peak.long <- dt.peak.long[featID %in% lpdt.overlap$featID]
      
      # ___plot for review pre-fix-------------------------------------
      
      # advance plot counter
      numPlot <- numPlot + 1
      
      plot.review <- WindowPlot(plot.review, numPlot,
                                dt.windows = lpdt.overlap,
                                lpdt.peak.long, dom.featID)
      
      # ___cluster prep -------------------------------------------
      
      # create copy of lpdt.peak.long
      lpdt.cluster <- copy(lpdt.peak.long)
      
      # create scaled columns
      lpdt.cluster[,mzScaled := pkMz / (2 * lpvar.mzDiff)]
      lpdt.cluster[,rtScaled := pkRT / minWindowRtDur]
      
      # maximum number of potential clusters (based on number of
      # overlapping windows)
      maxNumCluster <- max(numGroupOverlap * 1.8, 5)
      
      # create McLust cluster object
      d_clust <- Mclust(lpdt.cluster[,c("rtScaled", "mzScaled")],
                        G = 1:maxNumCluster,
                        modelNames = c("EII", "VII", "EEI",
                                       "VEI", "EVI", "VVI"),
                        verbose = FALSE)
      
      # plot(d_clust, what = "classification")
      
      # add cluster to table
      lpdt.cluster[, cluster := d_clust[["classification"]]]
      
      # __create cluster windows ---------------------------------------------
      
      # find min and max rt and mz by cluster
      lpdt.minMz <- lpdt.cluster[, min(pkMz, na.rm = TRUE),
                                 keyby = .(cluster)]
      colnames(lpdt.minMz)[2] <- "minMz"
      
      lpdt.maxMz <- lpdt.cluster[, max(pkMz, na.rm = TRUE),
                                 keyby = .(cluster)]
      colnames(lpdt.maxMz)[2] <- "maxMz"
      
      lpdt.medMz <- lpdt.cluster[, median(pkMz, na.rm = TRUE),
                                 keyby = .(cluster)]
      colnames(lpdt.medMz)[2] <- "medMz"
      
      lpdt.minRt <- lpdt.cluster[, min(pkRT, na.rm = TRUE),
                                 keyby = .(cluster)]
      colnames(lpdt.minRt)[2] <- "minRt"
      
      lpdt.maxRt <- lpdt.cluster[, max(pkRT, na.rm = TRUE),
                                 keyby = .(cluster)]
      colnames(lpdt.maxRt)[2] <- "maxRt"
      
      lpdt.medRt <- lpdt.cluster[, median(pkRT, na.rm = TRUE),
                                 keyby = .(cluster)]
      colnames(lpdt.medRt)[2] <- "medRt"
      
      # create a new data table for window information
      lpdt.clusterWindow <- copy(lpdt.minMz)
      
      # add other columns
      lpdt.clusterWindow[lpdt.maxMz, maxMz := maxMz]
      lpdt.clusterWindow[lpdt.medMz, medMz := medMz]
      lpdt.clusterWindow[lpdt.minRt, minRt := minRt]
      lpdt.clusterWindow[lpdt.maxRt, maxRt := maxRt]
      lpdt.clusterWindow[lpdt.medRt, medRt := medRt]
      
      rm(lpdt.minMz,
         lpdt.maxMz,
         lpdt.medMz,
         lpdt.minRt,
         lpdt.maxRt,
         lpdt.medRt)
      
      # create window RT bounds with small buffers
      lpdt.clusterWindow[, windRtMin := round(minRt - 0.0002,4)]
      lpdt.clusterWindow[, windRtMax := round(maxRt + 0.0002,4)]
      
      # enforce window minimum duration
      lpdt.clusterWindow[, windRtMin := pmin(windRtMin,
                                             round(medRt - minWindowRtDur/2,4))]
      lpdt.clusterWindow[, windRtMax := pmax(windRtMax,
                                             round(medRt + minWindowRtDur/2,4))]
      
      # create window mz bounds using extreme points median and mz based on ppm
      lpdt.clusterWindow[, windMzMin := pmin(round(minMz - 0.00002,5),
                                             (medMz - lpvar.mzDiff/2))]
      lpdt.clusterWindow[, windMzMax := pmax(round(minMz + 0.00002,5),
                                             (medMz + lpvar.mzDiff/2))]
      
      # __filter cluster windows by cluster population ----------------------
      
      # find number of points per cluster
      lpdt.clusterPop <- lpdt.cluster[, .N, by = cluster]
      
      # find clusters with fewer than minimum points
      lpdt.clustRm <- lpdt.clusterPop[N < minApexPerCluster]
      
      # if clusters have fewer than minimum points, remove them
      if (nrow(lpdt.clustRm) > 0) {
        
        # remove clusters with fewer than minimum points
        lpdt.clusterWindow <- lpdt.clusterWindow[!lpdt.clusterWindow$cluster %in%
                                                   lpdt.clustRm$cluster]
      }
      
      # __Filter cluster windows by mixing proportion ------------------------
      
      # min mixing proportion per cluster
      minMixProportion <- 0.2 / numGroupOverlap
      
      lpdt.clustMixPro <- data.table(cluster = 1:d_clust[["G"]],
                                     mixPro = d_clust[["parameters"]][["pro"]])
      
      # find clusters with lower than minimum mixing proportion
      lpdt.clustRm <- lpdt.clustMixPro[mixPro < minMixProportion]
      
      # if clusters have lower than minimum mixing proportion, remove them
      if (nrow(lpdt.clustRm) > 0) {
        
        # remove clusters with lower than minimum mixing proportion
        lpdt.clusterWindow <- lpdt.clusterWindow[!lpdt.clusterWindow$cluster %in%
                                                   lpdt.clustRm$cluster]
        
      }
      
      # __combine overlapping cluster windows --------------------------------
      
      numClusterRemain <- nrow(lpdt.clusterWindow)
      
      # create window table
      lpdt.window <- subset(lpdt.clusterWindow,
                            select = c("windMzMin", "windMzMax",
                                       "windRtMin", "windRtMax"))
      
      # initialize overlap counter to start while loop
      numOverlap <- 1
      
      while (numOverlap  > 0) { # run until lpdt.clusterWindow is gone
        
        # reset overlap counter
        numOverlap  <- 0
        
        # copy lpdt.window to loop through
        lpdt.clusterGroupShrink <- copy(lpdt.window)
        
        # reset lpdt.window
        lpdt.window <- data.table(windMzMin = numeric(),
                                  windMzMax = numeric(),
                                  windRtMin = numeric(),
                                  windRtMax = numeric())
        
        # add duration column
        lpdt.clusterGroupShrink[, windRtDur := windRtMax - windRtMin]
        
        # sort by duration
        setorder(lpdt.clusterGroupShrink, -windRtDur)
        
        # add rank
        lpdt.clusterGroupShrink[, rank := 1:nrow(lpdt.clusterGroupShrink)]
        
        while (nrow(lpdt.clusterGroupShrink) > 0) {
          
          # pull window values from table top row
          domClus.windMzMin <- lpdt.clusterGroupShrink$windMzMin[1]
          domClus.windMzMax <- lpdt.clusterGroupShrink$windMzMax[1]
          domClus.windRtMin <- lpdt.clusterGroupShrink$windRtMin[1]
          domClus.windRtMax <- lpdt.clusterGroupShrink$windRtMax[1]
          domClus.rank <- lpdt.clusterGroupShrink$rank[1]
          
          # subset dt.loopShrink to overlapping windows
          lpdt.overlapCluster <- subset(lpdt.clusterGroupShrink,
                                        windMzMin <= domClus.windMzMax)
          if (nrow(lpdt.overlapCluster) > 0) {
            lpdt.overlapCluster <- subset(lpdt.overlapCluster,
                                          windMzMax >= domClus.windMzMin)
          }
          if (nrow(lpdt.overlapCluster) > 0) {
            lpdt.overlapCluster <- subset(lpdt.overlapCluster,
                                          windRtMin <= domClus.windRtMax)
          }
          if (nrow(lpdt.overlapCluster) > 0) {
            lpdt.overlapCluster <- subset(lpdt.overlapCluster,
                                          windRtMax >= domClus.windRtMin)
            
          }
          if (nrow(lpdt.overlapCluster) > 1) {
            
            # increase overlap counter
            numOverlap  <- numOverlap + 1
            
            # remove features in lpdt.overlap from lpdt.clusterGroupShrink
            lpdt.clusterGroupShrink <- lpdt.clusterGroupShrink[!lpdt.clusterGroupShrink$rank %in%
                                                                 lpdt.overlapCluster$rank]
            
            # create table row of new merged window boundaries
            lpdt.windowClustMerge <- data.table(
              windMzMin = min(lpdt.overlapCluster$windMzMin),
              windMzMax = max(lpdt.overlapCluster$windMzMax),
              windRtMin = min(lpdt.overlapCluster$windRtMin),
              windRtMax = max(lpdt.overlapCluster$windRtMax)
            )
            
            # add to growing lpdt.window
            lpdt.window <- rbind(lpdt.window, lpdt.windowClustMerge)
            
          } else {# if no overlap
            
            # remove top row from table
            lpdt.clusterGroupShrink <- subset(lpdt.clusterGroupShrink,
                                              rank != domClus.rank)
            
            # add dominant window values to lpdt.window
            lpdt.windowClustMerge <- data.table(
              windMzMin = domClus.windMzMin,
              windMzMax = domClus.windMzMax,
              windRtMin = domClus.windRtMin,
              windRtMax = domClus.windRtMax
            )
            
            # add to growing lpdt.window
            lpdt.window <- rbind(lpdt.window, lpdt.windowClustMerge)
            
          } # end overlap if/else
          
        } # end cluster table shrink while loop
        
      } # end cluster merge while loop
      
      # _plot new windows -----------------------------------------------
      
      # advance plot counter
      numPlot <- numPlot + 1
      
      plot.review <- WindowPlot(plot.review, numPlot,
                                dt.windows = lpdt.window,
                                lpdt.peak.long, dom.featID)
      
      # _if no overlaps ------------------------------------------------
      
    } else {
      
      # if no overlaps create from single lpdt.overlap entry
      lpdt.window <- subset(lpdt.overlap,
                            select = c("windMzMin",
                                       "windMzMax",
                                       "windRtMin",
                                       "windRtMax"))
      
    } # end original group overlap if/else statement
    
    # _add new windows to dt.windowGrow ---------------------------------
    
    # add new windows to dt.windowGrow
    dt.windowGrow <- rbind(dt.windowGrow, lpdt.window)
    
    # update numWindowsRemain
    numWindowsRemain <- nrow(dt.loopShrink)
    
    percComplete <- round((numWindowsInitially - numWindowsRemain)/numWindowsInitially * 100,1)
    
    print(paste0("cluster loop completion: ", percComplete,"%"))
    
  } # end overlap while loop
  
  summaryText <-
    updateSummary("end cluster loop")
  
  # _PDF review plots -------------------------------------------
  
  if (plotWindows == TRUE) {
    
    if (length(plot.review) > 0) {
      
      pdfPlotsSingleFile(plot.review, "landscape",
                         outputFolderPath, "reviewOverlapPlots")
    }
  }
  
  # _Final merge overlap loop ------------------------------------------------------
  
  # to ensure no overlaps that weren't checked as previous grouping
  
  # create window table
  lpdt.window2 <- copy(dt.windowGrow)
  
  # track window only plots
  numPlotWindOnly <- 0
  
  # initialize plot list
  plot.reviewWindOnly <- list()
  
  # initialize overlap counter to start while loop
  numOverlap <- 1
  
  print("begin final merge loop")
  
  while (numOverlap  > 0) {
    
    # reset overlap counter
    numOverlap  <- 0
    
    # copy lpdt.window2 to loop through
    lpdt.finalMergeShrink <- copy(lpdt.window2)
    
    # reset lpdt.window2
    lpdt.window2 <- data.table(windMzMin = numeric(),
                               windMzMax = numeric(),
                               windRtMin = numeric(),
                               windRtMax = numeric())
    
    # add duration column
    lpdt.finalMergeShrink[, windRtDur := windRtMax - windRtMin]
    
    # sort by duration
    setorder(lpdt.finalMergeShrink, -windRtDur)
    
    # add rank
    lpdt.finalMergeShrink[, rank := 1:nrow(lpdt.finalMergeShrink)]
    
    while (nrow(lpdt.finalMergeShrink) > 0) {
      
      # pull window values from table top row
      domClus.windMzMin <- lpdt.finalMergeShrink$windMzMin[1]
      domClus.windMzMax <- lpdt.finalMergeShrink$windMzMax[1]
      domClus.windRtMin <- lpdt.finalMergeShrink$windRtMin[1]
      domClus.windRtMax <- lpdt.finalMergeShrink$windRtMax[1]
      domClus.rank <- lpdt.finalMergeShrink$rank[1]
      
      # subset dt.loopShrink to overlapping windows
      lpdt.overlapFinal <- subset(lpdt.finalMergeShrink,
                                  windMzMin <= domClus.windMzMax)
      if (nrow(lpdt.overlapFinal) > 0) {
        lpdt.overlapFinal <- subset(lpdt.overlapFinal,
                                    windMzMax >= domClus.windMzMin)
      }
      if (nrow(lpdt.overlapFinal) > 0) {
        lpdt.overlapFinal <- subset(lpdt.overlapFinal,
                                    windRtMin <= domClus.windRtMax)
      }
      if (nrow(lpdt.overlapFinal) > 0) {
        lpdt.overlapFinal <- subset(lpdt.overlapFinal,
                                    windRtMax >= domClus.windRtMin)
        
      }
      if (nrow(lpdt.overlapFinal) > 1) {
        
        # increase overlap counter
        numOverlap  <- numOverlap + 1
        
        # remove features in lpdt.overlap from lpdt.finalMergeShrink
        lpdt.finalMergeShrink <- lpdt.finalMergeShrink[!lpdt.finalMergeShrink$rank %in%
                                                         lpdt.overlapFinal$rank]
        
        # __plot windows pre-merge ----------------------------------
        
        # advance plot counter
        numPlotWindOnly <- numPlotWindOnly + 1
        
        plot.reviewWindOnly <- WindowPlot(plot.reviewWindOnly, numPlotWindOnly,
                                          dt.windows = lpdt.overlapFinal,
                                          lpdt.peak.long = NULL,
                                          dom.featID = domClus.windMzMax)
        
        # __merge -----------------------------------------------------------
        
        # create table row of new merged window boundaries
        lpdt.windowFinalMerge <- data.table(
          windMzMin = min(lpdt.overlapFinal$windMzMin),
          windMzMax = max(lpdt.overlapFinal$windMzMax),
          windRtMin = min(lpdt.overlapFinal$windRtMin),
          windRtMax = max(lpdt.overlapFinal$windRtMax)
        )
        
        # add to growing lpdt.window2
        lpdt.window2 <- rbind(lpdt.window2, lpdt.windowFinalMerge)
        
        # _plot post merge -------------------------------------------------
        
        # advance plot counter
        numPlotWindOnly <- numPlotWindOnly + 1
        
        plot.reviewWindOnly <- WindowPlot(plot.reviewWindOnly, numPlotWindOnly,
                                          dt.windows = lpdt.windowFinalMerge,
                                          lpdt.peak.long = NULL,
                                          dom.featID = domClus.windMzMax)
        
      } else {# if no overlap
        
        # remove top row from table
        lpdt.finalMergeShrink <- subset(lpdt.finalMergeShrink,
                                        rank != domClus.rank)
        
        # add dominant window values to lpdt.window2
        lpdt.windowFinalMerge <- data.table(
          windMzMin = domClus.windMzMin,
          windMzMax = domClus.windMzMax,
          windRtMin = domClus.windRtMin,
          windRtMax = domClus.windRtMax
        )
        
        # add to growing lpdt.window2
        lpdt.window2 <- rbind(lpdt.window2, lpdt.windowFinalMerge)
        
      } # end overlap if/else
      
    } # end shrink table while loop
    
    print(paste0("remaining window overlaps: ", numOverlap))
    
  } # end final merge while loop
  
  summaryText <-
    updateSummary("end final merge loop")
  
  # assign data table to post final merge
  dt.windowFinal <- copy(lpdt.window2)
  
  # _PDF review plots -------------------------------------------
  
  if (plotWindows == TRUE) {
    
    if (length(plot.reviewWindOnly) > 0) {
      
      pdfPlotsSingleFile(plot.reviewWindOnly, "landscape",
                         outputFolderPath, "reviewFinalMergePlots")
      
    }
  }
  
  tableName <- "Table_MzRtWindows"
  
  dt.windowFinal <- WrtWindowTable(dt.windowFinal, outputFolderPath, tableName)
  
} # end address overlap if statement

# context window table -------------------------------------------------------

# create table with context (more of chromatogram)
dt.mzMinePWVcontext <- copy(dt.windowFinal)

dt.mzMinePWVcontext[, windowID := paste0(windowID,"context")]
dt.mzMinePWVcontext[, windRtMin := ifelse(windRtMin - rtContextWindowBuffer < 0,
                                          0,windRtMin - rtContextWindowBuffer)]
dt.mzMinePWVcontext[, windRtMax := windRtMax + rtContextWindowBuffer]

# normal window table with context rows added
dt.mzMinePWVcontext <- rbind(dt.windowFinal, dt.mzMinePWVcontext,
                             use.names = FALSE)

# write using defined function
wrtTable(dt.mzMinePWVcontext, outputFolderPath,
         "Table_MZminePeakWindowViewerContext")

# Plot Window Stats -----------------------------------------------------------

dt.windowFinalStats <- copy(dt.windowFinal)

# add mzSpan, rtDur, midMz, midRt columns
dt.windowFinalStats[, midMz := (windMzMax + windMzMin)/2]
dt.windowFinalStats[, midRt := (windRtMax + windRtMin)/2]
dt.windowFinalStats[, mzSpan := windMzMax - windMzMin]
dt.windowFinalStats[, rtDur := windRtMax - windRtMin]

wrtTable(dt.windowFinalStats, outputFolderPath, "Table_MzRtWindows_addDims")

# _plot duration distribution ------------------------------------------------

plot.dur <-
  ggplot() +
  geom_density(aes(rtDur), data = dt.windowFinalStats) +
  geom_freqpoly(aes(rtDur, stat(density)), data = dt.windowFinalStats,
                bins = 100, color = "gray") +
  xlim(0,NA) +
  theme.std

ggsave(filename = "plot_windowDurationDistribution.pdf",
       plot = plot.dur,
       device = "pdf",
       path = outputFolderPath,
       width = 11,
       height = 8.5,
       units = "in")

# _plot mzSpan distribution ------------------------------------------------

plot.mzSpan <-
  ggplot() +
  geom_density(aes(mzSpan), data = dt.windowFinalStats) +
  geom_freqpoly(aes(mzSpan, stat(density)), data = dt.windowFinalStats,
                bins = 100, color = "gray") +
  xlim(0,NA) +
  theme.std

ggsave(filename = "plot_windowMZSpanDistribution.pdf",
       plot = plot.mzSpan,
       device = "pdf",
       path = outputFolderPath,
       width = 11,
       height = 8.5,
       units = "in")

# _plot RT distribution ------------------------------------------------------

plot.midRt <-
  ggplot() +
  geom_density(aes(midRt), data = dt.windowFinalStats) +
  geom_freqpoly(aes(midRt, stat(density)), data = dt.windowFinalStats,
                bins = 100, color = "gray") +
  xlim(0,NA) +
  theme.std

ggsave(filename = "plot_windowRtDistribution.pdf",
       plot = plot.midRt,
       device = "pdf",
       path = outputFolderPath,
       width = 11,
       height = 8.5,
       units = "in")

# _plot mzRt Scatter ----------------------------------------------------------

plot.mzRtScatter <-
  ggplot() +
  geom_point(aes(x = midRt, y = midMz), data = dt.windowFinalStats,
             shape = ".") +
  xlim(0,NA) +
  theme.std

ggsave(filename = "plot_mzRtScatter.pdf",
       plot = plot.mzRtScatter,
       device = "pdf",
       path = outputFolderPath,
       width = 11,
       height = 8.5,
       units = "in")


