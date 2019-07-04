# Description -----------------------------------------------------------------

# create training, validation and test sets from labels and corresponding
# image matrices

# User Input ------------------------------------------------------------------

# Entering file paths:
# for Windows operating system use \\ or / instead of \ to separate directories
# example:
# ExampleDir <- 'M:\\folder1\\folderWithFile'


# path to folder where output directory will be created
outPutParentDir <- "E:\\Example\\Example"

# data location table
# CSV file that includes the locations of the peak group matrices in column 1
# and the name of the corresponding label CSV file. The label CSV files should
# be in the output parent directory and standard format
CSV_dataLocation <- "dataLocation.csv"

# proportion of cases to be used for training and validation/calibration
# remaining will be used for testing
propTrain <- 0.5
propValid <- 0.25

# number of shuffled versions of training matrices to use
# for expanding the training data
# shuffling will be applied to match the minority class size to the
# majority class. This is the number of additional shuffles
numShuffle <- 1

# End user input ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

# advanced user input (probably don't need to change) ------------------------

# text to remove from the window identifier (e.g. ".jpg") in label table to
# match with matrix csv file names which don't include extension
removeText <- ".jpg"

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

OutputFolderName <- paste0("peakNeuralNetBuildData")

if (restart == TRUE) {

  outputFolderPath <- existingOutputPath

} else {

  # create output folder if it doesn't exist
  outputFolderPath <- paste0(outPutParentDir,"/",
                             OutputFolderName, startTimeStamp)

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

# import location table ----------------------------------------------------

filePath <- paste0(outPutParentDir, "/", CSV_dataLocation)

dt.location <- fread(file = filePath)

colnames(dt.location)[1] <- "dirPath"
colnames(dt.location)[2] <- "labelCSV"

# functions ------------------------------------------------------------------

# _import function ----------------------------------------------------------
# function for reading and converting to matrix
ReadToMatrix <- function(path) {

  # load
  data <- fread(file = path, header = FALSE)

  # replace zeros removed in export
  data[is.na(data)] <- 0

  # convert to matrix
  data <- as.matrix(data)

}

# _functions for expanding data by mirroring and shuffling ---------------------

# function for mirroring matrix rows
mirrorMatrix <- function(mat) {
  # mat is the matrix to be mirroed

  matMirror <- mat[, ncol(mat):1]

  return(matMirror)

}

# function for shuffling matrix rows but leaving the last row which
# indicates the location of the window in the same location
shuffleMatrixRows <- function(mat) {
  # mat is the matrix to be shuffled - the final row remains in place

  matRowShuffle <- mat[c(sample(nrow(mat) - 1), nrow(mat)),]

  return(matRowShuffle)

}

# function to expand training and validation data by mirroring and shuffling
# rows using the above functions
expandExamplesMirrorShuffle <- function(matList, numShuffle) {
  # matList is list of matrices
  # numShuffle is the number of times the original list is duplicated by
  # creating row-shuffled versions of the original matrices

  # create list of mirrored matrices
  matListMirror <- lapply(matList, mirrorMatrix)

  # start list of combined matrix lists
  matListOgAndMirror <- c(matList, matListMirror)

  # copy to add shuffled lists
  matListCombined <- matListOgAndMirror

  # create numShuffled shuffled versions of matListOgAndMirror and
  # add to matListCombined
  for (i in 1:numShuffle) {

    # create shuffled matrix version of original matrix list
    matListShuffle <- lapply(matListOgAndMirror, shuffleMatrixRows)

    # add to combined list
    matListCombined <- c(matListCombined, matListShuffle)

  } # end for-loop

  return(matListCombined)

}

# _match train class size -----------------------------------------------------

matchTrainClassSize <- function(dt.minority, array.minority,
                                dt.majority) {

  # more bad examples, add shuffled good examples

  increaseFactor <- nrow(dt.majority) / nrow(dt.minority)

  increaseFactorRndUp <- ceiling(increaseFactor)

  matListCombined <- array.minority

  for (i in 1:increaseFactorRndUp) {

    # create shuffled matrix version of original matrix list
    matListShuffle <- lapply(array.minority, shuffleMatrixRows)

    # add to combined list
    matListCombined <- c(matListCombined, matListShuffle)

  }

  # update label table with additional shuffle rows

  dt.trainCombined <- copy(dt.minority)

  for (i in 1:increaseFactorRndUp) {

    dt.trainShuffle <- copy(dt.minority)

    dt.trainCombined <- rbind(dt.trainCombined, dt.trainShuffle)

  }

  # trim to match size of majority class
  dt.trainCombined <- dt.trainCombined[1:nrow(dt.majority)]
  matListCombined <-  matListCombined[1:nrow(dt.majority)]

  fnOutputList <- list(dt.trainCombined, matListCombined)
  return(fnOutputList)

}




outPutIndivSets <- function(setNum, dt.location){

  peakMatrixParentDir <- dt.location$dirPath[setNum]

  # CSV file with window labels
  # needs to be located in the outPutParentDir
  CSV_labels <- dt.location$labelCSV[setNum]

  # import label table ------------------------------------------------------------

  # Column number of window identifier - must match with matrix file name
  colNumWindowID <- 2
  # column number viewed
  colNumViewed <- 3
  # column number with classification label, binary 1 for good 0 for bad
  colNumLabel <- 4
  # column number indicating window bound issues
  colNumWindowBad <- 5
  # column number indicating if peak is borderline
  colNumBorderline <- 6

  importFilePath <- paste0(outPutParentDir, "/", CSV_labels)
  importColNumVect <- c(colNumWindowID, colNumViewed, colNumLabel,
                        colNumWindowBad, colNumBorderline)

  dt.Labels <- fread(file = importFilePath, header = TRUE,
                     select = importColNumVect)

  # rename columns
  colnames(dt.Labels) <- c("windName", "viewed", "label",
                           "windowBad", "borderline")

  dt.Labels[, windName := gsub(removeText, "", windName)]

  # window review summary stats -------------------------------------------------------

  totalWindowsInList <- nrow(dt.Labels)
  totalWindowsViewed <- sum(dt.Labels$viewed)
  totalWindowsBorderline <- sum(dt.Labels$borderline)
  totalWindowsWindowBad <- sum(dt.Labels$windowBad)

  # subset to useable images, needed columns ----------------------------------------------------

  dt.Labels <- dt.Labels[viewed == 1][windowBad == 0][borderline == 0]

  dt.Labels <- unique(dt.Labels, by = "windName")

  dt.Labels <- subset(dt.Labels,
                      select = c("windName", "label"))

  # window review summary table ------------------------------------------------------------

  totalWindowsGoodUsable <- sum(dt.Labels$label)
  totalWindowsBadUsable <- nrow(dt.Labels) - totalWindowsGoodUsable

  vect.statNames <- c("totalWindowsInList",
                      "totalWindowsViewed",
                      "totalWindowsBorderline",
                      "totalWindowsWindowBad",
                      "totalWindowsGoodUsable",
                      "totalWindowsBadUsable")

  vect.stats <- c(totalWindowsInList,
                  totalWindowsViewed,
                  totalWindowsBorderline,
                  totalWindowsWindowBad,
                  totalWindowsGoodUsable,
                  totalWindowsBadUsable)

  dt.windowReviewSummary <- data.table(stat = vect.statNames,
                                       value = vect.stats)

  WrtTable(dt.windowReviewSummary, outputFolderPath,
           paste0("windowReviewSummaryTable",setNum))

  # Determine files to load -----------------------------------------------

  # create a list of peak window tables
  fileListPath <- list.files(peakMatrixParentDir, full.names = TRUE,
                             pattern = ".csv")

  # make table
  dt.peakFiles <- data.table(path = fileListPath)

  # create column with just window names
  dt.peakFiles[, windName := gsub("^.*/", "", path)]
  dt.peakFiles[, windName := gsub("\\.csv", "", windName)]

  # subset to windows with both matrix files and labels
  dt.windowsInBoth <- fintersect(dt.peakFiles[,list(windName)],
                                 dt.Labels[,list(windName)])
  dt.peakFiles <- dt.peakFiles[windName %in% dt.windowsInBoth$windName]

  rm(dt.windowsInBoth)

  # add labels to dt.peakFiles
  setkey(dt.peakFiles, windName)
  setkey(dt.Labels, windName)
  dt.peakFiles[dt.Labels, label := label]

  rm(dt.Labels)

  # create subsets for good and bad
  dt.LabelsGood <- dt.peakFiles[label == 1]
  dt.LabelsBad <- dt.peakFiles[label == 0]

  # find the number of cases of the minority class
  numRowMinorityClass <- min(nrow(dt.LabelsGood),
                             nrow(dt.LabelsBad))

  # determine number of rows of each class in each set
  numRowTrain <- ceiling(numRowMinorityClass * propTrain)
  numRowValidate <- ceiling(numRowMinorityClass * propValid)
  numRowTest <- numRowMinorityClass - numRowTrain - numRowValidate

  # create random vectors
  vect.shuffle.good <- sample(nrow(dt.LabelsGood))
  vect.shuffle.bad <- sample(nrow(dt.LabelsBad))

  # add shuffle vectors to tables and sort by it
  dt.LabelsGood[, shuffle := vect.shuffle.good]
  setkey(dt.LabelsGood, shuffle)

  dt.LabelsBad[, shuffle := vect.shuffle.bad]
  setkey(dt.LabelsBad, shuffle)

  validStart <- numRowTest + 1
  validEnd <- numRowTest + numRowValidate

  # create column indicating train validate test
  dt.LabelsGood[, useCat := ifelse(shuffle <= numRowTest, 3,
                                   ifelse(shuffle %between% c(validStart,validEnd),
                                          2,1))]

  dt.LabelsBad[, useCat := ifelse(shuffle <= numRowTest, 3,
                                  ifelse(shuffle %between% c(validStart,validEnd),
                                         2,1))]

  # combine for export
  dt.trainValidTest <- rbind(dt.LabelsGood, dt.LabelsBad,
                             use.names = TRUE, fill = FALSE,
                             idcol = NULL)

  # export for record
  outputFilePath <- paste0(outputFolderPath,"/","trainValidTestTable",setNum,".csv")
  fwrite(dt.trainValidTest, file = outputFilePath)


  # read and export test and validation ----------------------------------------------------

  # create table of test cases
  dt.test <- dt.trainValidTest[useCat == 3]
  dt.valid <- dt.trainValidTest[useCat == 2]

  # read in tables to list
  array.Test <- lapply(dt.test$path, ReadToMatrix)
  array.Valid <- lapply(dt.valid$path, ReadToMatrix)

  # convert to array
  array.Test <- abind(array.Test, along = 0, force.array = TRUE)
  array.Valid <- abind(array.Valid, along = 0, force.array = TRUE)

  # export arrays
  outputFilePath <- paste0(outputFolderPath,"/","TestData",setNum,".rds")
  saveRDS(array.Test, outputFilePath)

  outputFilePath <- paste0(outputFolderPath,"/","ValidationData",setNum,".rds")
  saveRDS(array.Valid , outputFilePath)

  # export tables
  outputFilePath <- paste0(outputFolderPath,"/","TestLabels",setNum,".csv")
  fwrite(dt.test, file = outputFilePath)

  outputFilePath <- paste0(outputFolderPath,"/","ValidationLabels",setNum,".csv")
  fwrite(dt.valid, file = outputFilePath)

  rm(array.Test)
  rm(array.Valid)

  logText <- UpdateLogText(logText,"test and validation data exported")
  UpdateLogFile(logFilePath, logText)


  # load training and validation sets ---------------------------------------

  # create table of train cases
  dt.trainGood <- dt.trainValidTest[useCat == 1][label == 1]
  dt.trainBad <- dt.trainValidTest[useCat == 1][label == 0]

  # read in tables to list
  array.TrainGood <- lapply(dt.trainGood$path, ReadToMatrix)
  array.TrainBad <- lapply(dt.trainBad$path, ReadToMatrix)

  # apply function to increase minority class
  if (nrow(dt.trainBad) > nrow(dt.trainGood)) {

    fnList <- matchTrainClassSize(dt.trainGood, array.TrainGood,
                                  dt.trainBad)

    dt.trainGood <- fnList[[1]]
    array.TrainGood <- fnList[[2]]

    rm(fnList)

  } else {

    if (nrow(dt.trainGood) > nrow(dt.trainBad)) {

      fnList <- matchTrainClassSize(dt.trainBad, array.TrainBad,
                                    dt.trainGood)

      dt.trainBad <- fnList[[1]]
      array.TrainBad <- fnList[[2]]

      rm(fnList)

    }
  }

  # combine good and bad train data ---------------------------------------------

  array.Train <- c(array.TrainGood, array.TrainBad)
  dt.train <- rbind(dt.trainGood, dt.trainBad)

  rm(dt.trainBad)
  rm(dt.trainGood)
  rm(array.TrainGood)
  rm(array.TrainBad)

  # Expand training  data (Apply functions) --------------------------------------------

  array.Train <- expandExamplesMirrorShuffle(array.Train, numShuffle)

  # Export training data --------------------------------------------------------

  # convert to array format
  array.Train <- abind(array.Train, along = 0, force.array = TRUE)

  # export array
  outputFilePath <- paste0(outputFolderPath,"/","TrainingData",setNum,".rds")
  saveRDS(array.Train, outputFilePath)

  # expand training label table to account for matrix example expansion
  dt.trainOgAndMirror <- rbind(dt.train, dt.train)

  dt.trainCombined <- copy(dt.trainOgAndMirror)

  for (i in 1:numShuffle) {

    dt.trainCombined <- rbind(dt.trainCombined, dt.trainOgAndMirror)

  }

  outputFilePath <- paste0(outputFolderPath,"/","TrainingLabels",setNum,".csv")
  fwrite(dt.trainCombined, file = outputFilePath)

}

# apply function create output for each input set --------------------------

setSeq <- 1:nrow(dt.location)

lapply(setSeq, outPutIndivSets, dt.location)

# create combined tables ---------------------------------------------------

combineSetTables <- function(outputFolderPath, tableName, numSets) {


  tableList <- list()

  for (i_tableNum in 1:numSets) {

    filePath <- paste0(outputFolderPath,"/", tableName, i_tableNum, ".csv")

    tableList[[i_tableNum]] <- fread(filePath)

    file.remove(filePath)

  }

  combinedTable <- rbindlist(tableList)

  return(combinedTable)

}

numSets <- nrow(dt.location)

dt.TestLabels <- combineSetTables(outputFolderPath, "TestLabels", numSets)
WrtTable(dt.TestLabels, outputFolderPath, "TestLabels")
rm(dt.TestLabels)

dt.ValidationLabels <- combineSetTables(outputFolderPath,
                                        "ValidationLabels", numSets)
WrtTable(dt.ValidationLabels, outputFolderPath, "ValidationLabels")
rm(dt.ValidationLabels)

dt.TrainingLabels <- combineSetTables(outputFolderPath,
                                      "TrainingLabels", numSets)
WrtTable(dt.TrainingLabels, outputFolderPath, "TrainingLabels")
rm(dt.TrainingLabels)

dt.trainValidTestTable <- combineSetTables(outputFolderPath,
                                           "trainValidTestTable", numSets)
WrtTable(dt.trainValidTestTable, outputFolderPath, "trainValidTestTable")
rm(dt.trainValidTestTable)

dt.windowReviewSummaryTable <- combineSetTables(outputFolderPath,
                                                "windowReviewSummaryTable",
                                                numSets)
WrtTable(dt.windowReviewSummaryTable, outputFolderPath,
         "windowReviewSummaryTable")
rm(dt.windowReviewSummaryTable)

# create combined image arrays -----------------------------------------------------

combineSetArrays <- function(outputFolderPath, arrayName, numSets) {

  arrayList <- list()

  # append remaining arrays
  for (i.arrayNum in 1:numSets) {

    filePath <- paste0(outputFolderPath,"/", arrayName, i.arrayNum, ".rds")

    arrayList[[i.arrayNum]] <- readRDS(file = filePath)

    file.remove(filePath)

  }

  array <- abind(arrayList, along = 1, force.array = TRUE)

  outputFilePath <- paste0(outputFolderPath,"/",arrayName,".rds")

  saveRDS(array, outputFilePath)

  return(NULL)

}

combineSetArrays(outputFolderPath,"TestData", numSets)
combineSetArrays(outputFolderPath,"ValidationData", numSets)
combineSetArrays(outputFolderPath,"TrainingData", numSets)

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