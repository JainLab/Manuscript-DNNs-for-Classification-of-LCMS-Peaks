# Description ------------------------------------

# creates a csv file for tracking peak review

# User Input ------------------------------------------------------------------

# Path to folder with peak image files
# for Windows operating system use \\ or / instead of \ to separate directories
# example:
# peakTableParentDir <- 'M:\\folder1\\foldwerWithFile'
imageDir <- "M:\\folder1\\foldwerWithFile"

# Path to folder were file will be created
outputDir <- "M:\\folder2\\foldwerWithFile"

# shuffle the output table?
shuffle <- TRUE

# End user input ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

# Library (packages)--------------------------------------------------------
library(data.table) # faster version of data frames

# Timing --------------------------------------------------------------
startTime <- Sys.time()
startTimeStamp <- format(startTime, "(%Y.%m.%d.%H%M)")

# create table --------------------------------------------------------------


# create a vector images
vect.imagePaths <- list.files(path = imageDir,
                              pattern = "\\.jpg",
                              full.names = TRUE)


# create a vector images
vect.imageNames <- list.files(path = imageDir,
                              pattern = "\\.jpg",
                              full.names = FALSE)


# number of images
numImages <- length(vect.imageNames)

# create table
dt.ReviewTrack <- data.table(rowNum = seq(from = 1, to = numImages),
                             image = vect.imageNames,
                             viewed = numeric(length = numImages),
                             good1bad0 = rep(1, times = numImages),
                             windowBad = rep(0, times = numImages),
                             borderline = rep(0, times = numImages),
                             path = vect.imagePaths)

# shuffle -------------------------------------------------------------------

if (shuffle == TRUE) {

  # create shuffle column
  dt.ReviewTrack[, shuffNum := sample(nrow(dt.ReviewTrack), replace = FALSE)]

  # order by shuffled column
  setorder(dt.ReviewTrack, shuffNum)

  # replace rowNum col with shuffle column
  dt.ReviewTrack[, rowNum := shuffNum]

  # remove shuffle column
  dt.ReviewTrack$shuffNum <- NULL

}

# write table -----------------------------------------------------

# file path to write to
ouputFilePath <- paste0(outputDir,"/","reviewTracking",startTimeStamp,".csv")

# write file
fwrite(dt.ReviewTrack,
       file = ouputFilePath)

# Endscript (record keeping) -----------------------------------------

# get script file name
scriptName <- basename(sys.frame(1)$ofile)

# remove the file extension ".R"
scriptName <- gsub("\\.R","",scriptName)

if (1 == 1) {
  print("script complete")
  print(Sys.time() - startTime)
}