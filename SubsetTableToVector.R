# Description -----------------------------------------------------------------

# Remove rows from one table based on a label column from another table

# User Input ------------------------------------------------------------------

# Path to folder with individual peak table directories
# for Windows operating system use \\ or / instead of \ to separate directories
# example:

# path to folder where output directory will be created
# put copies of the CSV files in this folder
outPutParentDir <- "E:\\Example1\\Example1"

# CSV file with column of subset list
CSV_subsetList <- "exampleFileName1.csv"

# column index number for filter vector
colIndexFilterVect <- 1

# Directory of table file to be subset
inputCSVDir <- "E:\\Example2\\Example2"

# CSV file of table to be subset
CSV_toSubset <- "exampleFileName2.csv"

# column index number corresponding to filter vector
colIndexFullVect <- 1

# End user input ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

# advanced user input (probably don't need to change) ------------------------

# Library (packages)--------------------------------------------------------
library(data.table) # faster version of data frames

library(lcmsMetab)

# Timing --------------------------------------------------------------
startTime <- Sys.time()
startTimeStamp <- format(startTime, "(%Y%m%dh%H%M)")

# set or create the output folder -------------------------------------------

OutputFolderName <- "SubsetTableToVector"

# create output folder if it doesn't exist
outputFolderPath <- paste0(outPutParentDir,"/",
                           OutputFolderName,startTimeStamp)

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

# import data -----------------------------------------------------------------

# import table to be subset
importFilePath <- paste0(inputCSVDir, "/", CSV_toSubset)
dt.main <- fread(file = importFilePath)

# import vector to subset to
importFilePath <- paste0(outPutParentDir, "/", CSV_subsetList)
dt.subsetList <- fread(file = importFilePath,
                 select = c(colIndexFilterVect),
                 col.names = c("filter"))

# subset ----------------------------------------------------------------------

# rename column for ease
saveName <- colnames(dt.main)[colIndexFullVect]
colnames(dt.main)[colIndexFullVect] <- "fullVect"

dt.main <- dt.main[fullVect %in% dt.subsetList$filter]

# revert to original column name
colnames(dt.main)[colIndexFullVect] <- saveName

# write subset table ---------------------------------------------------------

exportName <- paste0(gsub("\\.csv$", "",CSV_toSubset),"_subset",startTimeStamp)

WrtTable(dt.main, outputFolderPath, gsub("\\.csv$", "",exportName))

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