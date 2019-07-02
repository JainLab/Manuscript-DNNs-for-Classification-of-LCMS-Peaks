# Description -------------------------------

# this script takes peak shape information exported by MZmine
# into separate csv files and combines the files into a single
# table and then uses the parameters and the category
# (good or bad peak) and uses multiple logistic
# regeression to generate a scoring function

# User Input ########################################################################

# Path to folder different peak parameter files file - adjusted file will be
# output to a new folder created in this folder
# use \\ instead of \ to separate directories
# example:
# 'M:\\folder1\\folderWithFile'
parentDir <- "M:\\folder1\\folderWithFile"

# peak attribute files are created by exporting an MZmine peak list to csv

# for the common elements select:
# row m/z
# row retention time

# for the data file elements select:
# peak duration
# peak height
# peak area
# peak FWHM
# peak tailing factor
# peak asymmetry factor

# make sure no other check boxes are checked as the
# column order must match the above order for all
# files


# Inputs for model building ---------------------------------

# Build new model? (TRUE or FALSE)
buildModel <- FALSE

# peak characteristic filenames for generating
# the multiple logistic regression model

# peak characteristic file of peaks reviewed
# for training model
# Injection columns are ordered:
# duration, height, area, FWHM, tailing factor,
# asymmetry factor. The MZmine export order

# max duration (min) as set in deconvolution settings
# used to adjust erroneous FWHM outputs
maxDurSet <- 4

# 3 separate datasets are utilized. One for training the
# model (Training), one to check for overfitting and revising
# the model (validation), and one for evaluating the model (Test)
# The script combines all of the input data and partitions
# it into these sets based on the following proportions:
# .5 training : .25 validation : .25 test

# Each data set requires a charachteristic csv file for all
# peaks reviewed with the columns indicated above.
# The script uses a table in csv format with peaks that
# were retained as good out of those that were reviewed
# column 1 mz, column 2 RT
# These tables are input in the form of a csv table
# formatted as follows:
# column 1: the peak characteristic file name
# column 2: the retained good peaks mz-rt file name
# this table should have a header row
# input the name of this file below:
csvTableOfDataFiles <- "ModelBuildData.csv"

# generate sample size-performance plots?
# this requires building and testing the models
# multiple times with increasing proportions of
# the training set to help identify if additional
# data is likely to significantly improve the model
sampSizePerfEval <- FALSE

# if evaluating sample size, input the number of
# sample size iterations. Each iteration will
# add 1/n (entered value) observations to the
# training set from the total training set
numSampSizeIter <- NA

# Reuse previous partition - for reproducing output
# partition will be reused rather than determined
# from new shuffling
usePrevPartition <- TRUE

# if reusing existing partition - input CSV
# this is the file output as trainValidTest.csv
prevPartitionCSV <- "trainValidTest.csv"

# copy input data files to output folder?
# For record keeping - copies files from
# csvTableOfDataFiles into the output folder.
copyInputDataFiles <- FALSE

# Inputs for model application -----------------------------

# apply model? (TRUE/FALSE)
applyModel <- TRUE

# model file to apply
# if model is being built, this should be set to NA
loadModelFile <- "PeakPredictionModel.rda"

# Files for features to run model on to predict
# peak classification
# Injection columns are ordered:
# duration, height, area, FWHM, tailing factor,
# asymmetry factor. The MZmine export order

# peak characteristic file for prediction data
pkCharFilePred <- "MzmineExport_6CharForPrediction_Subset.csv"

# MZmine database output:
# The script will output two MZmine databases
# one with all peaks above the input probability
# cutoff and one with  user set number of rows
# This is faster than having to rerun the script
# MZmine databases can also be created from the
# output feature summary table

# probability cutoff
probCutoff <- 0

# peak number cutoff
peakNumCutoff <- 1

# End User Input Section ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

# Advanced user input (probably don't need to change) ------------------------------------------------

# Index number of mz column
mzColNum <- 1

# Index number of RT column
rtColNum <- 2

# Index number of first injection column
firstInjColNum <- 3

# max p value of added variable
maxAllowAblePVal <- 0.05

# Library (packages)--------------------------------------------------------
library(data.table) # faster version of data frames
library(ggplot2) # for plotting - included in tidyverse
library(gridExtra) # for arranging multiple plot outputs

library(randomForest) # for random forest model
library(pROC) # for receiver operator curves
library(lcmsMetab) # LCMS data processing tools

# script version ------------------------------------------------------

scriptVersion <- "PeakPredictor(2019.07.01)"

# Timing --------------------------------------------------------------
startTime <- Sys.time()
startTimeStamp <- format(startTime, "(%Y.%m.%d.%H%M)")

# create the output folder -------------------------------------------

OutputFolderName <- "PeakPredictor"

outputFolderPath <- CreateOutputDir(parentDir, OutputFolderName)

# Initialize summary file -------------------------------------------

logFilePath <- paste0(outputFolderPath,"/",
                      "PeakPredictorLogFile",
                      startTimeStamp, ".txt")
logText <- paste0("start time ", startTimeStamp)
logText <- UpdateLogText(logText,"Initialize Summary File")
UpdateLogFile(logFilePath, logText)


# input checks -------------------------------------------------

# built model and model building data should not both be included
if (buildModel == TRUE & !is.na(loadModelFile)) {

  stop("If \"buildModel\" is set to TRUE, \"loadModelFile\' should be set to NA")

}

# summarize inputs ------------------------------------------------

logText <-
  UpdateLogText(logText,"working Directory", parentDir,
                "mz column number", mzColNum,
                "RT column number", rtColNum,
                "first injection column number", firstInjColNum,
                "build new model?",buildModel,
                "maximum duration setting" , maxDurSet,
                "model file",loadModelFile,
                "peak characteristic prediction file", pkCharFilePred)
UpdateLogFile(logFilePath, logText)

# defined functions ------------------------------------------------

# define function for copying files to output folder
copyToOutput <- function(fileX, outputFolderPath) {
  outputFilePath <- paste0(outputFolderPath,"/",fileX)
  file.copy(from = fileX, to = outputFilePath)
}

# :: Build Model --------------------------------------------

if (buildModel) {

  if (usePrevPartition == FALSE) {

    # :_prepare data --------------------------------------------

    # use CreateBuildData function from lcmsMetab package
    dt.featSummary.build <- CreateBuildData(parentDir, csvTableOfDataFiles,
                                            outputFolderPath)

    # export dt.featSummary.build --------------------------------

    WrtTable(dt.featSummary.build, outputFolderPath, "featureSummary")

    # downsample, shuffle and partition build data ------------------------------------

    # split input data into good and bad peaks
    dt.featSummary.good <- subset(dt.featSummary.build,
                                  numClass == 1)
    dt.featSummary.bad <- subset(dt.featSummary.build,
                                 numClass == 0)

    rm(dt.featSummary.build)

    numRowMinorityClass <- min(nrow(dt.featSummary.good),
                               nrow(dt.featSummary.bad))

    if (numRowMinorityClass < 70) {
      stop("minority class has too few observations - review more cases. Recommend more than 100 in each class")
    }

    # determine number of rows of each class in each set
    numRowTrain <- ceiling(numRowMinorityClass * 0.50)
    numRowValidate <- ceiling(numRowMinorityClass * 0.25)
    numRowTest <- numRowMinorityClass - numRowTrain - numRowValidate

    # create random vectors
    vect.shuffle.good <- sample(nrow(dt.featSummary.good))
    vect.shuffle.bad <- sample(nrow(dt.featSummary.bad))

    # add shuffle vectors to tables and sort by it
    dt.featSummary.good[, shuffle := vect.shuffle.good]
    setkey(dt.featSummary.good, shuffle)

    dt.featSummary.bad[, shuffle := vect.shuffle.bad]
    setkey(dt.featSummary.bad, shuffle)

    # trim tables to same number of rows
    dt.featSummary.good <- dt.featSummary.good[1:numRowMinorityClass,]
    dt.featSummary.bad <- dt.featSummary.bad[1:numRowMinorityClass,]

    # create a vector for row use category 1 Train, 2 validate, 3 test
    vect.useCat <- c(rep(1, numRowTrain),
                     rep(2, numRowValidate),
                     rep(3, numRowTest))



    # add vect.useCat as column
    dt.featSummary.good[, useCat := vect.useCat]
    dt.featSummary.bad[, useCat := vect.useCat]

    # combine for export
    dt.featSummary.trainValidTest <- rbind(dt.featSummary.good, dt.featSummary.bad,
                                           use.names = TRUE, fill = FALSE,
                                           idcol = NULL)

    # export for record
    WrtTable(dt.featSummary.trainValidTest, outputFolderPath, "trainValidTest")
  } # end if statement for no previous partition

  # import previously partitioned data
  if (usePrevPartition == TRUE) {

    # Full path to previously partitioned data
    prevPartFullPath <- paste0(parentDir,"/",prevPartitionCSV)

    # read file
    dt.featSummary.trainValidTest <- fread(prevPartFullPath)

    # separate back into good and bad tables
    dt.featSummary.good <- subset(dt.featSummary.trainValidTest,
                                  numClass == 1)
    dt.featSummary.bad <- subset(dt.featSummary.trainValidTest,
                                  numClass == 0)


  }

  rm(dt.featSummary.trainValidTest)

  # remove shuffle column
  dt.featSummary.good[, shuffle := NULL]
  dt.featSummary.bad[, shuffle := NULL]

  # create subset tables based on use category
  dt.featSummary.train <- rbind(subset(dt.featSummary.good, useCat == 1),
                                subset(dt.featSummary.bad, useCat == 1))
  dt.featSummary.valid <- rbind(subset(dt.featSummary.good, useCat == 2),
                                subset(dt.featSummary.bad, useCat == 2))
  dt.featSummary.test <- rbind(subset(dt.featSummary.good, useCat == 3),
                               subset(dt.featSummary.bad, useCat == 3))

  # remove useCat column
  dt.featSummary.train[, useCat := NULL]
  dt.featSummary.valid[, useCat := NULL]
  dt.featSummary.test[, useCat := NULL]

  # :_review variables ----------------------------------------------

  # variable correlations (moved for data combination) -----------------------

  # create variable correlation matrix
  dt.corMat <- cor(dt.featSummary.train[,4:(ncol(dt.featSummary.train) - 2)])
  dt.corMat <- data.table(dt.corMat, keep.rownames = TRUE)

  # convert to long form
  dt.cor.long <- melt(dt.corMat, id.vars = "rn",
                      measure.vars = c(2:ncol(dt.corMat)),
                      variable.name = "corVar",
                      value.name = "corr",
                      na.rm = TRUE)

  # remove correlation of variables with themselves
  dt.cor.long <- subset(dt.cor.long,
                        rn != corVar)

  # add column for absolute value of correlation
  dt.cor.long[, absCorr := abs(corr)]

  # export using defined function
  WrtTable(dt.cor.long, outputFolderPath, "VariableCorrelations")

  logText <-
    UpdateLogText(logText,"export variable correlations")
  UpdateLogFile(logFilePath, logText)

  # review individual variables ---------------------------------

  # number of potential predictive variables
  # total columns minus non-prediction variables
  numPotPred <- ncol(dt.featSummary.train) - 5

  # create tracking vectors
  vect.varPVal <- numeric(length = numPotPred)
  vect.relAbsEff <- numeric(length = numPotPred)
  vect.mostCorVar <- character(length = numPotPred)
  vect.mostCor <- numeric(length = numPotPred)
  vect.varSd <- numeric(length = numPotPred)
  vect.VarAbsEffSz <- numeric(length = numPotPred)
  vect.avgAbsEr <- numeric(length = numPotPred)
  vect.missIDrate <- numeric(length = numPotPred)
  vect.avgAbsErrKFold <- numeric(length = numPotPred)

  # create list of empty lists for each page
  plotPage <- rep(list(list()), numPotPred)

  # formula start text
  lgRegFmla <- "class  ~ "

  # store current warning setting
  setWarn <- options()[["warn"]]

  # set warn option to 0 setting
  options(warn = 0)

  # for each variable column, get a linear model p-value
  # and generate a violin plot
  for (varNum in 1:numPotPred) {

    # variable column
    lpvar.colNum <- varNum + 3

    # variable column name
    lpvar.varName <- colnames(dt.featSummary.train)[lpvar.colNum]

    # extract variable column
    lpvar.varCol <- subset(dt.featSummary.train,
                           select = c(lpvar.colNum))

    # variable min
    lpvar.varMin <- min(lpvar.varCol)

    # variable max
    lpvar.varMax <- max(lpvar.varCol)

    # create table from variable and class columns
    lptbl.varCol <- subset(dt.featSummary.train,
                           select = c(lpvar.colNum,(ncol(dt.featSummary.train) - 1), ncol(dt.featSummary.train)))

    # rename variable column for easy use
    colnames(lptbl.varCol)[1] <- "var"

    # variable standard deviation
    lpvar.sd <- sd(lptbl.varCol$var, na.rm = TRUE)

    # table of means by group
    lptbl.grpMeans <- lptbl.varCol[, mean(var, na.rm = TRUE),
                                          keyby = .(class)]
    colnames(lptbl.grpMeans)[2] <- "avgVar"

    # absolute effect size
    lpvar.absEffsz <- abs(lptbl.grpMeans$avgVar[1] - lptbl.grpMeans$avgVar[2])

    # relative absolute effect size
    lpvar.relAbsEffsz <- lpvar.absEffsz / lpvar.sd

    # subset correlation table
    lptbl.cor <- subset(dt.cor.long,
                        rn == lpvar.varName)

    # subset to max absolute correlation
    lptbl.cor <- subset(lptbl.cor,
                        absCorr == max(absCorr))

    # get most correlated variable
    lpvar.mostCorVar <- as.character(lptbl.cor$corVar)

    # get correlation for most correlated
    lpvar.mostCor <- round(lptbl.cor$corr,4)
    lpvar.mostCorPrnt <- round(lpvar.mostCor,4)

    # text for formula
    lpvar.formText <- paste0(lgRegFmla,lpvar.varName)

    # create formula object
    lpvar.fmla <- as.formula(lpvar.formText)

    # create logistic regression model
    lpvar.lm <- glm(lpvar.fmla,
                    family = binomial,
                    data = dt.featSummary.train)

    # use predict to get a vector of the predicted probablilities of good
    lpvect.predCheck <- predict(lpvar.lm,
                              newdata = dt.featSummary.train,
                              type = "response")

    # add predicted probability
    lptbl.varCol[, predProb := lpvect.predCheck]

    # create column for error
    lptbl.varCol[, error := numClass - predProb]

    # create column for absolute error
    lptbl.varCol[, absError := abs(error)]

    # create column for misidentification (absolute error > 0.5)
    lptbl.varCol[, missID := ifelse(absError > 0.5, 1, 0)]

    # find average absolute error
    lpvar.avgAbsEr <- mean(lptbl.varCol$absError)

    # calculate missidentification rate
    lpvar.missIDrate <- mean(lptbl.varCol$missID)

    # get p-value
    lpvar.pval <- coef(summary(lpvar.lm))[2,4]
    # format for printing
    lpvar.pvalprnt <- sprintf("%1.2E",lpvar.pval)

    # k-fold cross validation
    # number of folds
    k <- 5
    # number of training observations
    n.trainObs <- nrow(dt.featSummary.train)
    # vector to randomly shuffle training observations
    vect.random.order <- sample(seq(1,n.trainObs),n.trainObs, replace = FALSE)
    # vector of sizes of cross-validation folds
    cv.size <- c(rep(floor(n.trainObs/k),(k - 1)),
                 (n.trainObs - sum(rep(floor(n.trainObs/k),(k - 1)))))
    # vector for starting points of folds
    cum.size <- c(0, cumsum(cv.size[-length(cv.size)]))
    # initialize vector for tracking fold mean absolute error
    lpvect.FoldavgAbsErr <- rep(0,k)

    for (i in 1:k) {
      # vector of indices of the random indices to be included in validation
      lpvect.foldInds <- (cum.size[i] + 1):(cum.size[i] + cv.size[i])
      # train on subset that does not include selected observations
      lpdt.train.subset <- dt.featSummary.train[-vect.random.order[lpvect.foldInds],]
      # create multiple logistic regression model
      lplm.kFold <- glm(lpvar.fmla,
                        family = binomial,
                        data = lpdt.train.subset)
      # create validation subset
      lpdt.FoldValidate <- dt.featSummary.train[vect.random.order[lpvect.foldInds],]
      # predict on validation subset
      predict.validateFold <- predict(lplm.kFold,
                                      newdata = lpdt.FoldValidate,
                                      type = "response")
      # vector of squared errors
      lpvect.FoldavgAbsErr[i] <- mean(abs(lpdt.FoldValidate$numClass - predict.validateFold))
    }

    # average absolute error for all k folds
    lpvar.avgAbsErrKFold <- mean(lpvect.FoldavgAbsErr)

    # add to tracking vectors
    vect.mostCorVar[varNum] <- lpvar.mostCorVar
    vect.mostCor[varNum] <- lpvar.mostCor
    vect.varSd[varNum] <- lpvar.sd
    vect.VarAbsEffSz[varNum] <- lpvar.absEffsz
    vect.relAbsEff[varNum] <- lpvar.relAbsEffsz
    vect.varPVal[varNum] <- lpvar.pval
    vect.avgAbsEr[varNum] <- lpvar.avgAbsEr
    vect.missIDrate[varNum] <- lpvar.missIDrate
    vect.avgAbsErrKFold[varNum] <- lpvar.avgAbsErrKFold

    # create plot title
    lpvar.infoTitle <- paste0(" P-value = ", lpvar.pvalprnt,
           "\n Rel Abs Effect = ", round(lpvar.relAbsEffsz, digits = 3),
           "\n Most correlated: ", lpvar.mostCorVar,
           " ", lpvar.mostCorPrnt,
           "\n AvgAbsError: ", round(lpvar.avgAbsEr, digits = 3))

    # violin plots
    plotPage[[varNum]][[1]] <-
      ggplot() +
      geom_violin(data = dt.featSummary.train,
                  aes_string(x = "class", y = lpvar.varName),
                  fill = NA, size = .75, draw_quantiles = c(0.5)) +
      # geom_point(data = dt.featSummary.train, aes_string(x = "class", y = lpvar.varName), shape = ".") +
      ggtitle(lpvar.varName) +
      ylab(lpvar.varName) +
      xlab('Classification') +
      ggThemeLcmsMetab()

    # Logistic Regression Plots

    # prediction values for creating curve
    lptbl.plotPred <- data.table(Xpred = seq(from = lpvar.varMin,
                                             to = lpvar.varMax,
                                             length.out = 200))
    colnames(lptbl.plotPred)[1] <- lpvar.varName

    # prediction values vectors
    lpvect.plotYpred <- predict(lpvar.lm,
                                newdata = lptbl.plotPred,
                                type = "response")

    # add to plot data table
    lptbl.plotPred[, probPred := lpvect.plotYpred]

    # logistic regression plots
    plotPage[[varNum]][[2]] <-
      ggplot() +
      # geom_point(data = dt.featSummary.train,
      # aes_string(x = lpvar.varName, y = "numClass"), shape = ".") +
      geom_line(data = lptbl.plotPred,
                aes_string(x = lpvar.varName, y = "probPred")) +
      geom_hline(yintercept = 1) +
      geom_hline(yintercept = 0) +
      ggtitle(lpvar.infoTitle) +
      ylab("predicted probability") +
      xlab(lpvar.varName) +
      ggThemeLcmsMetab()

  } # end variable examination for loop with plots

  # return warning option
  options(warn = setWarn)

  # create data table for variable summary
  dt.varSummary <- data.table(var = colnames(dt.featSummary.train)[4:(3 + numPotPred)],
                              mostCorVar = vect.mostCorVar,
                              mostCor = vect.mostCor,
                              absMostCor = abs(vect.mostCor),
                              StdDev = vect.varSd,
                              absEffSz =  vect.VarAbsEffSz,
                              relAbsEff = vect.relAbsEff,
                              pVal_1Regr = vect.varPVal,
                              avgAbsEr = vect.avgAbsEr,
                              missIDrate = vect.missIDrate,
                              avgAbsErrKFold = vect.avgAbsErrKFold)

  WrtTable(table = dt.varSummary, outputFolderPath,
           name = paste0("variableSummary", startTimeStamp))

  logText <-
    UpdateLogText(logText,"individual variable review completed")
  UpdateLogFile(logFilePath, logText)

  # PDF individual variable plots --------------------------------------------------

  PdfFilePath <- paste0(outputFolderPath,"/",
                        "IndividualVariablePlots",
                        startTimeStamp,".pdf")

  pdf(PdfFilePath, width = 10, height = 7.5, paper = "USr", onefile = TRUE)
  for (i in seq(length(plotPage))) {
    grid.arrange(grobs = plotPage[[i]], ncol = 2)

    percComplete <- i/length(plotPage) * 100

    print(paste0("individual variable review pdf completion ",
                 round(percComplete,1),"%"))
  }
  dev.off()

  logText <-
    UpdateLogText(logText,"individual variable review PDF file generated")
  UpdateLogFile(logFilePath, logText)

  # plot avgAbsError & missIDrate Individual variables ---------------------------------------------------------

  # add kFoldAvgAbsErRank as additional column
  dt.varSummary[, kFoldAvgAbsErRank := frank(-avgAbsErrKFold)]

  # sort by rank for plot labelling
  setkey(dt.varSummary, kFoldAvgAbsErRank)

  # convert avgAbsErRank to factor
  dt.varSummary$kFoldAvgAbsErRank <- as.factor(dt.varSummary$kFoldAvgAbsErRank)

  # save filenames in order as vector
  vect.varNames <- dt.varSummary$var

  # create empty plotlist
  plotVarModEr <- rep( list(list()), 3 )

  plotVarModEr[[1]][[1]] <-
    ggplot() +
    geom_point(data = dt.varSummary,
               aes(x = kFoldAvgAbsErRank, y = round(avgAbsErrKFold,3))) +
    ggtitle("Single Variable Model Error K-Fold Cross-Validation") +
    ylab("k-fold average absolute error") +
    xlab("prediction variable") +
    scale_x_discrete(labels = vect.varNames) +
    ggThemeLcmsMetab() +
    theme(axis.text.x = element_text(angle = 90))

  # create plot object average absolute error
  plotVarModEr[[2]][[1]] <-
    ggplot() +
    geom_point(data = dt.varSummary,
               aes(x = kFoldAvgAbsErRank, y = round(avgAbsEr,3))) +
    ggtitle("Single Variable Model Error") +
    ylab("average absolute error") +
    xlab("prediction variable") +
    scale_x_discrete(labels = vect.varNames) +
    ggThemeLcmsMetab() +
    theme(axis.text.x = element_text(angle = 90))

  # create plot object average absolute error
  plotVarModEr[[3]][[1]] <-
    ggplot() +
    geom_point(data = dt.varSummary,
               aes(x = kFoldAvgAbsErRank, y = round(missIDrate,3))) +
    ggtitle("Single Variable Model Misidentification Rate") +
    ylab("misidentification rate") +
    xlab("prediction variable") +
    scale_x_discrete(labels = vect.varNames) +
    ggThemeLcmsMetab() +
    theme(axis.text.x = element_text(angle = 90))

  # create pdf
  PdfPlotsSingleFile(plotVarModEr, "landscape", outputFolderPath, "PlotIndividualVarAvgAbsErMissID")

  logText <-
    UpdateLogText(logText,"plot complete: individual avg abs error and missidentification rate")
  UpdateLogFile(logFilePath, logText)

  # :_Build Models ---------------------------------------------------
  #
  # Setup Sample Size check -------------------------------------------

  if (sampSizePerfEval == TRUE) {
    vect.trainSize <- c(seq(from = ceiling(nrow(dt.featSummary.train)/numSampSizeIter),
                            to = nrow(dt.featSummary.train),
                            by = ceiling(nrow(dt.featSummary.train)/numSampSizeIter)),
                        nrow(dt.featSummary.train))

  } else {
    vect.trainSize <- seq(from = nrow(dt.featSummary.train),
                          to = nrow(dt.featSummary.train))
  }

  vect.acc.mlr <- numeric(length = length(vect.trainSize))
  vect.acc.rndFrst <- numeric(length = length(vect.trainSize))

  for (trainSize.i in 1:length(vect.trainSize)) {

    if (sampSizePerfEval == FALSE) {
      lpdt.featSummary.train <- copy(dt.featSummary.train)
    } else {

      trainSize <- vect.trainSize[trainSize.i]

      # create shuffle vector - duplicated to catch equal true and false peaks
      vect.shuffle.train <- sample(nrow(dt.featSummary.train)/2,
                                   size = nrow(dt.featSummary.train)/2,
                                   replace = FALSE)

      vect.shuffle.train <- rep(vect.shuffle.train, 2)

      # add shuffle vectors to tables and sort by it
      lpdt.featSummary.train <- copy(dt.featSummary.train)
      lpdt.featSummary.train[, shuffle := vect.shuffle.train]
      setkey(lpdt.featSummary.train, shuffle)

      # reduce to size inicated in vector
      lpdt.featSummary.train <- lpdt.featSummary.train[1:trainSize,]

      vect.trainSize[trainSize.i] <- nrow(lpdt.featSummary.train)

      # remove shuffle column
      lpdt.featSummary.train$shuffle <- NULL

    }


    # _multiple logistic regression model ----------------------------------------------

    # model curation training set (forward selection) ------------------------------------------------------------

    # set var as character
    dt.varSummary$var <- as.character(dt.varSummary$var)

    # find the variable that produces lowest average absolute error

    # find the minimum avgAbsErrKFold for all single variable models
    minAbsErKFold <- min(dt.varSummary$avgAbsErrKFold, na.rm = TRUE)

    # subset dt.varSummary to the variable that produced the
    # lowest absolute error
    dt.lowAbsEr <- subset(dt.varSummary, avgAbsErrKFold == minAbsErKFold)

    # create a data table of the remaining variables
    # to shrink as new variables are added to model
    dt.hiAbsErShrink <- subset(dt.varSummary, avgAbsErrKFold != minAbsErKFold,
                               select = c("var"))

    # create new vectors for tracking improvement of model
    # with each variable addition
    vect.varImprv <- character(length = numPotPred)
    vect.modelKFoldAvgAbsEr <- numeric(length = numPotPred)
    vect.modelMissIDrate <- numeric(length = numPotPred)
    vect.modelErPercImp <- numeric(length = numPotPred)
    vect.modelPvalImp <- numeric(length = numPotPred)
    vect.modelPvalMaxWald <- numeric(length = numPotPred)
    vect.modelIncludeVar <- logical(length = numPotPred)

    # add best single model variable as first entry in vectors
    vect.varImprv[1] <- dt.lowAbsEr$var[1]
    vect.modelKFoldAvgAbsEr[1] <- minAbsErKFold
    vect.modelMissIDrate[1] <- dt.lowAbsEr$missIDrate[1]
    vect.modelErPercImp[1] <- NA
    vect.modelPvalImp[1] <- dt.lowAbsEr$pVal_1Regr[1]
    vect.modelPvalMaxWald[1] <- dt.lowAbsEr$pVal_1Regr[1]
    vect.modelIncludeVar[1] <- TRUE

    # create starting model formula
    lpvar.formText <- paste0(lgRegFmla,dt.lowAbsEr$var[1])

    # initialize list of formulas and list of models
    list.fmla <- list()
    list.models <- list()

    # initialize lpvar.minkFoldAvgAbsErAdd
    lpvar.minkFoldAvgAbsErAdd <- minAbsErKFold

    # for loop to determine variable order of most effect (forward selection)
    for (varNum in 1:(numPotPred - 1)) {

      # get the number of remaining variables
      varRemaining <-  numPotPred - varNum

      # create tracking vectors for checking the addition
      # of each variable
      lpvect.varImprvChk <- character(length = varRemaining)
      lpvect.modelAvgAbsErChk <- numeric(length = varRemaining)
      lpvect.modelMissIDRateChk <- numeric(length = varRemaining)
      lpvect.pval.addVar <- numeric(length = varRemaining)
      lpvect.pval.maxWald <- numeric(length = varRemaining)
      lpvect.kFoldAvgAbsEr <- numeric(length = varRemaining)

      # create formula object
      lpvar.fmla <- as.formula(lpvar.formText)

      # add to list
      list.fmla[[varNum]] <- lpvar.fmla

      # create logistic regression model
      lpvar.lm <- glm(lpvar.fmla,
                      family = binomial,
                      data = lpdt.featSummary.train)

      # add to list
      list.models[[varNum]] <- lpvar.lm

      # store deviance and df.residual of existing model
      lpvar.deviance.exist <- lpvar.lm$deviance
      lpvar.dfRes.exist <- lpvar.lm$df.residual

      # make a model adding each remaining variable to
      # previous optimum variables
      for (varCheckNum in 1:varRemaining) {

        # get name of variable to test
        lpvar.varName <- dt.hiAbsErShrink$var[varCheckNum]

        # add to tracking vector
        lpvect.varImprvChk[varCheckNum] <- lpvar.varName

        # create new formula to test
        lpvar.formTextChk <- paste0(lpvar.formText," + ",lpvar.varName)

        # create formula object
        lpvar.fmla <- as.formula(lpvar.formTextChk)

        # create logistic regression model
        lpvar.lm <- glm(lpvar.fmla,
                        family = binomial,
                        data = lpdt.featSummary.train)

        # create logistic regression model summary object
        lpvar.lm.smry <- summary(lpvar.lm)

        # get max wald p-value of model
        lpvar.maxpvalWald <- max(coef(lpvar.lm.smry)[,4])

        # log ratio test p value (compare with existing simpler model)
        lpvar.pvalLRT <- with(lpvar.lm, pchisq(lpvar.deviance.exist - deviance,
                                               lpvar.dfRes.exist - df.residual,
                                               lower.tail = FALSE))

        # use predict to get a vector of the predicted probablilities of good
        lpvect.predCheck <- predict(lpvar.lm,
                                    newdata = lpdt.featSummary.train,
                                    type = "response")

        # add predicted probability
        lptbl.varCol[, predProb := lpvect.predCheck]

        # create column for error
        lptbl.varCol[, error := numClass - predProb]

        # create column for absolute error
        lptbl.varCol[, absError := abs(error)]

        # find average absolute error
        lpvar.avgAbsEr <- mean(lptbl.varCol$absError)

        # create column for misidentification (absolute error > 0.5)
        lptbl.varCol[, missID := ifelse(absError > 0.5, 1, 0)]

        # calculate missidentification rate
        lpvar.missIDrate <- mean(lptbl.varCol$missID)

        # k-fold cross validation
        # number of folds
        k <- 5
        # number of training observations
        n.trainObs <- nrow(lpdt.featSummary.train)
        # vector to randomly shuffle training observations
        vect.random.order <- sample(seq(1,n.trainObs),
                                    n.trainObs,
                                    replace = FALSE)
        # vector of sizes of cross-validation folds
        cv.size <- c(rep(floor(n.trainObs/k),(k - 1)),
                     (n.trainObs - sum(rep(floor(n.trainObs/k),(k - 1)))))
        # vector for starting points of folds
        cum.size <- c(0, cumsum(cv.size[-length(cv.size)]))
        # initialize vector for tracking fold mean absolute error
        lpvect.FoldavgAbsErr <- rep(0,k)

        for (i in 1:k) {
          # vector of indices of the random indices to be included in validation
          lpvect.foldInds <- (cum.size[i] + 1):(cum.size[i] + cv.size[i])
          # train on subset that does not include selected observations
          lpdt.train.subset <- lpdt.featSummary.train[-vect.random.order[lpvect.foldInds],]
          # create multiple logistic regression model
          lplm.kFold <- glm(lpvar.fmla,
                            family = binomial,
                            data = lpdt.train.subset)
          # create validation subset
          lpdt.FoldValidate <- lpdt.featSummary.train[vect.random.order[lpvect.foldInds],]
          # predict on validation subset
          predict.validateFold <- predict(lplm.kFold,
                                          newdata = lpdt.FoldValidate,
                                          type = "response")
          # vector of squared errors
          lpvect.FoldavgAbsErr[i] <- mean(abs(lpdt.FoldValidate$numClass - predict.validateFold))
        }

        # average absolute error for all k folds
        lpvar.avgAbsErrKFold <- mean(lpvect.FoldavgAbsErr)

        # add to tracking vectors
        lpvect.modelAvgAbsErChk[varCheckNum] <- lpvar.avgAbsEr
        lpvect.pval.addVar[varCheckNum] <- lpvar.pvalLRT
        lpvect.pval.maxWald[varCheckNum] <- lpvar.maxpvalWald
        lpvect.modelMissIDRateChk[varCheckNum] <- lpvar.missIDrate
        lpvect.kFoldAvgAbsEr[varCheckNum] <- lpvar.avgAbsErrKFold

      } # end for-loop check remaining variables

      # create table from checked variable values
      lptbl.varImprvChk <- data.table(varImprv = lpvect.varImprvChk,
                                      avgAbsError = lpvect.modelAvgAbsErChk,
                                      missIDrate = lpvect.modelMissIDRateChk,
                                      pvalAddVarLRT = lpvect.pval.addVar,
                                      pvalMaxWald = lpvect.pval.maxWald,
                                      kFoldAvgAbsEr = lpvect.kFoldAvgAbsEr)

      # set previous minimum avgAbsEr to track improvement
      lpvar.minkFoldAvgAbsErAdd.old <- lpvar.minkFoldAvgAbsErAdd

      # find minimum avgAbsEr for added variable
      lpvar.minkFoldAvgAbsErAdd <- min(lptbl.varImprvChk$kFoldAvgAbsEr)


      # calculate percent reduction
      lpvar.percErReduction <-
        (lpvar.minkFoldAvgAbsErAdd.old - lpvar.minkFoldAvgAbsErAdd)/lpvar.minkFoldAvgAbsErAdd.old * 100


      # find corresponding variable and p-value of improvement
      lptbl.minErVarImprvChk <- subset(lptbl.varImprvChk,
                                       kFoldAvgAbsEr == lpvar.minkFoldAvgAbsErAdd)

      lpvar.varAdd <- lptbl.minErVarImprvChk$varImprv[1]
      lpvar.pvalAddVar <- lptbl.minErVarImprvChk$pvalAddVarLRT[1]
      lpvar.pvalMaxWald <- lptbl.minErVarImprvChk$pvalMaxWald[1]
      lpvar.missIDrateLowEr <- lptbl.minErVarImprvChk$missIDrate[1]

      # add variable to formula
      lpvar.formText <- paste0(lpvar.formText, " + ", lpvar.varAdd)

      # add to tracking vectors
      vect.varImprv[varNum + 1] <- lpvar.varAdd
      vect.modelKFoldAvgAbsEr[varNum + 1] <-  lpvar.minkFoldAvgAbsErAdd
      vect.modelMissIDrate[varNum + 1] <- lpvar.missIDrateLowEr
      vect.modelErPercImp[varNum + 1] <- lpvar.percErReduction
      vect.modelPvalImp[varNum + 1] <- lpvar.pvalAddVar
      vect.modelPvalMaxWald[varNum + 1] <- lpvar.pvalMaxWald

      # check if variable will be included in final model
      lpvar.include <- ifelse(max(lpvar.percErReduction <= 0, na.rm = TRUE),
                              FALSE, TRUE)

      # add to inclusion tracking vector
      vect.modelIncludeVar[varNum + 1] <- lpvar.include

      # remove variable from dt.hiAbsErShrink
      dt.hiAbsErShrink <- subset(dt.hiAbsErShrink,
                                 var != lpvar.varAdd)

      percComplete <- round(varNum/(numPotPred - 1) * 100, 1)

      print(paste0("model curation loop completion ",
                   round(percComplete, 1),"%"))

      # end loop if k-fold cross validation suggests no improvement
      if (lpvar.percErReduction <= 0) {

        break

      }

    } # end forward selection for loop

    # create table from tracking vectors
    dt.modelAddVar <- data.table(order = as.factor(seq(from = 1, to = numPotPred)),
                                 var =  vect.varImprv,
                                 kFoldAvgAbsError = vect.modelKFoldAvgAbsEr,
                                 missIDrate = vect.modelMissIDrate,
                                 ErPercImp = vect.modelErPercImp,
                                 pValAddVar = vect.modelPvalImp,
                                 pValMaxWald = vect.modelPvalMaxWald,
                                 include = vect.modelIncludeVar)

    # subset to cut unused rows
    dt.modelAddVar <- subset(dt.modelAddVar,
                             kFoldAvgAbsError != 0)

    # Validation and revision -------------------------------------

    # create vector of variables added to validate
    vect.varValidate <- subset(dt.modelAddVar,
                               include == TRUE,
                               select = c("var"))

    vect.varValidate <- unlist(vect.varValidate)

    # get number of models to validate
    numValidate <- length(vect.varValidate)

    # remove models that will not meeting training criteria
    # from list.models
    list.models <- list.models[1:numValidate]

    # define a function for collecting validation values for model
    modelValidate <- function(model, dt.validationData) {

      # model is the logistic regression model
      # dt.validationData is the validation data set

      # get number of observations in training data set
      numTrainObs <- nrow(model[["data"]])

      # get model formula as text string
      fmla <- Reduce(paste, deparse(model[["terms"]]))

      # get model coefficients
      coefficients <- paste(round(model[["coefficients"]],1), sep = ",", collapse = " ")

      # use predict to get a vector of the predicted probablilities of good
      vect.predCheck <- predict(model,
                                newdata = dt.validationData,
                                type = "response")

      # add predicted probability
      dt.validationData[, predProb := vect.predCheck]

      # create column for error
      dt.validationData[, error := numClass - predProb]

      # create column for absolute error
      dt.validationData[, absError := abs(error)]

      # create column for misidentification (absolute error > 0.5)
      dt.validationData[, missID := ifelse(absError > 0.5, 1, 0)]

      # find average absolute error
      avgAbsEr <- mean(dt.validationData$absError)

      # calculate missidentification rate
      missIDrate <- mean(dt.validationData$missID)

      # create data table of validation values
      dt.validationValues <- data.table(numObs = numTrainObs,
                                        fmla = fmla,
                                        coefficients = coefficients,
                                        avgAbsEr = avgAbsEr,
                                        missIDrate = missIDrate)

      return(dt.validationValues)

    }

    # apply function to list.models
    dt.validation <- lapply(list.models, modelValidate,
                            dt.featSummary.valid)

    # combine list of tables into one table
    dt.validation <- rbindlist(dt.validation)

    # add variable column
    dt.validation[, var := vect.varValidate]

    # create a vector of avgAbsError differences
    vect.difference <- diff(dt.validation$avgAbsEr)
    vect.difference <- c(NA, vect.difference)

    # create a vector of percent error improvement
    vect.percErImprv <- -vect.difference /
      (dt.validation$avgAbsEr - vect.difference) * 100

    # add as column
    dt.validation[, percErImprv := vect.percErImprv]

    # add order column
    dt.validation[, order := seq(nrow(dt.validation))]

    # find first added variable with no error improvement
    dt.noImprove <- subset(dt.validation,
                           percErImprv < 0)


    # variable number to include
    if (nrow(dt.noImprove) == 0) {

      numVarInclude <- nrow(dt.validation)

    } else {

      dt.noImprove <- subset(dt.noImprove,
                             order == min(dt.noImprove$order))

      numVarInclude <- dt.noImprove$order[1] - 1

    }

    # add inclusion column to table
    dt.validation[, include := c(rep(TRUE, numVarInclude),
                                 rep(FALSE, nrow(dt.validation) - numVarInclude))]

    # build final mult Logit Regr model -----------------------------------------

    # select model formula for model before validation showed
    # no improvement
    finalFmla <- list.fmla[[numVarInclude]]

    # trim added columns from dt.featSummary.valid
    dt.featSummary.valid <- subset(dt.featSummary.valid,
                                   select = colnames(lpdt.featSummary.train))

    # create model
    logRegModel.final <- glm(finalFmla, family = binomial,
                       data = lpdt.featSummary.train,
                       model = FALSE,
                       y = FALSE)

    # test model multiple logistic regression ---------------------------------------

    # use defined function to get test values
    dt.test <- modelValidate(logRegModel.final, dt.featSummary.test)

    # create table of probability values and class
    dt.MLogRegTest <- subset(dt.featSummary.test,
                             select = c("featID","numClass", "class"))

    # use predict to get a vector of the predicted probablilities of good
    vect.predCheck <- predict(logRegModel.final,
                              newdata = dt.featSummary.test,
                              type = "response")

    # add predicted probability
    dt.MLogRegTest[, predProb := vect.predCheck]

    # create column for error
    dt.MLogRegTest[, error := numClass - predProb]

    # create column for absolute error
    dt.MLogRegTest[, absError := abs(error)]

    # create column for squared error
    dt.MLogRegTest[, SqError := error^2]

    # find average squared error
    sumSqError <- sum(dt.MLogRegTest$SqError)
    avgSqError <- mean(dt.MLogRegTest$SqError)
    avgAbsError <- mean(dt.MLogRegTest$absError)
    sdAbsError <- sd(dt.MLogRegTest$absError)
    missIDrate <- round(dt.test$missIDrate[1],4)
    acc.mlr <- 1 - missIDrate

    # find margins of error (Moe)
    avgAbsErrorMoe <- 1.96 * sdAbsError / sqrt(nrow(dt.MLogRegTest))
    missIDrateMoe <- 1.96 *  sqrt(missIDrate * (1 - missIDrate) / nrow(dt.MLogRegTest))

    # calculate correlation
    evalCorr <- cor(x = dt.MLogRegTest$numClass, y = dt.MLogRegTest$predProb)

    # add accuracy to tracking vector
    vect.acc.mlr[trainSize.i] <- acc.mlr

    # create and write score tables (mult Logi Reg) --------------------------------
    vect.predCheck <- predict(logRegModel.final,
                              newdata = dt.featSummary.train,
                              type = "response")

    dt.MLogRegTrain <- data.table(featID = dt.featSummary.train$featID,
                                  class = dt.featSummary.train$class,
                                  predProb = vect.predCheck)

    vect.predCheck <- predict(logRegModel.final,
                              newdata = dt.featSummary.valid,
                              type = "response")

    dt.MLogRegValid <- data.table(featID = dt.featSummary.valid$featID,
                                  class = dt.featSummary.valid$class,
                                  predProb = vect.predCheck)

    WrtTable(dt.MLogRegTest, outputFolderPath, "ScoresTestMultLogReg")
    WrtTable(dt.MLogRegTrain, outputFolderPath, "ScoresTrainMultLogReg")
    WrtTable(dt.MLogRegValid, outputFolderPath, "ScoresValidationMultLogReg")

    # plot prob vs prob rank color on numClass? (too many points)

    # _random forest model -----------------------------------------

    # convert class to factor for random for random forest to
    # run in classification mode instead of regression mode
    lpdt.featSummary.train$classFact <- as.factor(lpdt.featSummary.train$class)
    dt.featSummary.valid$classFact <- as.factor(dt.featSummary.valid$class)
    dt.featSummary.test$classFact <- as.factor(dt.featSummary.test$class)

    # create new formula for random forest
    fmla.rndFr <- ""

    for (i.var in 1:(length(colnames(lpdt.featSummary.train)) - 6)) {
      lpvar <- colnames(lpdt.featSummary.train)[i.var + 3]

      if (i.var == 1) {
        fmla.rndFr <- paste0(lpvar)
      } else {
        fmla.rndFr <- paste0(fmla.rndFr, " + ", lpvar)
      }
    }
    fmla.rndFr <- paste0("classFact ~ ",fmla.rndFr)
    fmla.rndFr <- as.formula(fmla.rndFr)

    # Random forest calibrate (node size optimize) -----------------------------------------------

    maxNumSizeCheckIter <- 50
    numSizeCheckIter <- 1
    nodeSize <- ceiling(nrow(lpdt.featSummary.train)/2)
    vect.nodeSize <- c(nodeSize)

    # create vector of node sizes to test
    while (nodeSize > 1 & numSizeCheckIter < 50) {

      numSizeCheckIter <- numSizeCheckIter + 1

      nodeSize <- ceiling(nodeSize/2)

      vect.nodeSize[numSizeCheckIter] <- nodeSize

    } # end node size loop

    vect.rfAvgAbsErCalibrate <- numeric(length = length(vect.nodeSize))
    vect.rfAvgAbsErTrain <- numeric(length = length(vect.nodeSize))

    # test each node size in vector
    for (i_nodeSize in 1:length(vect.nodeSize)) {

      lpvar.nodeSize <- vect.nodeSize[i_nodeSize]

      # create model with training data
      lpmod.rndFrMod <- randomForest(fmla.rndFr,
                               data = lpdt.featSummary.train,
                               ntree = 1000,
                               nodesize = lpvar.nodeSize,
                               importance = TRUE)

      # test model with calibration set
      lpdt.rndFrstCalib <- as.data.table(predict(lpmod.rndFrMod,
                                              dt.featSummary.valid,
                                              type = "prob"))
      colnames(lpdt.rndFrstCalib) <- c("badProb", "goodProb")

      # add calibration set labels
      lpdt.rndFrstCalib[, Good1Bad0 := dt.featSummary.valid$numClass]

      # get absolute error
      lpdt.rndFrstCalib[, absError := abs(Good1Bad0 - goodProb)]

      # compute average absolute error
      lpvar.avgAbsErCalibrate <- mean(lpdt.rndFrstCalib$absError)

      # add to tracking vector
      vect.rfAvgAbsErCalibrate[i_nodeSize] <- lpvar.avgAbsErCalibrate

      # do the same for training set for comparison
      lpdt.rndFrstTrain <- as.data.table(predict(lpmod.rndFrMod,
                                                 lpdt.featSummary.train,
                                                 type = "prob"))
      colnames(lpdt.rndFrstTrain) <- c("badProb", "goodProb")
      lpdt.rndFrstTrain[, Good1Bad0 := lpdt.featSummary.train$numClass]
      lpdt.rndFrstTrain[, absError := abs(Good1Bad0 - goodProb)]
      lpvar.avgAbsErTrain <- mean(lpdt.rndFrstTrain$absError)
      vect.rfAvgAbsErTrain[i_nodeSize] <- lpvar.avgAbsErTrain

    } # end node size test loop

    dt.rfNodeSizeTest <- data.table(nodeSize = vect.nodeSize,
                                    avgAbsErTrain = vect.rfAvgAbsErTrain,
                                    avgAbsErCalibrate = vect.rfAvgAbsErCalibrate)

    # find optimum node size
    rfOptNodeSize <- subset(dt.rfNodeSizeTest,
                            avgAbsErCalibrate ==
                              min(dt.rfNodeSizeTest$avgAbsErCalibrate))
    rfOptNodeSize <- rfOptNodeSize$nodeSize[1]

    # random forest final model build -------------------------------------
    rndFrMod <- randomForest(fmla.rndFr,
                             data = lpdt.featSummary.train,
                             nodesize = rfOptNodeSize,
                             importance = TRUE)

    # test model: Random Forest -----------------------------------------

    dt.rndFrstTest <- as.data.table(predict(rndFrMod,
                                            dt.featSummary.test,
                                            type = "prob"))
    colnames(dt.rndFrstTest) <- c("Bad","Good")

    dt.rndFrstTest$predTest <- predict(rndFrMod, dt.featSummary.test, type = "class")

    # add observed class
    dt.rndFrstTest$obsClass <- dt.featSummary.test$classFact
    dt.rndFrstTest$featID <- dt.featSummary.test$featID

    # get model accuracy
    acc.rndFrst <- round(mean(dt.rndFrstTest$predTest == dt.featSummary.test$classFact), 4)
    missIDrate.rndFrst <- 1 - acc.rndFrst

    # add accuracy to tracking vector
    vect.acc.rndFrst[trainSize.i] <- acc.rndFrst

    print(paste0("sample size review loop ",
                 trainSize.i, " of ", length(vect.trainSize)))

  } # end sample size eval loop

  # Scores validation set
  dt.rndFrstValid <- as.data.table(predict(rndFrMod,
                                          dt.featSummary.valid,
                                          type = "prob"))
  colnames(dt.rndFrstValid) <- c("Bad","Good")
  # add observed class
  dt.rndFrstValid$obsClass <- dt.featSummary.valid$class
  dt.rndFrstValid$featID <- dt.featSummary.valid$featID

  # Scores training set
  dt.rndFrstTrain <- as.data.table(predict(rndFrMod,
                                           dt.featSummary.train,
                                           type = "prob"))
  colnames(dt.rndFrstTrain) <- c("Bad","Good")
  # add observed class
  dt.rndFrstTrain$obsClass <- dt.featSummary.train$class
  dt.rndFrstTrain$featID <- dt.featSummary.train$featID

  # Variable importance ---------------------------------------------

  importance(rndFrMod)

  # create variable importance table
  dt.rndFrImportance <- data.table(rndFrMod[["importance"]],
                                   keep.rownames = TRUE)

  # write the table
  WrtTable(dt.rndFrImportance, outputFolderPath, "RandomForestVariableImportance")


  PdfFilePath <- paste0(outputFolderPath, "/",
                        "randomForestVariableImportance",
                        startTimeStamp, ".pdf")

  pdf(PdfFilePath, width = 10, height = 7.5, paper = "USr", onefile = TRUE)
  varImpPlot(rndFrMod)
  dev.off()


  # write tables (random forest) ----------------------------------------

  WrtTable(dt.rfNodeSizeTest, outputFolderPath, "RandomForestNodeSizeCalibration")
  WrtTable(dt.rndFrstTest, outputFolderPath, "ScoresTestRandomForest")
  WrtTable(dt.rndFrstTrain, outputFolderPath, "ScoresTrainRandomForest")
  WrtTable(dt.rndFrstValid, outputFolderPath, "ScoresValidationRandomForest")

  # plot sample size evaluation --------------------------

  if (sampSizePerfEval == TRUE) {

    # create data table for plotting
    dt.sampSizePerfEval <- data.table(trainSize = vect.trainSize,
                                      acc.mlr = vect.acc.mlr,
                                      acc.rndFrst = vect.acc.rndFrst)

    WrtTable(dt.sampSizePerfEval,outputFolderPath,"SampleSizeEvalTable")

    # create empty plot lists
    plotSampSizePerfEval <- rep( list(list()), 1 )

    # create plot object
    plotSampSizePerfEval[[1]][[1]] <-
      ggplot() +
      geom_line(data = dt.sampSizePerfEval,
                 aes(x = trainSize, y = acc.mlr), linetype = 1) +
      geom_line(data = dt.sampSizePerfEval,
                aes(x = trainSize, y = acc.rndFrst), linetype = 2) +
      ggtitle("Performance vs. Sample Size") +
      ylab("model accuracy") +
      xlab("Size of Training Set") +
      ggThemeLcmsMetab()

    # create pdf
    PdfPlotsSingleFile(plotSampSizePerfEval, "landscape",
                       outputFolderPath, "performanceVsSampleSize")

  }


  # write tables model curation and validation --------------------------------------------

  # export table using defined function
  WrtTable(dt.modelAddVar, outputFolderPath, "ModelCurationTable")

  # write table
  WrtTable(dt.validation, outputFolderPath, "ModelValidationTable")

  # plot: added effect of new variables --------------------------------------

  # create empty plot lists
  plotVarModEr <- rep(list(list()), 5)

  # create plot object
  plotVarModEr[[1]][[1]] <-
    ggplot() +
    geom_point(data = dt.modelAddVar,
               aes(x = order, y = round(kFoldAvgAbsError,6), shape = include)) +
    ggtitle("Multi Variable Model K-Fold Mean Absolute Error") +
    ylab("average absolute error") +
    xlab("additional prediction variable") +
    scale_x_discrete(labels = dt.modelAddVar$var) +
    ggThemeLcmsMetab() +
    theme(axis.text.x = element_text(angle = 90))

  plotVarModEr[[2]][[1]] <-
    ggplot() +
    geom_point(data = dt.modelAddVar,
               aes(x = order, y = round(missIDrate,6), shape = include)) +
    ggtitle("Misidentification Rate") +
    ylab("misidentification rate") +
    xlab("additional prediction variable") +
    scale_x_discrete(labels = dt.modelAddVar$var) +
    ggThemeLcmsMetab() +
    theme(axis.text.x = element_text(angle = 90))

  # create plot object
  plotVarModEr[[3]][[1]] <-
    ggplot() +
    geom_point(data = dt.modelAddVar,
               aes(x = order, y = round(ErPercImp,1), shape = include)) +
    ggtitle("Multi Variable Error Reduction") +
    ylab("percent error reduction") +
    xlab("additional prediction variable") +
    scale_x_discrete(labels = dt.modelAddVar$var) +
    ggThemeLcmsMetab() +
    theme(axis.text.x = element_text(angle = 90))

  # create plot object
  plotVarModEr[[4]][[1]] <-
    ggplot() +
    geom_point(data = dt.modelAddVar,
               aes(x = order, y = round(pValAddVar,6), shape = include)) +
    ggtitle("Multi Variable Model p-value of Model Improvement") +
    ylab("p value") +
    xlab("additional prediction variable") +
    scale_x_discrete(labels = dt.modelAddVar$var) +
    ggThemeLcmsMetab() +
    theme(axis.text.x = element_text(angle = 90))

  # create plot object
  plotVarModEr[[5]][[1]] <-
    ggplot() +
    geom_point(data = dt.modelAddVar,
               aes(x = order, y = round(pValMaxWald,6), shape = include)) +
    ggtitle("Multi Variable Model Max Variable Wald Test p Value") +
    ylab("Max Variable Wald Test p Value") +
    xlab("additional prediction variable") +
    scale_x_discrete(labels = dt.modelAddVar$var) +
    ggThemeLcmsMetab() +
    theme(axis.text.x = element_text(angle = 90))

  # create pdf
  PdfPlotsSingleFile(plotVarModEr, "landscape",
                     outputFolderPath, "MultiVarAvgAbsErImprovement")



  # plot: Validation Mult Logi Regr ---------------------------------------------

  # convert label vector to character
  vect.varValidate <- as.character(vect.varValidate)

  # convert order column to factor
  dt.validation[, order := as.factor(order)]

  # create empty plot lists
  plotValidation <- rep( list(list()), 3 )

  # create plot object
  plotValidation[[1]][[1]] <-
    ggplot() +
    geom_point(data = dt.validation,
               aes(x = order, y = round(avgAbsEr,6), shape = include)) +
    ggtitle("Multi Variable Model Error") +
    ylab("average absolute error") +
    xlab("additional prediction variable") +
    scale_x_discrete(labels = vect.varValidate) +
    ggThemeLcmsMetab() +
    theme(axis.text.x = element_text(angle = 90))

  plotValidation[[2]][[1]] <-
    ggplot() +
    geom_point(data = dt.validation,
               aes(x = order, y = round(missIDrate,6), shape = include)) +
    ggtitle("Misidentification Rate") +
    ylab("misidentification rate") +
    xlab("additional prediction variable") +
    scale_x_discrete(labels = vect.varValidate) +
    ggThemeLcmsMetab() +
    theme(axis.text.x = element_text(angle = 90))

  # create plot object
  plotValidation[[3]][[1]] <-
    ggplot() +
    geom_point(data = dt.validation,
               aes(x = order, y = round(percErImprv,1), shape = include)) +
    geom_hline(yintercept = 0) +
    ggtitle("Multi Variable Error Reduction") +
    ylab("percent error reduction") +
    xlab("additional prediction variable") +
    scale_x_discrete(labels = vect.varValidate) +
    ggThemeLcmsMetab() +
    theme(axis.text.x = element_text(angle = 90))

  # create pdf
  PdfPlotsSingleFile(plotValidation, "landscape",
                     outputFolderPath, "MultLogiRegrValidationPlots")

  # plot Calibration Random Forest -------------------------------------

  # create empty plot lists
  plotRFCalibration <- rep( list(list()), 1 )

  # create plot object
  plotRFCalibration[[1]][[1]] <-
    ggplot() +
    geom_point(data = dt.rfNodeSizeTest,
               aes(x = nodeSize, y = round(avgAbsErCalibrate,6), shape = "Calibrate")) +
    geom_point(data = dt.rfNodeSizeTest,
               aes(x = nodeSize, y = round(avgAbsErTrain,6), shape = "Train")) +
    ggtitle("Random Forest Node Size Calibration") +
    ylab("average absolute error") +
    xlab("node size") +
    ggThemeLcmsMetab()

  # create pdf
  PdfPlotsSingleFile(plotRFCalibration, "landscape",
                     outputFolderPath, "RF_NodeSizeCalibrationPlot")

  # plot histogram ---------------------------------------------------

  titleText <- paste0("Results of Model Run on Test Data \n",
                      "Misidentification Rate: ", round(missIDrate * 100, 2),
                      "\u00B1", round(missIDrateMoe * 100, 2) , "% \n",
                      "Average Absolute Probability Error: ", round(avgAbsError, 3),
                      "\u00B1", round(avgAbsErrorMoe, 3))

  plotModelHist <- list(list())

  plotModelHist[[1]][[1]] <-
    ggplot(dt.MLogRegTest, aes(predProb, fill = class)) +
    geom_histogram(binwidth = 0.025, center = 0.0125) +
    scale_fill_manual(values = c("#999999", "#000000")) +
    ggtitle(titleText) +
    ylab('training peak count') +
    xlab('Probability Good') +
    ggThemeLcmsMetab()

  titleText <- paste0("Random Forest")

  plotModelHist[[2]] <- list()

  plotModelHist[[2]][[1]] <-
    ggplot(dt.rndFrstTest, aes(Good, fill = obsClass)) +
    geom_histogram(binwidth = 0.025, center = 0.0125) +
    scale_fill_manual(values = c("#999999", "#000000")) +
    ggtitle(titleText) +
    ylab('training peak count') +
    xlab('Probability Good') +
    ggThemeLcmsMetab()

  # create pdf
  PdfPlotsSingleFile(plotModelHist, "landscape",
                     outputFolderPath, "ModelPerformanceHist")

  # plot distribution density -------------------------------------

  plotDensity <- rep( list(list()), 2)

  vect.mLogReg.predProb.Bad <- subset(dt.MLogRegTest, class == FALSE)$predProb
  vect.mLogReg.predProb.Good <- subset(dt.MLogRegTest, class == TRUE)$predProb

  titleText <- "Multiple Logistic Regression: Good vs. Bad Distributions"

  plotDensity[[1]][[1]] <-
    ggplot() +
    geom_density(aes(vect.mLogReg.predProb.Bad, stat(count))) +
    geom_density(aes(vect.mLogReg.predProb.Good, stat(count))) +
    ggtitle(titleText) +
    xlab('Probability Good') +
    ggThemeLcmsMetab()

  # Next for random forest model

  vect.rndFrst.predProb.Bad <- subset(dt.rndFrstTest, obsClass == FALSE)$Good
  vect.rndFrst.predProb.Good <- subset(dt.rndFrstTest, obsClass == TRUE)$Good

  titleText <- "Random Forest: Good vs. Bad Distributions"

  plotDensity[[2]][[1]] <-
    ggplot() +
    geom_density(aes(vect.rndFrst.predProb.Bad, stat(count))) +
    geom_density(aes(vect.rndFrst.predProb.Good, stat(count))) +
    ggtitle(titleText) +
    xlab('Probability Good') +
    ggThemeLcmsMetab()

  # create pdf
  PdfPlotsSingleFile(plotDensity, "landscape",
                     outputFolderPath, "DensityPlots")

  # plot ROC  -------------------------------------------------------

  roc.multLog <- roc(response = dt.MLogRegTest$numClass,
                     predictor = dt.MLogRegTest$predProb,
                     quiet = TRUE)

  auc.multLog <- round(auc(roc.multLog),4)

  optThresh.multLog <- round(coords(roc.multLog, "best",
                              ret = "threshold",
                              best.method = "closest.topleft"),3)

  # create empty plot lists
  plotROC <- rep( list(list()), 3 )

  # plot title
  titleText <- paste0("ROC multiple logistic regression \n",
                      "Test AUC: ",auc.multLog,"\n",
                      "Optimum Threshold: ",optThresh.multLog)

  # create plot object
  plotROC[[1]][[1]] <-
    ggplot() +
    geom_step(aes(x = rev(roc.multLog$specificities),
                  y = rev(roc.multLog$sensitivities))) +
    geom_segment(aes(x = 0, y = 1, xend = 1,yend = 0),
                 alpha = 0.5) +
    scale_x_reverse(name = "Specificity (true neg / all neg)",
                    limits = c(1,0),
                    expand = c(0.001,0.001)) +
    scale_y_continuous(name = "Sensitivity (true pos/ all pos)",
                       limits = c(0,1),
                       expand = c(0.001, 0.001)) +
    ggtitle(titleText) +
    ggThemeLcmsMetab() +
    coord_equal()

  # random forest model

  roc.rndFrst <- roc(response = as.numeric(dt.rndFrstTest$obsClass),
                     predictor = dt.rndFrstTest$Good,
                     quiet = FALSE)

  auc.rndFrMod.test <- round(auc(roc.rndFrst),3)

  optThresh.RF <- round(coords(roc.rndFrst, "best",
                              ret = "threshold",
                              best.method = "closest.topleft"),3)

  titleText <- paste0("ROC random forest \n", "AUC: ",auc.rndFrMod.test,"\n",
                      "Optimum Threshold: ",optThresh.RF)

  plotROC[[2]][[1]] <-
    ggplot() +
    geom_step(aes(x = rev(roc.rndFrst$specificities),
                  y = rev(roc.rndFrst$sensitivities))) +
    geom_segment(aes(x = 0, y = 1, xend = 1,yend = 0),
                 alpha = 0.5) +
    scale_x_reverse(name = "Specificity (true neg / all neg)",
                    limits = c(1,0),
                    expand = c(0.001,0.001)) +
    scale_y_continuous(name = "Sensitivity (true pos/ all pos)",
                       limits = c(0,1),
                       expand = c(0.001, 0.001)) +
    ggtitle(titleText) +
    ggThemeLcmsMetab() +
    coord_equal()

  # Compare Mult-Logistic Regression and Random Forest

  titleText <- paste0("ROC - Both Models")

  plotROC[[3]][[1]] <-
    ggplot() +
    geom_step(aes(x = rev(roc.rndFrst$specificities),
                  y = rev(roc.rndFrst$sensitivities))) +
    geom_step(aes(x = rev(roc.multLog$specificities),
                  y = rev(roc.multLog$sensitivities))) +
    geom_segment(aes(x = 0, y = 1, xend = 1,yend = 0),
                 alpha = 0.5) +
    scale_x_reverse(name = "Specificity (true neg/ all neg)",
                    limits = c(1,0),
                    expand = c(0.001,0.001)) +
    scale_y_continuous(name = "Sensitivity (true pos/ all pos)",
                       limits = c(0,1),
                       expand = c(0.001, 0.001)) +
    ggtitle(titleText) +
    ggThemeLcmsMetab() +
    coord_equal()

  # create pdf
  PdfPlotsSingleFile(plotROC, "landscape",
                     outputFolderPath, "ROC")

  # plot proportion above threshold -----------------------------

  # function for computing the proportion remaining at each threshold
  PropRemainCol <- function(dt, scoreColName) {

    # dt is a datatable object
    # coreColName is the name of column that contains the scores

    # change score column name for easier syntax
    colNumScore <- match(scoreColName, colnames(dt))
    colnames(dt)[colNumScore] <- "score"

    # sort by score
    setkey(dt, score)

    # add rank column
    dt[, rank := seq_along(dt$score)]

    # add proportion remaining column
    dt[, propRemain :=
         (nrow(dt) - rank) / nrow(dt)]

    dt$rank <- NULL

    # reset score column name
    colnames(dt)[colNumScore] <- scoreColName

    return(dt)

  }

  # split tables by class
  dt.scoresTrueMultLogiReg <- subset(dt.MLogRegTest, class == TRUE)
  dt.scoresFalseMultLogiReg <- subset(dt.MLogRegTest, class == FALSE)

  dt.scoresTrueMultLogiReg <- PropRemainCol(dt.scoresTrueMultLogiReg, "predProb")
  dt.scoresFalseMultLogiReg <- PropRemainCol(dt.scoresFalseMultLogiReg, "predProb")

  plotThreshInclusion <- rep( list(list()), 2)

  titleText <- paste0("Multiple Logistic Regression: Good vs. Bad - Proportion Above Threshold",
                      "\n", "Optimum Threshold: ",optThresh.multLog)

  plotThreshInclusion[[1]][[1]] <-
    ggplot() +
    geom_step(data = dt.scoresTrueMultLogiReg,
              aes(x = predProb, y = propRemain)) +
    geom_step(data =  dt.scoresFalseMultLogiReg,
              aes(x = predProb, y = propRemain),
              linetype = 2) +
    geom_vline(xintercept = optThresh.multLog) +
    ggtitle(titleText) +
    ylab("Proportion of Cases Above Threshold") +
    xlab('Probability Good') +
    ggThemeLcmsMetab()

  # do the same for random forest

  # split tables by class
  dt.scoresTrueRF <- subset(dt.rndFrstTest, obsClass == TRUE)
  dt.scoresFalseRF <- subset(dt.rndFrstTest, obsClass == FALSE)

  dt.scoresTrueRF <- PropRemainCol(dt.scoresTrueRF, "Good")
  dt.scoresFalseRF <- PropRemainCol(dt.scoresFalseRF, "Good")

  titleText <- paste0("Random Forest: Good vs. Bad - Proportion Above Threshold",
                      "\n", "Optimum Threshold: ",optThresh.RF)

  plotThreshInclusion[[2]][[1]] <-
    ggplot() +
    geom_step(data = dt.scoresTrueRF,
              aes(x = Good, y = propRemain)) +
    geom_step(data =  dt.scoresFalseRF,
              aes(x = Good, y = propRemain),
              linetype = 2) +
    geom_vline(xintercept = optThresh.RF) +
    ggtitle(titleText) +
    ylab("Cumulative Proportion of cases") +
    xlab('Probability Good') +
    ggThemeLcmsMetab()

  PdfPlotsSingleFile(plotThreshInclusion, "landscape",
                     outputFolderPath, "ThresholdInclusion")

  # add model information to variable summary table ----------------

  # create model summary object
  modelSummary <- summary(logRegModel.final)

  # create model coefficients table
  dt.modCoeff <- data.table(modelSummary[["coefficients"]],
                            keep.rownames = TRUE)

  # add columns to dt.varSummary
  setkey(dt.modCoeff, rn)
  setkey(dt.varSummary, var)
  dt.varSummary <- merge(dt.varSummary, dt.modCoeff,
                         by.x = "var", by.y = "rn")


  # write variable summary table ------------------------------------

  # export table
  WrtTable(dt.varSummary, outputFolderPath, "VariableSummary")

  # copy model build data tables (is this a good approach?) -----------------------------------------

  if (copyInputDataFiles == TRUE) {


    # create output folder if it doesn't exist
    copyInputTablesPath <- paste0(outputFolderPath,"/",
                                  "InputTables",startTimeStamp)

    if (!dir.exists(copyInputTablesPath)) {
      print("create new output folder")
      print(copyInputTablesPath)
      dir.create(copyInputTablesPath)
    }

    # copy table listing data tables first
    copyToOutput(csvTableOfDataFiles, copyInputTablesPath)

    # read csvTableOfDataFiles
    dt.buildTables <- fread(csvTableOfDataFiles,
                            col.names = c("peakCharFile", "retainedPeakFile"))

    # create vector of names of other tables to copy
    vect.copyFiles <- c(dt.buildTables$peakCharFile,
                        dt.buildTables$retainedPeakFile)

    # use lapply to copy all files in list
    lapply(vect.copyFiles, copyToOutput, copyInputTablesPath)

    print("input data copied")
    logText <-
      UpdateLogText(logText,"data files used for model copied to output folder")
    UpdateLogFile(logFilePath, logText)

  }

  # export model -------------------------------------------------------

  # save model and identification objects

  # create copy of objects to be renamed with new name
  logisticRegModel <- logRegModel.final
  modelExternalSummaryText <- logText
  modelExternalFeatSummaryTable <- dt.featSummary.train
  modelExternalScriptVersion <- scriptVersion

  # list of objects to save
  saveList <- c("logisticRegModel",
                "modelExternalSummaryText",
                "modelExternalFeatSummaryTable",
                "modelExternalScriptVersion",
                "rndFrMod")
  # file path and name for file to be created
  saveFilePath <- paste0(outputFolderPath,"/",
                         "PeakPredictionModel",
                         startTimeStamp,".rda")
  # save to .rda file
  save(list = saveList, file = saveFilePath)

  logText <-
    UpdateLogText(logText,"model saved to output folder")
  UpdateLogFile(logFilePath, logText)


} # end buildModel if statement

# :: Apply Model -----------------------------------------------------

if (applyModel == TRUE) {

  # load model (if applicable) ---------------------------------------

  # load only if buildmodel is set to FALSE
  if (buildModel == FALSE) {

    modelFilePath <- paste0(parentDir,"/",loadModelFile)

    # load the model file indicated in user input section
    load(modelFilePath, verbose = TRUE)

    logText <-
      UpdateLogText(logText,"model loaded")
    UpdateLogFile(logFilePath, logText)

    # check that script version for model is the same
    # as the current version
    if (scriptVersion != modelExternalScriptVersion) {

      logText <-
        UpdateLogText(logText,"model version different from script version")
      UpdateLogFile(logFilePath, logText)

      stop("model was created from a different version of the script,
           please generate a new model from the training data using
           the newest version of the script")

    }

    # copy model file to the output folder
    copyToOutput(loadModelFile, outputFolderPath)

  }

  if (buildModel == TRUE) {

    # assign the built model to the logisticRegModel variable
    logisticRegModel <- logRegModel.final

  }

  # import prediction data -------------------------------------------------

  pkCharFilePredPath <- paste0(parentDir, "/", pkCharFilePred)

  print("begin import 6 attribute file")

  dt.AllCharPred <- fread(pkCharFilePredPath, sep = ",",
                          header = TRUE)

  print("complete import 6 attribute file")

  # Copy file to output folder
  file.copy(from = pkCharFilePredPath, to = outputFolderPath)

  # create feature summary table for prediction data -------------------------

  dt.featSummary.pred <- PeakPredictionTable(dt.AllCharPred)

  print("complete prediction table creation")

  logText <-
    UpdateLogText(logText,"prediction feature summary table created")
  UpdateLogFile(logFilePath, logText)

  # Apply the multi-logistic-regression model ---------------------------------------------------------

  # use predict to get a vector of the predicted probablilities of good
  vect.predProb <- predict(logisticRegModel,
                           newdata = dt.featSummary.pred,
                           type = "response")

  logText <-
    UpdateLogText(logText,"predicted probabilities determined")
  UpdateLogFile(logFilePath, logText)

  # rename columns for easy reference
  colnames(dt.featSummary.pred)[2] <- "mz"
  colnames(dt.featSummary.pred)[3] <- "RT"

  # create prediction table for multiple logistic regression model
  dt.featMultLogRegrPred <- subset(dt.featSummary.pred,
                                   select = c("featID",
                                              "mz", "RT"))

  # add prediction vector to table
  dt.featMultLogRegrPred[, predProb := vect.predProb]

  #  create printable predicted probability
  dt.featMultLogRegrPred[, predProbRnd :=
                        formatC(vect.predProb, format = "e", digits = 2)]

  # add to dt.featSummary.pred
  dt.featSummary.pred[, multLogiRegrScr := vect.predProb]

  # add rank column to summary table ---------------------------------------

  # order by predProb
  setorder(dt.featMultLogRegrPred, -predProb)

  # create rank column
  dt.featMultLogRegrPred[, probRank := 1:nrow(dt.featMultLogRegrPred)]

  # add data and file extension to the file name
  OutputFileName <- paste0("PredictionMultiLogisticRegr",
                           startTimeStamp)

  # write table
  WrtTable(dt.featMultLogRegrPred, outputFolderPath, OutputFileName)

  # plot histogram predicted values  ------------------------------------------

  # create function for vertical line positions
  # set_breaks = function(limits) {
  #   seq(limits[1], limits[2], by = .1)
  # }
  #
  vect.breaks <- seq(0, 1, by = .1)

  # create plot list
  plotPrediction <- list()
  plotPrediction[[1]] <- list()

  # create plot object
  plotPrediction[[1]][[1]] <-
    ggplot() +
    scale_x_discrete(limits = vect.breaks, labels = vect.breaks, expand = c(0,0)) +
    geom_histogram(data = dt.featMultLogRegrPred,
                aes(predProb), binwidth = 0.025, center = 0.0125) +
    ggtitle("Multiple Logistic Regression Score Distribution") +
    ylab("Peak Count by Bin") +
    xlab('Probability Peak is Acceptable') +
    ggThemeLcmsMetab()

  # generate pdf
  PdfPlotsSingleFile(plotPrediction, "landscape", outputFolderPath,
                     "PredictionHistogramMultLogiRegr")

  logText <-
    UpdateLogText(logText,"predicted probability values histogram created")
  UpdateLogFile(logFilePath, logText)

  # export MZmine databases multiple logistic regression -------------------------------------------------

  # create peak number cutoff table
  dt.MZmineDbPkNumCut <- subset(dt.featMultLogRegrPred,
                   probRank <= peakNumCutoff,
                   select = c("mz","RT", "predProbRnd", "probRank"))


  # add data and file extension to the file name
  OutputFileName <- paste0(outputFolderPath,"/",
                           "MZmineDB_",peakNumCutoff,"peaks_MultiLogRegr",
                           startTimeStamp,".csv")

  #Write the adjusted dataframe to csv
  fwrite(dt.MZmineDbPkNumCut, file = OutputFileName, row.names = TRUE)

  logText <-
    UpdateLogText(logText,"Number of peaks MZmine database file written")
  UpdateLogFile(logFilePath, logText)

  # create probability cutoff table
  dt.MZmineDbProbCut <- subset(dt.featMultLogRegrPred,
                               predProb >= probCutoff,
                               select = c("mz","RT", "predProbRnd", "probRank"))

  # add data and file extension to the file name
  OutputFileName <- paste0(outputFolderPath,"/",
                           "MZmineDB_",probCutoff,
                           "probabilityCutoff_MultiLogRegr",
                           startTimeStamp,".csv")

  #Write the adjusted dataframe to csv
  fwrite(dt.MZmineDbProbCut, file = OutputFileName, row.names = TRUE)

  # apply Random Forest Model ---------------------------------------

  dt.rndFrstPred <- as.data.table(predict(rndFrMod,
                                          dt.featSummary.pred,
                                          type = "prob"))
  colnames(dt.rndFrstPred) <- c("Bad","Good")

  # add columns from prediction data table
  dt.rndFrstPred$featID <- dt.featSummary.pred$featID
  dt.rndFrstPred$mz <- dt.featSummary.pred$"mz"
  dt.rndFrstPred$RT <- dt.featSummary.pred$"RT"

  dt.featSummary.pred$RFScr <- dt.rndFrstPred$Good

  # export feature summary table ----------------------------------------

  # add data and file extension to the file name
  OutputFileName <- paste0("PredictionFeatureSummaryTable",
                           startTimeStamp)

  # Write the summary table to csv
  WrtTable(dt.featSummary.pred, outputFolderPath, OutputFileName)

  logText <-
    UpdateLogText(logText,"predicted feature summary table written")
  UpdateLogFile(logFilePath, logText)

  # Export random forest prediction MZmine databases ---------------------------------------

  dt.MZmineDbProbCutRndFrst <- subset(dt.rndFrstPred,
                                      Good >= probCutoff,
                                      select = c("mz","RT", "Good"))

  # add data and file extension to the file name
  OutputFileName <- paste0(outputFolderPath,"/",
                           "MZmineDB_",probCutoff,
                           "probabilityCutoffRndFrst",
                           startTimeStamp,".csv")

  #Write the adjusted dataframe to csv
  fwrite(dt.MZmineDbProbCutRndFrst, file = OutputFileName, row.names = TRUE)

  # plot histogram predicted values  ------------------------------------------

  # create function for vertical line positions
  # set_breaks = function(limits) {
  #   seq(limits[1], limits[2], by = .1)
  # }
  #
  vect.breaks <- seq(0, 1, by = .1)

  # create plot list
  plotPrediction <- list()
  plotPrediction[[1]] <- list()

  # create plot object
  plotPrediction[[1]][[1]] <-
    ggplot() +
    scale_x_discrete(limits = vect.breaks, labels = vect.breaks, expand = c(0,0)) +
    geom_histogram(data = dt.rndFrstPred,
                   aes(Good), binwidth = 0.025, center = 0.0125) +
    ggtitle("Predicted Probability of Acceptable Distribution") +
    ylab("Peak Count by Bin") +
    xlab('Probability Peak is Acceptable') +
    ggThemeLcmsMetab()

  # generate pdf
  PdfPlotsSingleFile(plotPrediction, "landscape", outputFolderPath,
                     "PredictionHistogramRndForest")

} # end apply model if statement

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
  UpdateLogText(logText,"end script")
UpdateLogFile(logFilePath, logText)

if (1 == 1) {
print("script complete")
print(runTime(startTime))
}
