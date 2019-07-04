# Description -----------------------------------------------------------------

# use peak group matrices to create and apply neural network classifier

# User Input ------------------------------------------------------------------

# path to directory where training and test data is stored
inputDataParentDir <- "E:\\Example\\Example"

# Directory where output folder will be created
outPutParentDir <- "E:\\Example2\\Example2"

# End User input ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

# Advanced user input (typically don't need to change) ------------------------

# Training data .rds file
RdsTraining <- "TrainingData.rds"

# Training data label table csv
CsvTrainingLabels <- "TrainingLabels.csv"

# Column of training labels
colNumTrainLabels <- 3

# Validation/calibration data .rds file
RdsValid <- "ValidationData.rds"

# Validation/calibration data label table csv
CsvValidLabels <- "ValidationLabels.csv"

# Column of Validation/calibration labels
colNumValidLabels <- 3

# Test data .rds file
RdsTest <- "TestData.rds"

# Validation/calibration data label table csv
CsvTestLabels <- "TestLabels.csv"

# Column of Validation/calibration labels
colNumTestLabels <- 3

# Library (packages)--------------------------------------------------------
library(data.table) # faster version of data frames

library(ggplot2) # for plotting
library(gridExtra) # for organizing plots
library(pROC) # for receiver operator curves

library(keras) # for neural networks

library(lcmsMetab) # LCMS data processing tools

# Timing --------------------------------------------------------------
startTime <- Sys.time()
startTimeStamp <- format(startTime, "(%Y%m%dh%H%M)")

# set or create the output folder -------------------------------------------

OutputFolderName <- "peakNeuralNetModel"

# create output folder if it doesn't exist
outputFolderPath <- paste0(outPutParentDir,"/",
                           OutputFolderName, startTimeStamp)

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

# import training data
importFilePath <- paste0(inputDataParentDir,"/", RdsTraining)
images_train <- readRDS(file = importFilePath)

# import validation data
importFilePath <- paste0(inputDataParentDir,"/", RdsValid)
images_valid <- readRDS(file = importFilePath)

# import validation data
importFilePath <- paste0(inputDataParentDir,"/", RdsTest)
images_test <- readRDS(file = importFilePath)

# load labels as vectors
importFilePath <- paste0(inputDataParentDir,"/", CsvTrainingLabels)
LabelsTrain <- fread(file = importFilePath, select = c(colNumTrainLabels))
LabelsTrain <- LabelsTrain[[1]]

importFilePath <- paste0(inputDataParentDir,"/", CsvValidLabels)
LabelsValid <- fread(file = importFilePath, select = c(colNumValidLabels))
LabelsValid <- LabelsValid[[1]]

importFilePath <- paste0(inputDataParentDir,"/", CsvTestLabels)
LabelsTest <- fread(file = importFilePath, select = c(colNumTestLabels))
LabelsTest <- LabelsTest[[1]]

# import full test data table
importFilePath <- paste0(inputDataParentDir,"/", CsvTestLabels)
dt.testSet <- fread(file = importFilePath)


# Explore (for dev - comment out) ----------------------------------------------

# dim(images_train)
#
# image(images_train[200,,], useRaster = TRUE, axes = FALSE,
#       col = grey(seq(0, 1, length = 256)))


# reshape the data -----------------------------------------------------------

imageHeight <- dim(images_train)[2]
imageWidth <- dim(images_train)[3]

images_train <- array_reshape(images_train, c(nrow(images_train),
                                              imageHeight, imageWidth, 1))
images_valid <- array_reshape(images_valid, c(nrow(images_valid),
                                              imageHeight, imageWidth, 1))
images_test <- array_reshape(images_test, c(nrow(images_test),
                                            imageHeight, imageWidth, 1))

# Build the model (move to separate file) ---------------------------------

# _setup the layers -------------------------------------------------

model <- keras_model_sequential()
model %>%
  # first layer flattens the matrix to a single vector of 4960 pixels (28 x 28)
  layer_conv_2d(filters = 32, kernel_size = c(3,3), activation = 'relu',
                input_shape = c(31, imageWidth, 1)) %>%
  layer_max_pooling_2d(pool_size = c(2,2), strides = c(2,2)) %>%
  layer_dropout(rate = 0.2) %>%
  layer_conv_2d(filters = 16, kernel_size = c(3,3), activation = 'relu') %>%
  layer_max_pooling_2d(pool_size = c(2,2), strides = c(2,2)) %>%
  layer_dropout(rate = 0.2) %>%
  layer_flatten() %>%
  # dense layers are fully connected all 128 neurons connected to each other
  layer_dense(units = 64, activation = 'relu') %>%
  layer_dropout(rate = 0.5) %>%
  # The last layer is densely connected with a single output node. Using the
  # sigmoid activation function, this value is a float between 0 and 1,
  # representing a probability, or confidence level.
  # layer_dense(units = 1, activation = "sigmoid")
  # Project onto a single unit output layer, and squash it with a sigmoid
  layer_dense(1) %>%
  layer_activation("sigmoid")

summary(model)


# _compile the model ----------------------------------------------------------

# binary crossentropy used for loss function

model %>% compile(
  optimizer = 'adam',
  loss = 'binary_crossentropy',
  metrics = list('accuracy')
)

# _train the model -------------------------------------------------------

# Display training progress by printing a single dot for each completed epoch.
print_dot_callback <- callback_lambda(
  on_epoch_end = function(epoch, logs) {
    if (epoch %% 80 == 0) cat("\n")
    cat(".")
  }
)

# The patience parameter is the amount of epochs to check for improvement.
early_stop <- callback_early_stopping(monitor = "val_loss", patience = 5)


history <- model %>% fit(
  x = images_train,
  y = LabelsTrain,
  epochs = 40,
  batch_size = 512,
  validation_data = list(images_valid, LabelsValid),
  verbose = 1,
  callbacks = list(early_stop, print_dot_callback)
)

# Visualize training --------------------------------------------------------

# plot(history, metrics = "acc", smooth = FALSE)
plot(history)

# plot training ------------------------------------------------------------

# create data table from training tracking
dt.trainingStats <- data.table(epoch = 1:length(history[["metrics"]][["loss"]]),
                               loss = history[["metrics"]][["loss"]],
                               val_loss = history[["metrics"]][["val_loss"]],
                               acc = history[["metrics"]][["acc"]],
                               val_acc = history[["metrics"]][["val_acc"]])

WrtTable(dt.trainingStats, outputFolderPath, "modelTrainingStats")

plot.trainingStats <-
  plot(history) +
  ggThemeLcmsMetab()

ggsave("plotTrainingStats.pdf", plot = plot.trainingStats, device = "pdf",
       path = outputFolderPath,
       width = 8.5,
       height = 11,
       units = "in")


# Save model -----------------------------------------------------------------

exportFilePath <- paste0(outputFolderPath,"/","model", startTimeStamp,".h5")
save_model_hdf5(model, filepath = exportFilePath)


# apply the model to the test set --------------------------------------------

# create data table of predicted probabilities
dt.testPrediction  <- predict_proba(model, images_test)
dt.testPrediction <- data.table(dt.testPrediction)
colnames(dt.testPrediction)[1] <- "predProb"

# add pred probability to test set data table
dt.testSet[, predProb := dt.testPrediction$predProb]

# ROC ------------------------------------------------------------------

roc <- roc(response = dt.testSet$label,
                   predictor = dt.testSet$predProb,
                   quiet = TRUE)

auc <- round(auc(roc),4)

# optimum threshold --------------------------------------------------------

optThresh <- round(coords(roc, "best",
                                  ret = "threshold",
                                  best.method = "closest.topleft"),3)

# optimum threshold table
dt.optThresh <- copy(dt.testSet)

dt.optThresh[, tp := ifelse(label == 1 & predProb >= optThresh, 1, 0)]
dt.optThresh[, fp := ifelse(label == 1 & predProb < optThresh, 1, 0)]
dt.optThresh[, tn := ifelse(label == 0 & predProb < optThresh, 1, 0)]
dt.optThresh[, fn := ifelse(label == 0 & predProb >= optThresh, 1, 0)]

WrtTable(dt.optThresh, outputFolderPath, "optimumThresholdTable")

# confusion matrix values
optThreshTruePos <- sum(dt.optThresh$tp)
optThreshTrueNeg <- sum(dt.optThresh$tn)
optThreshFalsePos <- sum(dt.optThresh$fp)
optThreshFalseNeg <- sum(dt.optThresh$fn)

# calculate accuracy at optimum threshold
optThreshAccuracy <- (optThreshTruePos + optThreshTrueNeg) /
  (optThreshTruePos + optThreshTrueNeg + optThreshFalsePos + optThreshFalseNeg)

# compute Sensitivity at optimum threshold
optThreshSensitivity <- coords(roc, optThresh,
                               input = "threshold",
                               ret = "sensitivity")

# compute 1-specificity (False positive rate) at optimimum threshold
optThreshFPR <- coords(roc, optThresh,
                       input = "threshold",
                       ret = "1-specificity")


# plot ROC ---------------------------------------------------------------

# plot title
titleText <- paste0("ROC Deep Neural Network \n",
                    "Test AUC: ",auc,"\n",
                    "Optimum Threshold: ", optThresh, "- Acc at Optimum: ",
                    round(optThreshAccuracy,3)*100,"%")

# create plot object
plotROC <-
  ggplot() +
  geom_step(aes(x = rev(roc$specificities),
                y = rev(roc$sensitivities))) +
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

ggsave("plotROC.pdf", plot = plotROC, device = "pdf",
       path = outputFolderPath,
       width = 11,
       height = 8.5,
       units = "in")


# plot histogram --------------------------------------------------------------

# plot title
titleText <- paste0("Performance Histogram Deep Neural Network \n",
                    "Test AUC: ",auc,"\n",
                    "Optimum Threshold: ",optThresh)

# plot object
plotModelHist <-
  ggplot(dt.testSet, aes(predProb, fill = as.factor(label))) +
  geom_histogram(binwidth = 0.025, center = 0.0125) +
  scale_fill_manual(values = c("#999999", "#000000")) +
  ggtitle(titleText) +
  ylab('training peak count') +
  xlab('Probability Good') +
  ggThemeLcmsMetab()

# save plot
ggsave("plotHist.pdf", plot = plotModelHist, device = "pdf",
       path = outputFolderPath,
       width = 11,
       height = 8.5,
       units = "in")

# plot cumulative threshold ------------------------------------------------

# split tables by class
dt.scoresTrue <- subset(dt.testSet,
                        label == 1)
dt.scoresFalse <- subset(dt.testSet,
                         label == 0)

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

dt.scoresTrue <- PropRemainCol(dt.scoresTrue, "predProb")
dt.scoresFalse <- PropRemainCol(dt.scoresFalse, "predProb")

titleText <- paste0("Proportion Above Score Threshold", "\n",
                    "Optimum Threshold: ", optThresh, " (tpr = ",
                    round(optThreshSensitivity,3), ", fpr = ",
                    round(optThreshFPR,3),")")

# create plot object
plotThreshInclusion <-
  ggplot() +
  geom_step(data = dt.scoresTrue,
            aes(x = predProb, y = propRemain)) +
  geom_step(aes(x = c(dt.scoresFalse$predProb,1),
                c(y = dt.scoresFalse$propRemain,0)),
            size = 2) +
  geom_vline(xintercept = optThresh, linetype = 2) +
  # geom_hline(yintercept = optThreshSensitivity) +
  # geom_hline(yintercept = optThreshFPR) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  ggtitle(titleText) +
  ylab("Proportion Above Threshold") +
  xlab("Score Threshold") +
  ggThemeLcmsMetab() +
  theme(axis.text.x = element_text(angle = 90))

# export
ggsave("plotTprFpr.pdf", plot = plotThreshInclusion, device = "pdf",
       path = outputFolderPath,
       width = 11,
       height = 8.5,
       units = "in")

# check on mislabels ----------------------------------------------------------

# add average absolute error
dt.testSet[, avgAbsEr := round(abs(label - predProb),4)]

# write table
WrtTable(dt.testSet, outputFolderPath, "testSetPredictionTable")

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