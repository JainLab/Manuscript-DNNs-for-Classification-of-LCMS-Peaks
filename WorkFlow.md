# Typical Workflow for Utilizing the Scripts

## Using R Scripts

- Download and install the most recent version of R for your operating system (https://www.r-project.org/)

- We recommend using the popular free IDE RStudio (https://www.rstudio.com/) for running the scripts, though other IDEs can be used or the interpreter may be run from the command line.
- The scripts are heavily commented to help interested users understand the algorithm.
- Each script has 3 sections that a user should review:
  - Each script begins with a *Description* indicating briefly the purpose and function of the script
  - Each script includes a *User Input* section where the variables that need to be adjusted by the user are defined. Each variable has a brief description to help users assign the appropriate value
  - Each script includes a *Library* section indicating the R packages that need to downloaded and installed for the script to run. 

## Common Workflow

- Process your data using your preferred MZmine 2  (http://mzmine.github.io/) workflow so that you get an aligned peak list. The settings included in the manuscript supplement could be used as a starting point, but they would need to be adjusted to perform well with your method.

- Note that one benefit of using the machine learning scripts is that you can use less restrictive settings in MZmine 2's *Chromatogram deconvolution* algorithms, because a well trained machine learning model is very effective in removing false peaks

- Both scripts perform better when the data is well aligned in the retention time dimension. Retention time drift can be corrected and new adjusted mzXMLs can be created using the method described in:

  Watrous, J. D.; Henglin, M.; Claggett, B.; Lehmann, K. A.; Larson, M. G.; Cheng, S.; Jain, M.
  Visualization, Quantification, and Alignment of Spectral Drift in Population Scale Untargeted Metabolomics Data. *Analytical Chemistry* **2017**, *89* (3), 1399–1404. https://doi.org/10.1021/acs.analchem.6b04337

## Application of the Simple Machine Learning Script

- Export the MZmine 2 aligned peak list to a CSV file using the following settings:
  for the common elements select:
    row m/z
    row retention time
  
  for the data file elements select:
    peak duration
    peak height
    peak area
    peak FWHM
    peak tailing factor
    peak asymmetry factor

    make sure no other check boxes are checked as the column order must match the above order for all files

- Make sure you have the following R packages installed:

  data.table, ggplot2, gridExtra, randomForest, pROC, lcmsMetab

- Open SimpleMachineLearningPeakPedictor.R in R

- Scroll Down to [# User Input #] and follow the instructions for each input variable. The same script is used for building and applying models. Use the inputs to indicate the application



## Application of the Deep Neural Network Scripts

- Generate peak windows from MZmine 2 output. Different methods can be used for developing the window boundaries, using the apex values, or retention time start and end and m/z min and m/z max. The further the separation between isobaric peaks, the easier it will be to generate appropriate windows that do not split peaks or contain multiple peaks.
- Use mzXMLtoTables.R script to create MS1 tables from the mzXML files.
- Use GenStudyWndowApexIntensityTable.R to extract the max intensities from the peak windows
- Use peakGroupImageMatrixAndPlot.R to create peak image plots and matrices
- Manually review peak window plots to build training set for model generation.
  - First use createPeakReviewTrackingTable.R to create a review table
  - Then use peakReviewApp.R to review the plot images. This script opens a GUI for viewing and classifying peak groups using the plot images created by peakGroupImageMatrixAndPlot.R. When reviewing peaks, each image is classified as "Good" or acceptable by default. If the plot image shows a peak which is not acceptable it should be indicated by pressing the space bar or the [*Bad (space)*] button before advancing to the next image. Peaks that are not clearly good or bad but near the borderline should not be used to train the classifier. Flag these images with the [*(un)Flag Borderline (b)*] button.  There is also a button for flagging images where the window edges are not properly capturing the peaks, either capturing multiple peaks or splitting a single peak.
  - The more peaks that are reviewed, the more examples will be available to train the classifier and the better the classifier will perform. One option is to review peaks, train the classifier, review the performance from the plots, and if the performance is not adequate, review additional peak plots. 
- Prepare the data for training and testing the deep neural network model using NeuralNetworkDataPrep.R
- Use buildNeuralNetworkModel.R to train and test the neural network. This script will export the deep neural network model to an *.h5 file so it can be applied to other datasets. This script also exports a number of plots that will help you evaluate the performance of the model and determine what threshold you want to use to determine if a peak will remain in your dataset
- Apply the model using ApplyPeakPredictionModel.R. This script will provide you with a score determined by the model.
- Use filterOutBadPeaks.R to remove windows with scores below the selected threshold. You can also filter using the manual classification completed with peakReviewApp.R
- Use SubsetTableToVector.R to remove the windows classified as unacceptable from the peak intensity table generated earlier using GenStudyWindowApexIntensityTable.R. If windows were adjusted after that table was generated (perhaps to correct a window with multiple peaks), use GenStudyWindowApexIntensityTable.R again to create a new table of peak intensities for each mzXML file.

- From this point you may post-process or analyze the data as you might do if you had exported a peak height table from MZmine 2.