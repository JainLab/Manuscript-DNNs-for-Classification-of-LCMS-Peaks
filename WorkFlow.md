# Typical Workflow for Utilizing the Scripts

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

- Scroll Down to [# User Input #] and follow the instructions for each input variable. The same script is used for building and applying models.



## Application of the Deep Neural Network Scripts

- Generate peak windows from MZmine 2 output. 
- Use mzXMLtoTables.R script to create MS1 tables from the mzXML files
- Use GenStudyWndowApexIntensityTable.R to extract the max intensities from the peak windows
- use peakGroupImageMatrixAndPlot.R to create peak image plots and matrices
- Manually review peak window plots to build training set for model generation. First use createPeakReviewTrackingTable.R to create a review table, then use peakReviewApp.R to review the plot images

