# TCSS499 - Application of Natural Language Processing and Machine Learning to Radiology Reports

Most part of the programs are written as functional programming.



In Analysis.R:

* Lines from 6 ~ 33: load custom functions, key packages, and parameters.
* Lines from 36 ~ 505: helper functions to run the program
* Lines from 506 ~ 521: Run the program
* Lnes from 523 ~ 633: Analyze results
* Lines from 636 ~ 668: Generate graphs



Here are the list of functions
* Load data:
  * loadData()
  * loadDataSubset()
* Create document-term matrix:
  * createCorpus()
  * createWideDTM()
* Helper functions for machine learning algorithms:
  * featureSelection()
  * brierScore()
  * classification()
* Functions for machine learning algorithms:
  * knn_fxn()
  * rf_fxn()
  * gbm_fxn()
  * xgboost_fxn()
  * adaboost_fxn()
* Run the each of the machine learning algorithms above and collect result data:
  * collectData()
* Analyze Result:
  * contingency_table_analysis()
  * brier_score_analysis()
  * auc_analysis()


Key Features: 
* 1-2 gram model
* Threshold value for feature selection: 0.5

_______________________________________________________________________________________________________________________________________
Results in Google Folder

Images:
* Wordcloud_Pneumonia.tiff
  * Shows the most frequent terms in negative and postivie sides.
* Heatmap_Pneumonia.tiff
  * Heatmap that depicts how many times each of the top grams appears in each report.
* Barchart_Pneumonia.tiff
  * Barchar that shows the percent of reports that contain each gram.
* Histogram_GramCounts.tiff
  * Histogram that depicts the distribution of gram usage across reports.
* Brier_Scores_Graph.png
  * Graph that shows the result of brier scores on different number of reports. The lower the brier scores, the better the model is at predicting.
    The graph is horizontal bar graph, where the x-axis shows the brier scores, and y-axis shows the 5 different models. 
* AUC_Graph.png
  * Graph that shows the result of AUC on different number of reports. The higher the AUC, the better the model is at predicting. 
    The graph is horizontal bar graph, where the x-axis shows the AUC, and y-axis shows the 5 different models. 

Excel:
* Results.xlsx
  * Results shows the result of each run. Here are the list of information that the file contains:
    * Number of reports
    * Randome Seed
    * Patient ID
    * Report ID
    * Path
    * Result of each models (1 indicates that the model classified the report correctly, 0 indicates that the model classified the report incorrectly)
    * TN/TF/FP/FN for each models
* Contingency_Table.xlsx
  * It contains the average TN/TF/FP/FN rates for each run. It also contains average sensitivity and specificity values.
* Brier_scores.xlsx
  * It contains the results of brier scores on different number of reports for each model.
* AUC.xlsx
  * It contains the results of AUC on different number of reports for each model.
