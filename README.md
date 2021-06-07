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
