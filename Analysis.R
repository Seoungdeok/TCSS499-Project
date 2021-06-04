# The code in this document assumes that the reports have been stored in a csv
# file, as created using the code in the file "StructureData.R". These reports
# came from the MIMIC-CXR database. The file "mimic-cxr-2.0.0-negbio.csv" came
# from the MIMIC-CXR-JPG database.

# Load our custom functions.
source("fxns.R")

# Load key packages.
library(magrittr) # Pipe operator.
library(ggplot2) # Plotting.
library(cowplot)
theme_set(theme_cowplot())
library(tm) # Text mining.
library(kNN)
library(randomForest)
library(pROC)
library(gbm)
library(dplyr)
library(xlsx)

# Parameters for the analysis.
stem_words <- TRUE # Should words be stemmed or not.
remove_stopwords <- TRUE # Should stop words be removed or not.
words.to.be.rm <- c(# Comma separated list of words or phrases that should be removed.
  stopwords("en"),
  "final report"
)
sort(words.to.be.rm)
words.to.be.rm <- removeWords(words.to.be.rm, c("no", "not", "nor", "neither"))
gram_size <- 1:2 # Gram size. You can set this to a single integer or a range, e.g. 1:4





######### Load Data #########
loadData <- function() {
  # Set parameters.
  data_folder <- "data" # Folder containing file_map.csv and mimic-cxr-2.0.0-negbio.csv
  file_map_path <- file.path(data_folder, "file_map.csv")
  report_labels_path <- file.path(data_folder, "mimic-cxr-2.0.0-negbio.csv.gz")
  
  ## Create a results directory.
  #if(dir.exists("results")){
  #  unlink("results", recursive = TRUE, force = TRUE)
  #}
  #dir.create("results")
  
  # Load the data.
  file_map <- read_FileMap(file_map_path)
  report_labels <- readr::read_csv(report_labels_path) %>%
    dplyr::mutate_at(., dplyr::vars(-subject_id, -study_id), function(x){factor(
      x,
      levels = c(-1,0,1),
      labels = c("n","u","p")
    )}) %>%# Now -1=n, 0=u, 1=p, for negative, unspecified, and positive claims.
    dplyr::rename(# Make them match file_map.
      patient_id = "subject_id",
      report_id = "study_id"
    ) %>%
    dplyr::mutate(
      patient_id = paste0("p", patient_id),
      report_id = paste0("s", report_id)
    )
  
  # Merge the data sets.
  merged <- dplyr::inner_join(
    file_map, report_labels,
    by = c("patient_id", "report_id")
  )
  return(merged)
}
data <- loadData()

# Set the random seed, then randomly select some reports.
loadDataSubset <- function(RANDOM_SEED, num_reports) {
  data_subset <- dplyr::filter(data, !is.na(Pneumonia) & (Pneumonia == "n" | Pneumonia == "p")) # Select only reports for which Pneumonia is not NA and Pneumonia is "p" or "n".
  data_subset$Pneumonia <- factor(data_subset$Pneumonia) # Drop "u" level from Pneumonia factor.
  set.seed(RANDOM_SEED)
  keepers <- sample(1L:nrow(data_subset), min(num_reports, nrow(data_subset))) # Select 1000 documents randomly
  data_subset <- data_subset[keepers,]
  return(data_subset)
}

######### Apply tm package #########
# Create and clean up the corpus.
createCorpus <- function(data_subset) {
  corpus <- VCorpus(VectorSource(data_subset$report)) %>%
    tm_map(tolower) %>% # Make all lower case.
    tm_map(function(x){gsub("[[:digit:]]","",x)}) %>% # Remove numbers
    tm_map(function(x){gsub("'"," ",x)}) %>% # Apostrophes to spaces.
    tm_map(function(x){gsub("_+[[:alpha:]]+","_",x)}) %>% # Remove PII indicators.
    tm_map(function(x){gsub("aren't","are not",x)}) %>% # "aren't" to "are not"
    tm_map(function(x){gsub("isn't","is not",x)}) %>% # "isn't" to "is not"
    tm_map(function(x){gsub("can't","can not",x)}) %>% # "can't" to "can not"
    tm_map(function(x){gsub("haven't","have not",x)}) %>% # "haven't" to "have not"
    tm_map(removePunctuation) %>% # Remove punctuation, including _
    (function(x){# Remove stop words if remove_stopwords is TRUE.
      if(remove_stopwords){
        return(tm_map(x, removeWords, words.to.be.rm))
        words.to.be.rm
      }else{
        return(x)
      }
    }) %>%
    tm_map(stripWhitespace) %>% # Remove whitespace in the middle.
    tm_map(trimws) %>% # Remove leading and trailing spaces.
    (function(x){# Stem words if stem_words is TRUE.
      if(stem_words){
        return(tm_map(x, stemDocument))
      }else{
        return(x)
      }
    })
  return(corpus)
}

# Create a DTM by n-gram count table.
createWideDTM <- function(corpus) {
  dtm_wide <- 1L:length(corpus) %>% # For each document in the corpus...
    lapply(
      .,
      function(j){
        x <- unlist(corpus[j]) # Retrieve the document and convert to a string.
        splits <- strsplit(x, " ") %>% # Split it on spaces and...
          unlist()# Convert it to a vector where each element is a single word.
        lapply(# For each gram size specified in gram_size...
          gram_size,
          function(i){
            splits %>%# Retrieve all possible n grams from the document...
              NLP::ngrams(., i) %>%
              lapply(., paste, collapse = " ") %>%
              unlist()
          }
        ) %>%
          unlist() %>% # Convert the resulting list of n-grams to a single vector...
          table() %>% # And count the number of times each gram appears.
          as.matrix() %>% # Convert the count table to a matrix...
          as.data.frame(stringsAsFactors = FALSE) %>% # Convert to a data.frame
          tibble::rownames_to_column("gram") %>% # Move the row names (the grams) to a new column named "gram".
          dplyr::mutate(DOCUMENT_ID = j) %>% # Add the document id to a new column named "DOCUMENT_ID".
          dplyr::rename(count = "V1") %>% # Rename the column "V1" (number of times each gram appeared) to "count"
          tibble::as_tibble() # Convert to a tibble because they are easier to work with than a data.frame.
      }
    ) %>%
    dplyr::bind_rows() %>%# The above tibbles of gram counts for each document are stacked on top of each other...
    tidyr::pivot_wider(names_from = gram, values_from = count, values_fill = 0L) %>% # And then changed from tall format to wide format, adding zeros to all missing values.
    dplyr::arrange(DOCUMENT_ID) # Sort by DOCUMENT_ID to match the order in merged_subset.  
}


# Apply feature selection to DTM
featureSelection <- function(data_subset, dtm_wide, threshold) {
  # Select words that appear at least 10 times.
  dtm <- dtm_wide[,colSums(as.matrix(dtm_wide)) > 10]
  
  # Merge with Pneumonia data.
  dtm <- dplyr::mutate(dtm, Pneumonia = data_subset$Pneumonia)
  
  # Split by Pneumonia == "p" and Pneumonia == "n".
  dtm_positive <- dplyr::filter(dtm, Pneumonia == "p")
  dtm_negative <- dplyr::filter(dtm, Pneumonia == "n")
  
  # Remove columns called "DOCUMENT_ID" and "Pnemonia", and sum up the columns.
  grammed_p_freq <- colSums(as.matrix(dtm_positive[,c(-1,-ncol(dtm_positive))]))
  grammed_n_freq <- colSums(as.matrix(dtm_negative[,c(-1,-ncol(dtm_negative))]))
  
  # Identify terms that are frequent in “p” but not “n”, and vice versa. Then, sort them in descending order.
  grammed_freq_sum <- grammed_p_freq + grammed_n_freq
  grammed_freq_diff <- abs(grammed_p_freq - grammed_n_freq)
  grammed_freq_percent <- grammed_freq_diff/grammed_freq_sum
  freq_percent_sorted <- sort(grammed_freq_percent, decreasing = TRUE)
  
  # Select words
  selected_grams <- names(freq_percent_sorted[freq_percent_sorted > threshold])
  return(selected_grams)
}





######### Machine Learning #########
# Calculate the brier score
brierScore <- function(forecastProb, actualResult) {
  summation <- sum((actualResult-forecastProb)^2)
  return(summation/length(actualResult))
}

# Classify reports as TP, TN, FP, or FN.
classification <- function(predicted_result, actual_result) {
  returnVal <- c()
  for(i in 1:length(predicted_result)) {
    if (predicted_result[i] == "p" && actual_result[i] == "p") {
      returnVal <- append(returnVal, "TP");
    } else if (predicted_result[i] == "n" && actual_result[i] == "p") {
      returnVal <- append(returnVal, "FN");
    } else if (predicted_result[i] == "p" && actual_result[i] == "n") {
      returnVal <- append(returnVal, "FP");
    } else if (predicted_result[i] == "n" && actual_result[i] == "n") {
      returnVal <- append(returnVal, "TN");
    }
  }
  return(returnVal); 
}

# Apply KNN algorithm 
knn_fxn <- function(data_subset, dtm, RANDOM_SEED){
  # Normalization - scale numerical features from 0 to 1 (not including target - species)
  normalize <- function (x) {
    return ( (x - min(x)) / (max(x) - min(x)))
  }
  dtm_normalized <- as.data.frame(lapply(dtm, normalize))
  
  # Train set, test set
  set.seed(RANDOM_SEED)
  index_row <- sample(2, size=nrow(dtm), replace=TRUE, prob=c(0.70,0.30))
  train_data <- dtm_normalized[index_row==1,]
  test_data <- dtm_normalized[index_row==2,]
  train_label <- data_subset$Pneumonia[index_row==1]
  test_label <- data_subset$Pneumonia[index_row==2]
  
  # Results
  predicted_result <- class::knn(train=train_data, test=test_data, cl=train_label, k = 5,prob=TRUE) # cl = class
  correct_result <- ifelse(test_label == predicted_result, 1, 0) # 0 = wrong result, 1 = correct result
  
  # Analysis
  TP <- table(test_label, predicted_result)[2,2]
  TN <- table(test_label, predicted_result)[1,1]
  FP <- table(test_label, predicted_result)[1,2]
  FN <- table(test_label, predicted_result)[2,1]
  sensitivity <- TP/(TP+FN)
  specificity  <- TN/(TN+FP)
  class <- classification(predicted_result, test_label)
  
  brier_score <- brierScore(attr(predicted_result,"prob"), ifelse(test_label == "p", 1, 0))
  
  
  multi_return <- list("Predicted_Result" = predicted_result, "Correct_Result" = correct_result, "Test_Label" = test_label,
                       "TP" = TP, "TN" = TN, "FP" = FP, "FN" = FN, "Sensitivity" = sensitivity, "Specificity" = specificity,
                       "Class" = class, "Brier_Score" = brier_score)
  return(multi_return)
}

# Apply random forest algorithm.
rf_fxn <- function(data_subset, dtm, RANDOM_SEED) {
  # Merge document-term matrix with Pneumonia data.
  dtm_with_pnuemonia <- dplyr::mutate(dtm, Pneumonia = data_subset$Pneumonia) # Check if merged correctly - dtm[,c("Pneumonia")]
  
  # Train set, test set
  set.seed(RANDOM_SEED)
  index_row <- sample(2, size=nrow(dtm), replace=TRUE, prob=c(0.70,0.30))
  train_data <- dtm_with_pnuemonia[index_row==1,]
  test_data <- dtm_with_pnuemonia[index_row==2,]
  train_label <- data_subset$Pneumonia[index_row==1]
  test_label <- data_subset$Pneumonia[index_row==2]
  
  # Converting variable names to all be legal (white space is not allowed)
  names(train_data) <- make.names(names(train_data))
  names(test_data) <- make.names(names(test_data))
  
  # Perform training
  pneumonia_classifier <- randomForest(Pneumonia~., data=train_data)
  
  # Results
  predicted_result <- predict(pneumonia_classifier, test_data[,-ncol(test_data)])
  correct_result <- ifelse(test_label == predicted_result, 1, 0)
  
  # Analysis
  TP <- table(test_label, predicted_result)[2,2]
  TN <- table(test_label, predicted_result)[1,1]
  FP <- table(test_label, predicted_result)[1,2]
  FN <- table(test_label, predicted_result)[2,1]
  sensitivity <- TP/(TP+FN)
  specificity  <- TN/(TN+FP)
  class <- classification(predicted_result, test_label)
  
  predicted_result_prob <- predict(pneumonia_classifier, test_data[,-ncol(test_data)],type="prob")[,2]
  brier_score <- brierScore(predicted_result_prob, ifelse(test_label == "p", 1, 0))
  
  multi_return <- list("Predicted_Result" = predicted_result, "Correct_Result" = correct_result, "train_data" = train_data, "pneumonia_classifier" = pneumonia_classifier,
                       "TP" = TP, "TN" = TN, "FP" = FP, "FN" = FN, "Sensitivity" = sensitivity, "Specificity" = specificity,
                       "Class" = class, "Brier_Score" = brier_score)
  return(multi_return)
}

# Apply gradient boost algorithm
gbm_fxn <- function(data_subset, dtm, RANDOM_SEED) {
  library(gbm)
  
  # Merge document-term matrix with Pneumonia data (0 or 1).
  dtm_with_pnuemonia <- dplyr::mutate(dtm, Pneumonia = ifelse(data_subset$Pneumonia == "p", 1, 0)) # Check if merged correctly - dtm[,c("Pneumonia")]
  
  # Train set, test set
  set.seed(RANDOM_SEED)
  index_row <- sample(2, size=nrow(dtm), replace=TRUE, prob=c(0.70,0.30))
  train_data <- dtm_with_pnuemonia[index_row==1,]
  test_data <- dtm_with_pnuemonia[index_row==2,]
  train_label <- data_subset$Pneumonia[index_row==1]
  test_label <- data_subset$Pneumonia[index_row==2]
  
  # Converting variable names to all be legal (white space is not allowed)
  names(train_data) <- make.names(names(train_data))
  names(test_data) <- make.names(names(test_data))
  
  # Perform training
  model = gbm(Pneumonia ~ ., data = train_data, distribution = "bernoulli", 
              n.trees = 5000, interaction.depth = 4, shrinkage = 0.01)
  
  # Results
  predicted_result_prob <- predict(model, test_data[,-ncol(test_data)], n.trees = 5000, type = "response")
  predicted_result <- ifelse(predicted_result_prob > 0.5, "p", "n")
  correct_result <- ifelse(test_label == predicted_result, 1, 0)
  
  # Analysis
  TP <- table(test_label, predicted_result)[2,2]
  TN <- table(test_label, predicted_result)[1,1]
  FP <- table(test_label, predicted_result)[1,2]
  FN <- table(test_label, predicted_result)[2,1]
  sensitivity <- TP/(TP+FN)
  specificity  <- TN/(TN+FP)
  class <- classification(predicted_result, test_label)
  
  brier_score <- brierScore(predicted_result_prob, ifelse(test_label == "p", 1, 0))
  
  multi_return <- list("Predicted_Result" = predicted_result, "Correct_Result" = correct_result, "Test_Label" = test_label, "predicted_result_prob" = predicted_result_prob,
                       "TP" = TP, "TN" = TN, "FP" = FP, "FN" = FN, "Sensitivity" = sensitivity, "Specificity" = specificity,
                       "Class" = class, "Brier_Score" = brier_score)
  return(multi_return)
}

# Apply xgboost algorithm
xgboost_fxn <- function(data_subset, dtm, RANDOM_SEED) {
  library(xgboost)
  
  # Train set, test set
  set.seed(RANDOM_SEED)
  index_row <- sample(2, size=nrow(dtm), replace=TRUE, prob=c(0.70,0.30))
  train_data <- dtm[index_row==1,]
  test_data <- dtm[index_row==2,]
  train_label <- data_subset$Pneumonia[index_row==1]
  test_label <- data_subset$Pneumonia[index_row==2]
  
  # Perform training
  params <- list(set.seed = RANDOM_SEED, eval_metric = "auc", objective = "binary:logistic")
  model <- xgboost(data = as.matrix(train_data), label = ifelse(train_label == "p", 1, 0), params = params, nrounds = 20, verbose = 1)
  
  # Results
  predicted_result_prob <- predict(model, as.matrix(test_data))
  predicted_result <- ifelse(predicted_result_prob > 0.5, "p", "n")
  correct_result <- ifelse(test_label == predicted_result, 1, 0)
  
  # Analysis
  TP <- table(test_label, predicted_result)[2,2]
  TN <- table(test_label, predicted_result)[1,1]
  FP <- table(test_label, predicted_result)[1,2]
  FN <- table(test_label, predicted_result)[2,1]
  sensitivity <- TP/(TP+FN)
  specificity  <- TN/(TN+FP)
  class <- classification(predicted_result, test_label)
  
  brier_score <- brierScore(predicted_result_prob, ifelse(test_label == "p", 1, 0))
  
  multi_return <- list("Predicted_Result" = predicted_result, "Correct_Result" = correct_result, "Test_Label" = test_label, "predicted_result_prob" = predicted_result_prob,
                       "TP" = TP, "TN" = TN, "FP" = FP, "FN" = FN, "Sensitivity" = sensitivity, "Specificity" = specificity,
                       "Class" = class, "Brier_Score" = brier_score)
  return(multi_return)
}

# Apply adaboost algorithm
adaboost_fxn <- function(data_subset, dtm, RANDOM_SEED) {
  library(caret)
  library(JOUSBoost)
  
  # Merge document-term matrix with Pneumonia data (0 or 1).
  dtm_with_pnuemonia <- dplyr::mutate(dtm, Pneumonia = ifelse(data_subset$Pneumonia == "p", 1, -1)) # Check if merged correctly - dtm[,c("Pneumonia")]
  
  # Train set, test set
  set.seed(RANDOM_SEED)
  index_row <- sample(2, size=nrow(dtm), replace=TRUE, prob=c(0.70,0.30))
  train_data <- dtm_with_pnuemonia[index_row==1,]
  test_data <- dtm_with_pnuemonia[index_row==2,]
  train_label <- data_subset$Pneumonia[index_row==1]
  test_label <- data_subset$Pneumonia[index_row==2]
  
  # Perform training
  model = adaboost(X = as.matrix(train_data[-ncol(train_data)]), y = train_data$Pneumonia, tree_depth = 4, n_rounds = 200, verbose = TRUE)
  options(scipen = 999)
  
  # Results
  predicted_result_prob <- predict(model, as.matrix(test_data[-ncol(train_data)]),type="prob")
  predicted_result <- ifelse(predicted_result_prob > 0.5, "p", "n")
  correct_result <- ifelse(test_label == predicted_result, 1, 0)
  
  # Analysis
  TP <- table(test_label, predicted_result)[2,2]
  TN <- table(test_label, predicted_result)[1,1]
  FP <- table(test_label, predicted_result)[1,2]
  FN <- table(test_label, predicted_result)[2,1]
  sensitivity <- TP/(TP+FN)
  specificity  <- TN/(TN+FP)
  class <- classification(predicted_result, test_label)
  
  brier_score <- brierScore(predicted_result_prob, ifelse(test_label == "p", 1, 0))
  
  multi_return <- list("Predicted_Result" = predicted_result, "Correct_Result" = correct_result, "Test_Label" = test_label, "predicted_result_prob" =predicted_result_prob,
                       "TP" = TP, "TN" = TN, "FP" = FP, "FN" = FN, "Sensitivity" = sensitivity, "Specificity" = specificity,
                       "Class" = class, "Brier_Score" = brier_score, "Test_Label" = test_label)
  return(multi_return)
}

# Run each of the machine learning algorithms and collect result data
collectData <- function(data_subset, dtm, RANDOM_SEED, num_reports) {
  set.seed(RANDOM_SEED)
  index_row <- sample(2, size=nrow(dtm), replace=TRUE, prob=c(0.70,0.30))
  
  ### KNN
  knn_result <- knn_fxn(data_subset, dtm, RANDOM_SEED)
  knn_count <- 0
  for(var in knn_result$Correct_Result){
    if(var == 1) {
      knn_count = knn_count + 1; 
    }
  }
  
  ### Random Forest
  rf_result <- rf_fxn(data_subset, dtm, RANDOM_SEED)
  rf_count <- 0
  for(var in rf_result$Correct_Result){
    if(var == 1) {
      rf_count = rf_count + 1; 
    }
  }
  
  
  ### GBM boost
  gbm_result <- gbm_fxn(data_subset, dtm, RANDOM_SEED)
  gbm_count <- 0
  for(var in gbm_result$Correct_Result){
    if(var == 1) {
      gbm_count = gbm_count + 1; 
    }
  }
  
  #### xgboost
  xgboost_result <- xgboost_fxn(data_subset, dtm, RANDOM_SEED)
  xgboost_count <- 0
  for(var in xgboost_result$Correct_Result){
    if(var == 1) {
      xgboost_count = xgboost_count + 1; 
    }
  }
  
  
  #### Adaboost
  adaboost_result <- adaboost_fxn(data_subset, dtm, RANDOM_SEED)
  adaboost_count <- 0
  for(var in adaboost_result$Correct_Result){
    if(var == 1) {
      adaboost_count = adaboost_count + 1; 
    }
  }
  
  # Collect Results
  patient_id <- dplyr::select(data_subset[index_row==2,], patient_id)
  report_id <- dplyr::select(data_subset[index_row==2,], report_id)
  path <- dplyr::select(data_subset[index_row==2,], path)
  
  df_results <- data.frame("Num_Reports" = num_reports, "RANDOME_SEED" = RANDOM_SEED, "Patient_ID" = patient_id, "Report_ID" = report_id, "Path" = path, 
                           "KNN" = knn_result$Correct_Result, "RandomForest" = rf_result$Correct_Result,"GbmBoost" = gbm_result$Correct_Result, "XGBoost" = xgboost_result$Correct_Result, "AdaBoost" = adaboost_result$Correct_Result,
                           "Class_KNN" = knn_result$Class, "Class_RandomForest" = rf_result$Class, "Class_GbmBoost" = gbm_result$Class, "Class_XGBoost" = xgboost_result$Class, "Class_AdaBoost" = adaboost_result$Class)
  df_correctResultCount <- data.frame("Num_Reports" = num_reports, "RANDOME_SEED" = RANDOM_SEED, KNN_Count = knn_count, RandomForest_Count = rf_count, Gbm_Count = gbm_count, XGBoost_Count = xgboost_count, AdaBoost_Count = adaboost_count)
  df_brierScore <- data.frame("Num_Reports" = num_reports, "RANDOME_SEED" = RANDOM_SEED, "KNN_BrierScore" = knn_result$Brier_Score, "RandomForest_BrierScore" = rf_result$Brier_Score, "GbmBoost_BrierScore" = gbm_result$Brier_Score, "XGBoost_BrierScore" = xgboost_result$Brier_Score, "AdaBoost_BrierScore" = adaboost_result$Brier_Score)
  
  # AUC
  df_AUC <- data.frame("Num_Reports" = num_reports, "RANDOME_SEED" = RANDOM_SEED, 
                       "KNN_AUC" = auc(knn_result$Test_Label, attributes(knn_result$Predicted_Result)$prob), 
                       "RandomForest_AUC" = auc((rf_result$train_data)$Pneumonia, (rf_result$pneumonia_classifier)$votes[,1]), 
                       "GbmBoost_AUC" = auc(gbm_result$Test_Label, gbm_result$predicted_result_prob), 
                       "XGBoost_AUC" = auc(xgboost_result$Test_Label, xgboost_result$predicted_result_prob), 
                       "AdaBoost_AUC" = auc(adaboost_result$Test_Label, adaboost_result$predicted_result_prob))
  
  
  # Contingency Table
  df_contingencyTable_knn <- data.frame("Num_Reports" = num_reports, "RANDOME_SEED" = RANDOM_SEED, "Model" = "KNN",
                                        "TP" = knn_result$TP, "TN" = knn_result$TN, "FP" = knn_result$FP, "FN" = knn_result$FN, 
                                        "Sensitivity" = knn_result$Sensitivity, "Specificity" = knn_result$Specificity)
  df_contingencyTable_rf <- data.frame("Num_Reports" = num_reports, "RANDOME_SEED" = RANDOM_SEED, "Model" = "Random Forest",
                                       "TP" = rf_result$TP, "TN" = rf_result$TN, "FP" = rf_result$FP, "FN" = rf_result$FN, 
                                       "Sensitivity" = rf_result$Sensitivity, "Specificity" = rf_result$Specificity)
  df_contingencyTable_gbm <- data.frame("Num_Reports" = num_reports, "RANDOME_SEED" = RANDOM_SEED, "Model" = "GbmBoost",
                                        "TP" = gbm_result$TP, "TN" = gbm_result$TN, "FP" = gbm_result$FP, "FN" = gbm_result$FN, 
                                        "Sensitivity" = gbm_result$Sensitivity, "Specificity" = gbm_result$Specificity)
  df_contingencyTable_xgboost <- data.frame("Num_Reports" = num_reports, "RANDOME_SEED" = RANDOM_SEED, "Model" = "XGBoost",
                                            "TP" = xgboost_result$TP, "TN" = xgboost_result$TN, "FP" = xgboost_result$FP, "FN" = xgboost_result$FN, 
                                            "Sensitivity" = xgboost_result$Sensitivity, "Specificity" = xgboost_result$Specificity)
  df_contingencyTable_adaboost <- data.frame("Num_Reports" = num_reports, "RANDOME_SEED" = RANDOM_SEED, "Model" = "AdaBoost",
                                             "TP" = adaboost_result$TP, "TN" = adaboost_result$TN, "FP" = adaboost_result$FP, "FN" = adaboost_result$FN, 
                                             "Sensitivity" = adaboost_result$Sensitivity, "Specificity" = adaboost_result$Specificity)
  df_contingencyTable <- dplyr::bind_rows(df_contingencyTable_knn, df_contingencyTable_rf, df_contingencyTable_gbm, df_contingencyTable_xgboost, df_contingencyTable_adaboost)
  
  multi_return <- list("Results" = df_results, "Correct_Results_Count" = df_correctResultCount, "Contingency_Table" = df_contingencyTable, "Brier_Scores" = df_brierScore, "AUC" = df_AUC)
  #multi_return <- dplyr::mutate(df_results, df_correctResultCount, df_brierScore)
  
  return(multi_return)
  
  #writexl::write_xlsx(ResultsDataFrame, "Results.xlsx")
  
}


# All the combinations of parameters that are possible.
param_tibble <- expand.grid(RANDOM_SEED = c(123456789, 987654321, 192837465, 918273645, 111111111, 222222222), 
                            num_reports = c(500, 1000, 2000, 3000))
#out <- lapply(
#  1L:nrow(param_tibble),
#  function(i){
    data_subset <- loadDataSubset(123456789, 1000)
    dtm_wide <-  createCorpus(data_subset) %>%
      createWideDTM()
    selected_grams <- featureSelection(data_subset, dtm_wide, 0.5)
    remove_terms <- c("ap view chest",
                      "ap later view",
                      "lung volum low",
                      "opac may reflect",
                      "opac may repres",
                      "overt pulmonari edema",
                      "upright ap view",
                      "opac lung base",
                      "correct clinic set",
                      "left basilar opac",
                      "consolid concern pneumonia",
                      "portable view chest",
                      "compar prior studi")
    selected_grams <- selected_grams [! selected_grams %in% remove_terms]
    dtm <- dtm_wide[,selected_grams]
    result <- adaboost_fxn(data_subset, dtm, 123456789)
    
    #KNN
    roc(result$Test_Label, result$predicted_result_prob, legacy.axes=TRUE, percent=TRUE,
        xlab="False Positive Percentage", ylab="True Positive Percentage", col="#4daf4a", lwd=4, print.auc=TRUE)
    
    result <- collectData(data_subset, dtm, param_tibble$RANDOM_SEED[i], param_tibble$num_reports[i])
    return(result)
  }
#)

# Run all the combinations
out <- lapply(
  1L:nrow(param_tibble),
  function(i){
    data_subset <- loadDataSubset(param_tibble$RANDOM_SEED[i], param_tibble$num_reports[i])
    dtm_wide <- createCorpus(data_subset) %>%
      createWideDTM()
    selected_grams <- featureSelection(data_subset, dtm_wide, 0.5)
    dtm <- dtm_wide[,selected_grams]
    result <- collectData(data_subset, dtm, param_tibble$RANDOM_SEED[i], param_tibble$num_reports[i])
    return(result)
  }
)

################### Analyze Result #################
df_results <- NULL
df_contingency_table <- NULL
df_brier_scores <- NULL
df_auc <- NULL
df_summarized_contingency_table <- NULL
df_summarized_brier_scores <- NULL

for(i in 1:length(out)) {
  df_results <- dplyr::bind_rows(df_results, out[[i]]$Results) 
  df_contingency_table <- dplyr::bind_rows(df_contingency_table, out[[i]]$Contingency_Table) 
  df_brier_scores <- dplyr::bind_rows(df_brier_scores, out[[i]]$Brier_Scores) 
  df_auc <- rbind(df_auc, out[[i]]$AUC) 
}

xlsx::write.xlsx(df_results, file="Results.xlsx", row.names=FALSE) 

#### Contingency Table Analysis #####
contingency_table_analysis <- function(df_contingency_table, num_reports) {
  df_contingency_table_knn <- filter(df_contingency_table, Num_Reports==num_reports, Model == "KNN")
  df_contingency_table_rf <- filter(df_contingency_table, Num_Reports==num_reports, Model == "Random Forest")
  df_contingency_table_gbm <- filter(df_contingency_table, Num_Reports==num_reports, Model == "GbmBoost")
  df_contingency_table_xgboost <- filter(df_contingency_table, Num_Reports==num_reports, Model == "XGBoost")
  df_contingency_table_adaboost <- filter(df_contingency_table, Num_Reports==num_reports, Model == "AdaBoost")
  
  df_summarized_contingency_table_knn <- dplyr::summarise(df_contingency_table_knn, mean(TP), mean(TN), mean(FP), mean(FN), mean(Sensitivity), mean(Specificity))
  df_summarized_contingency_table_knn <- cbind(Num_Reports = num_reports, Model = "KNN", df_summarized_contingency_table_knn)
  
  df_summarized_contingency_table_rf <- dplyr::summarise(df_contingency_table_rf, mean(TP), mean(TN), mean(FP), mean(FN), mean(Sensitivity), mean(Specificity))
  df_summarized_contingency_table_rf <- cbind(Num_Reports = num_reports, Model = "Random Forest", df_summarized_contingency_table_rf)
  
  df_summarized_contingency_table_gbm <- dplyr::summarise(df_contingency_table_gbm, mean(TP), mean(TN), mean(FP), mean(FN), mean(Sensitivity), mean(Specificity))
  df_summarized_contingency_table_gbm <- cbind(Num_Reports = num_reports, Model = "GbmBoost", df_summarized_contingency_table_gbm)
  
  df_summarized_contingency_table_xgboost <- dplyr::summarise(df_contingency_table_xgboost, mean(TP), mean(TN), mean(FP), mean(FN), mean(Sensitivity), mean(Specificity))
  df_summarized_contingency_table_xgboost <- cbind(Num_Reports = num_reports, Model = "XGBoost", df_summarized_contingency_table_xgboost)
  
  df_summarized_contingency_table_adaboost <- dplyr::summarise(df_contingency_table_adaboost, mean(TP), mean(TN), mean(FP), mean(FN), mean(Sensitivity), mean(Specificity))
  df_summarized_contingency_table_adaboost <- cbind(Num_Reports = num_reports, Model = "AdaBoost", df_summarized_contingency_table_adaboost)
  
  df_summarized_contingency_table <- dplyr::bind_rows(df_summarized_contingency_table_knn, df_summarized_contingency_table_rf, df_summarized_contingency_table_gbm, 
                                                      df_summarized_contingency_table_xgboost, df_summarized_contingency_table_adaboost) 
}

df_summarized_contingency_table_500 <- contingency_table_analysis(df_contingency_table, 500)
df_summarized_contingency_table_1000 <- contingency_table_analysis(df_contingency_table, 1000)
df_summarized_contingency_table_2000 <- contingency_table_analysis(df_contingency_table, 2000)
df_summarized_contingency_table_3000 <- contingency_table_analysis(df_contingency_table, 3000)
df_summarized_contingency_table <- dplyr::bind_rows(df_summarized_contingency_table_500, df_summarized_contingency_table_1000, df_summarized_contingency_table_2000, df_summarized_contingency_table_3000)
df_summarized_contingency_table <- df_summarized_contingency_table %>% 
  rename(
    Mean.TP = "mean(TP)",
    Mean.TN = "mean(TN)",
    Mean.FP = "mean(FP)",
    Mean.FN = "mean(FN)",
    Mean.Sensitivity = "mean(Sensitivity)",
    Mean.Specificity = "mean(Specificity)"
  )

xlsx::write.xlsx(df_summarized_contingency_table, file="Contingency_Table.xlsx", row.names=FALSE) 

#### Brier Scores Analysis #####
brier_score_analysis <- function (df_brier_scores, num_reports) {
  df_filtered_brier_scores <- filter(df_brier_scores, Num_Reports==num_reports)
  df_summarized_brier_scores <- dplyr::summarise(df_filtered_brier_scores, mean(KNN_BrierScore), mean(RandomForest_BrierScore), mean(GbmBoost_BrierScore), mean(XGBoost_BrierScore), mean(AdaBoost_BrierScore))
  df_summarized_brier_scores <- cbind(Num_Reports = num_reports, df_summarized_brier_scores)
}

df_summarized_brier_scores_500 <- brier_score_analysis(df_brier_scores, 500)
df_summarized_brier_scores_1000 <- brier_score_analysis(df_brier_scores, 1000)
df_summarized_brier_scores_2000 <- brier_score_analysis(df_brier_scores, 2000)
df_summarized_brier_scores_3000 <- brier_score_analysis(df_brier_scores, 3000)
df_summarized_brier_scores <- dplyr::bind_rows(df_summarized_brier_scores_500, df_summarized_brier_scores_1000, df_summarized_brier_scores_2000, df_summarized_brier_scores_3000)
df_summarized_brier_scores <- df_summarized_brier_scores %>% 
  rename(
    Mean.KNN_BrierScore = "mean(KNN_BrierScore)",
    Mean.RandomForest_BrierScore = "mean(RandomForest_BrierScore)",
    Mean.GbmBoost_BrierScore = "mean(GbmBoost_BrierScore)",
    Mean.XGBoost_BrierScore = "mean(XGBoost_BrierScore)",
    Mean.AdaBoost_BrierScore = "mean(AdaBoost_BrierScore)"
  )
df_summarized_brier_scores_final <- df_summarized_brier_scores[,-1]
rownames(df_summarized_brier_scores_final) <- df_summarized_brier_scores[,1]

xlsx::write.xlsx(df_summarized_brier_scores, file="Brier_Scores.xlsx", row.names=FALSE) 

#### AUC Analysis #####
auc_analysis <- function (df_auc, num_reports) {
  df_filtered_auc <- filter(df_auc, Num_Reports==num_reports)
  df_summarized_auc <- dplyr::summarise(df_filtered_auc, mean(KNN_AUC), mean(RandomForest_AUC), mean(GbmBoost_AUC), mean(XGBoost_AUC), mean(AdaBoost_AUC))
  #df_summarized_auc_sd <- dplyr::summarise(df_filtered_auc, sd(KNN_AUC), sd(RandomForest_AUC), sd(GbmBoost_AUC), sd(XGBoost_AUC), sd(AdaBoost_AUC))
  df_summarized_auc <- cbind(Num_Reports = num_reports, df_summarized_auc)
}

df_summarized_auc_500 <- auc_analysis(df_auc, 500)
df_summarized_auc_1000 <- auc_analysis(df_auc, 1000)
df_summarized_auc_2000 <- auc_analysis(df_auc, 2000)
df_summarized_auc_3000 <- auc_analysis(df_auc, 3000)
df_summarized_auc <- dplyr::bind_rows(df_summarized_auc_500, df_summarized_auc_1000, df_summarized_auc_2000, df_summarized_auc_3000)
df_summarized_auc <- df_summarized_auc %>% 
  rename(
    Mean.KNN_AUC = "mean(KNN_AUC)",
    Mean.RandomForest_AUC = "mean(RandomForest_AUC)",
    Mean.GbmBoost_AUC = "mean(GbmBoost_AUC)",
    Mean.XGBoost_AUC = "mean(XGBoost_AUC)",
    Mean.AdaBoost_AUC = "mean(AdaBoost_AUC)"
  )
df_summarized_auc_final <- df_summarized_auc[,-1]
rownames(df_summarized_auc_final) <- df_summarized_auc[,1]

xlsx::write.xlsx(df_summarized_auc, file="AUC.xlsx", row.names=FALSE) 


#### Graphs ####
library(reshape2)
ggplot(df_summarized_contingency_table, aes(x=Num_Reports, y=Mean.Sensitivity)) + geom_bar(stat="identity") + facet_grid(.~Model)

# Brier Scores
df_summarized_brier_scores_final$num_reports  <- row.names(df_summarized_brier_scores_final)
df_summarized_brier_scores_final$num_reports <- factor(df_summarized_brier_scores_final$num_reports,levels = c("500", "1000", "2000", "3000"))
df_summarized_brier_scores_final.molten <- melt(df_summarized_brier_scores_final, value.name="Brier_Scores", variable.name="Model", na.rm=TRUE)
ggplot(df_summarized_brier_scores_final.molten, aes(x=Model, y=Brier_Scores, fill=num_reports)) + geom_bar(stat="identity", position="dodge") + coord_flip()

# AUC
df_summarized_auc_final$num_reports  <- row.names(df_summarized_auc_final)
df_summarized_auc_final$num_reports <- factor(df_summarized_auc_final$num_reports,levels = c("500", "1000", "2000", "3000"))
df_summarized_auc_final.molten <- melt(df_summarized_auc_final, value.name="AUC", variable.name="Model", na.rm=TRUE)
ggplot(df_summarized_auc_final.molten, aes(x=Model, y=AUC, fill=num_reports)) + geom_bar(stat="identity", position="dodge") + coord_flip()



d=data.frame(drink=c("coffee","tea","water"), mean=c(3,6,2), lower=c(2.6,5.6,1.8), upper=c(3.5,6.3,2.8))
ggplot() + 
  geom_errorbarh(data=d, mapping=aes(y=drink, x=upper, xmin=upper, xmax=lower), height=0.2, size=1, color="blue") + 
  geom_point(data=d, mapping=aes(y=drink, x=mean), size=4, shape=21, fill="white") 


df_auc
df_auc_final <- df_auc[1:24,-1]
df_auc_final <- df_auc_final[,-1]
#rownames(df_auc_final) <- df_auc[1:24,1]
df_auc_final$num_reports <- df_auc[1:24,1]
df_auc_final$num_reports <- factor(df_auc_final$num_reports,levels = c("500", "1000", "2000", "3000"))
df_auc_final.molten <- melt(df_auc_final, id="num_reports",value.name="AUC", variable.name="Model", na.rm=TRUE)
library(Hmisc)
ggplot(df_auc_final.molten, aes(x=Model, y=AUC, fill=num_reports)) + facet_wrap(~ num_reports) +
  stat_summary(geom="bar", fun = "mean") + stat_summary(geom = "errorbar", fun.data = "mean_cl_normal") + coord_flip()


df_brier_scores
df_brier_scores_final <- df_brier_scores[1:24,-1]
df_brier_scores_final <- df_brier_scores_final[,-1]
#rownames(df_auc_final) <- df_auc[1:24,1]
df_brier_scores_final$num_reports <- df_brier_scores[1:24,1]
df_brier_scores_final$num_reports <- factor(df_brier_scores_final$num_reports,levels = c("500", "1000", "2000", "3000"))
df_brier_scores_final.molten <- melt(df_brier_scores_final, id="num_reports",value.name="Brier_Score", variable.name="Model", na.rm=TRUE)
library(Hmisc)
ggplot(df_brier_scores_final.molten, aes(x=Model, y=Brier_Score, fill=num_reports)) + facet_wrap(~ num_reports) +
  stat_summary(geom="bar", fun = "mean") + stat_summary(geom = "errorbar", fun.data = "mean_cl_normal") + coord_flip()



######### ROC curve #########
rocCurve <- function(){
  
  data_subset <- loadDataSubset(123456789, 3000)
  dtm_wide <-  createCorpus(data_subset) %>% createWideDTM()
  selected_grams <- featureSelection(data_subset, dtm_wide)
  dtm <- dtm_wide[,selected_grams]
  
  # Train set, test set
  set.seed(123456789)
  index_row <- sample(2, size=nrow(dtm), replace=TRUE, prob=c(0.70,0.30))
  train_data <- dtm[index_row==1,]
  test_data <- dtm[index_row==2,]
  train_label <- data_subset$Pneumonia[index_row==1]
  test_label <- data_subset$Pneumonia[index_row==2]
  names(train_data) <- make.names(names(train_data))
  ### KNN
  knn_result <- knn_fxn(data_subset, dtm, 123456789)
  
  ### Random Forest
  rf_result <- rf_fxn(data_subset, dtm, 123456789)
  
  ### GBM boost
  gbm_result <- gbm_fxn(data_subset, dtm, 123456789)
  
  #### xgboost
  xgboost_result <- xgboost_fxn(data_subset, dtm, 123456789)
  
  #### Adaboost
  adaboost_result <- adaboost_fxn(data_subset, dtm, 123456789)
  
  # ROC CURVE
  par(pty="s")
  # dtm_test_label - pass in the known classification, positive or negative, for each sample.
  # attributes(prediction)$prob - esitemated probabilites that each sample is positive
  # plot=TRUE - tell the roc() function to draw the graph, not just calculate all of the numbers used to draw the graph.
  
  #KNN
  roc(test_label, attributes(knn_result$Predicted_Result)$prob, plot=TRUE, legacy.axes=TRUE, percent=TRUE,
      xlab="False Positive Percentage", ylab="True Positive Percentage", col="#4daf4a", lwd=4, print.auc=TRUE)
  #Radom Forest
  plot.roc((rf_result$train_data)$Pneumonia, (rf_result$pneumonia_classifier)$votes[,1], percent=TRUE, col="#377eb8", lwd=4, 
           print.auc=TRUE, add=TRUE, print.auc.y=40)
  # gradient boosting
  plot.roc(test_label, gbm_result$predicted_result_prob, percent=TRUE, col="Orange", lwd=4, 
           print.auc=TRUE, add=TRUE, print.auc.y=45)
  # xgboost boosting
  plot.roc(test_label, xgboost_result$predicted_result_prob, percent=TRUE, col="Violet", lwd=4, 
           print.auc=TRUE, add=TRUE, print.auc.y=35)
  # Adaboost boosting
  plot.roc(test_label, adaboost_result$predicted_result_prob, percent=TRUE, col="Grey", lwd=4, 
           print.auc=TRUE, add=TRUE, print.auc.y=30)
  legend("bottomright", legend=c("KNN", "Random Forest", "Gradient Boosting", "Xgboost", "Adaboost"), col=c("#4daf4a", "#377eb8", "Orange", "Violet", "Grey"), lwd=4)
}
rocCurve()


library(PRROC)
###
test_label<-ifelse(test_label == "p", 1, 0)

knn_scores <- data.frame("Predictions" = attributes(knn_result$Predicted_Result)$prob, "Labels" = test_label)
knn_pr <- pr.curve(scores.class0=knn_scores[knn_scores$Labels=="1",]$Predictions,
                   scores.class1=knn_scores[knn_scores$Labels=="0",]$Predictions,
                   curve=T)
str(knn_pr)
head(knn_pr$curve)
plot(knn_pr)

###
rf_test_label<-ifelse((rf_result$train_data)$Pneumonia == "p", 1, 0)
rf_scores <- data.frame("Predictions" = (rf_result$pneumonia_classifier)$votes[,1], "Labels" = rf_test_label)
rf_pr <- pr.curve(scores.class0=rf_scores[rf_scores$Labels=="1",]$Predictions,
                  scores.class1=rf_scores[rf_scores$Labels=="0",]$Predictions,
                  curve=T)


###
gbm_scores <- data.frame("Predictions" = gbm_result$predicted_result_prob, "Labels" = test_label)
gbm_pr <- pr.curve(scores.class0=gbm_scores[gbm_scores$Labels=="1",]$Predictions,
                   scores.class1=gbm_scores[gbm_scores$Labels=="0",]$Predictions,
                   curve=T)
plot(gbm_pr)

###
xgboost_scores <- data.frame("Predictions" = xgboost_result$predicted_result_prob, "Labels" = test_label)
xgboost_pr <- pr.curve(scores.class0=xgboost_scores[xgboost_scores$Labels=="1",]$Predictions,
                       scores.class1=xgboost_scores[xgboost_scores$Labels=="0",]$Predictions,
                       curve=T)
plot(xgboost_pr)

###
adaboost_scores <- data.frame("Predictions" = adaboost_result$predicted_result_prob, "Labels" = test_label)
adaboost_pr <- pr.curve(scores.class0=adaboost_scores[adaboost_scores$Labels=="1",]$Predictions,
                        scores.class1=adaboost_scores[adaboost_scores$Labels=="0",]$Predictions,
                        curve=T)
plot(adaboost_pr)

