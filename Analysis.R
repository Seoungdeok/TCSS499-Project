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

# Set parameters.
RANDOM_SEED <- 123456789
data_folder <- "data" # Folder containing file_map.csv and mimic-cxr-2.0.0-negbio.csv
file_map_path <- file.path(data_folder, "file_map.csv")
report_labels_path <- file.path(data_folder, "mimic-cxr-2.0.0-negbio.csv.gz")
num_reports <- 1000 # Reports to load, max of 227835 (total number of reports).

# Create a results directory.
if(dir.exists("results")){
  unlink("results", recursive = TRUE, force = TRUE)
}
dir.create("results")

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

# Look at the first few rows of each data set.
head(file_map)
nrow(file_map)

head(report_labels)
nrow(report_labels)

# Merge the data sets.
merged <- dplyr::inner_join(
  file_map, report_labels,
  by = c("patient_id", "report_id")
)
head(merged)
nrow(merged)


# Set the random seed, then randomly select some reports.
set.seed(RANDOM_SEED)
merged_subset <- dplyr::filter(merged, !is.na(Pneumonia)) # Select only reports for which Pneumonia is not NA.
merged_subset <- dplyr::filter(merged_subset, Pneumonia == "n" | Pneumonia == "p")  # Select only reports for which Pneumonia is "p" or "n".
merged_subset$Pneumonia <- factor(merged_subset$Pneumonia) # Drop "u" level from Pneumonia factor.
str(merged_subset) 
keepers <- sample(1L:nrow(merged_subset), min(num_reports, nrow(merged_subset))) # Select 1000 documents randomly
merged_subset <- merged_subset[keepers,]

# Create a vector of the reports.
reports_subset <- merged_subset$report
head(reports_subset)

######### Apply tm package ####################################################
# Parameters for the analysis.
stem_words <- TRUE # Should words be stemmed or not.
remove_stopwords <- TRUE # Should stop words be removed or not.
stop_words <- c(# Comma separated list of words or phrases that should be removed.
  stopwords("en"),
  "final report"
)

gram_size <- 2:3 # Gram size. You can set this to a single integer or a range, e.g. 1:4

# Create and clean up the corpus.
my_corpus <- VCorpus(VectorSource(reports_subset)) %>%
  tm_map(tolower) %>% # Make all lower case.
  tm_map(function(x){gsub("[[:digit:]]","",x)}) %>% # Remove numbers
  tm_map(function(x){gsub("'"," ",x)}) %>% # Apostrophes to spaces.
  tm_map(function(x){gsub("_+[[:alpha:]]+","_",x)}) %>% # Remove PII indicators.
  tm_map(removePunctuation) %>% # Remove punctuation, including _
  (function(x){# Remove stop words if remove_stopwords is TRUE.
    if(remove_stopwords){
      return(tm_map(x, removeWords, stop_words))
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

# Create a document by n-gram count table.
grammed_wide <- 1L:length(my_corpus) %>% # For each document in the corpus...
  lapply(
    .,
    function(j){
      x <- unlist(my_corpus[j]) # Retrieve the document and convert to a string.
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
grammed_wide
# Create a tall version of the data too, mainly useful for plotting purposes.
grammed_tall <- grammed_wide %>%
  tidyr::pivot_longer(-one_of("DOCUMENT_ID"), names_to = "gram", values_to = "count")
grammed_tall

# Count the number of reports that contain each gram.
doc_gram_presence <- grammed_tall %>%
  dplyr::mutate(
    presence = factor(
      count > 0L,
      levels = c(FALSE, TRUE),
      labels = c("absent", "present")
    )
  ) %>%
  dplyr::group_by(gram) %>%
  dplyr::summarise(
    present_count = sum(presence == "present"),
    present_percent = mean(presence == "present")*100,
    .groups = "drop"
  ) %>%
  dplyr::arrange(present_count)

# Create a histogram that depicts the distribution of gram usage across reports.
# Most of the counts for this histogram are on the left, indicating that most
# grams are used in only a small number of reports. However, there are some
# grams that are found in most reports.
plt <- ggplot()+
  geom_histogram(
    data = doc_gram_presence,
    mapping = aes(present_percent-1/(2*num_reports)*100),
    breaks = seq(0,100,2),
    fill = "black", colour = "white"
  )+
  xlab("Reports (%)")+ylab("Count")+
  coord_cartesian(xlim = c(0,100))+
  scale_x_continuous(breaks = seq(0,100,10), expand = expansion(0,0))+
  theme(
    plot.margin = margin(t = 1, r = 1, unit = "cm")
  )+
  scale_y_log10(limits = c(1,NA), expand = expansion(c(0,0), c(0,0.5)))

cowplot::save_plot(
  file.path("results","Histogram_GramCounts.tiff"),
  plt,
  base_width = 6.5, base_height = 4
)

# Next, we will break this down by disease status.
# First, merge data sets.
grammed_tall_merged <- dplyr::inner_join(
  merged_subset %>%
    dplyr::mutate(DOCUMENT_ID = dplyr::row_number()),
  grammed_tall,
  by = "DOCUMENT_ID"
)
grammed_tall_merged

# Next, select how many of the most common grams should be used in plotting. I
# recommend keeping this relatively small.
num_top <- 20
top_n_grams <- tail(doc_gram_presence, num_top)$gram # Retrieve the top grams. Using tail() because it's in ascending order.

# Plot a heatmap that depicts how many times each of the top grams appears in
# each report.
plt <- ggplot()+
  geom_tile(
    data = grammed_tall_merged %>%
      dplyr::filter(gram %in% top_n_grams) %>%
      dplyr::mutate(
        gram = factor(gram, levels = top_n_grams),
        DOCUMENT_ID = factor(DOCUMENT_ID)
      ),
    mapping = aes(DOCUMENT_ID, gram, fill = count)
  )+
  scale_fill_gradient(low = "white", high = "black")+
  xlab("Report")+labs(fill = "Count")+
  scale_x_discrete(expand = expansion(0,0))+
  scale_y_discrete(expand = expansion(0,0))+
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )+
  facet_grid(cols = vars(Pneumonia), scales = "free_x", space = "free_x")

cowplot::save_plot(
  file.path("results","Heatmap_Pneumonia.tiff"),
  plt,
  base_width = 6.5, base_height = 4
)

# Next, we will plot the percent of reports that contain each gram. Again, we
# will only consider the top grams.
# First, merge data sets and calculate percentages.
merged_doc_gram_presence <- grammed_tall_merged %>%
  dplyr::mutate(
    presence = factor(
      count > 0L,
      levels = c(FALSE, TRUE),
      labels = c("absent", "present")
    )
  ) %>%
  dplyr::group_by(gram, Pneumonia) %>%
  dplyr::summarise(
    present_count = sum(presence == "present"),
    present_percent = mean(presence == "present")*100,
    .groups = "drop"
  ) %>%
  tidyr::drop_na()

plt <- ggplot()+
  geom_col(
    data = merged_doc_gram_presence %>%# Filter to just the top grams.
      dplyr::filter(gram %in% top_n_grams) %>%
      dplyr::mutate(
        gram = factor(gram, levels = top_n_grams)
      ),
    mapping = aes(present_percent, gram),
    fill = "black"# Bar color.
  )+
  facet_wrap(~Pneumonia, nrow = 1)+# Create a panel for each pneumonia group.
  coord_cartesian(xlim = c(0,100))+# Make the axes pretty.
  scale_x_continuous(expand = expansion(0,0))+
  scale_y_discrete(expand = expansion(0,0))+
  theme(# Add grid lines and change spacing.
    panel.grid.major.x = element_line(colour = "gray"),
    panel.spacing.x = unit(1.2, "lines"),
    plot.margin = margin(r = 1, unit = "cm")
  )+
  panel_border(colour = "black")+# Add a border to each panel.
  xlab("Reports (%)")# X axis label.

cowplot::save_plot(
  file.path("results","Barchart_Pneumonia.tiff"),
  plt,
  base_width = 6.5, base_height = 6
)

# install.packages("ggwordcloud") Finally, we will create a word cloud. First,
# customize the tibble created above by adding to each row random rotations (90
# degree or 0 degree rotation) and 1 of 7 random colors (generated via the
# rainbow function).
temp <- merged_doc_gram_presence %>%
  dplyr::filter(gram %in% top_n_grams) %>%
  dplyr::mutate(
    gram = factor(gram, levels = top_n_grams),
    rotation = 90*sample(c(0,1), dplyr::n(), replace = TRUE),
    colour = sample(rainbow(7), dplyr::n(), replace = TRUE)
  )
plt <- ggplot()+
  ggwordcloud::geom_text_wordcloud_area(
    data = temp,
    aes(
      label = gram,
      size = present_percent,
      angle = rotation,
      colour = colour
    )
  )+
  facet_wrap(~Pneumonia, nrow = 1)+# Panel for each pneumonia group.
  theme_minimal()

cowplot::save_plot(
  file.path("results","Wordcloud_Pneumonia.tiff"),
  plt,
  base_width = 6.5, base_height = 4
)

# Feature selection.
#
# I am not sure if this is the best approach. This is what I do with
# RNA-sequencing data and I think it should work here. We should explore other
# options.
# To install edgeR...
install.packages('BiocManager')
BiocManager::install('edgeR')

np_rows <- merged_subset$Pneumonia %in% c("n","p")
counts <- grammed_wide[np_rows,] %>%
  dplyr::select(-DOCUMENT_ID) %>%
  as.matrix()
pneumonia_status <- factor(
  as.character(merged_subset$Pneumonia[np_rows]),
  levels = c("n","p")
)
design <- model.matrix(~pneumonia_status)
dispersion <- edgeR::estimateCommonDisp(t(counts))
fit <- edgeR::glmFit(
  t(counts),
  design = design,
  robust = TRUE,
  dispersion = dispersion,
  abundance.trend = FALSE
)
test_of_fit <- edgeR::glmLRT(fit)
test_results <- edgeR::topTags(test_of_fit, Inf)

# Put into a consistent format.
test_results <- tibble::rownames_to_column(
  as.data.frame(test_results), "Gram"
)
test_results <- tibble::as_tibble(test_results)
test_results <- dplyr::select(
  test_results,
  -logCPM, -LR
)




##################KNN###################
# Returns one row tibble. Columns of this going to be K, Threshold, TP, TN, FP, FN, Sensitivity, Specificity.
knn_fxn <- function(K, Threshold){
  # Code for running the KNN model.
  # Select words that appear at least 10 times.
  grammed_popular <- grammed_wide[,colSums(as.matrix(grammed_wide)) > 10]
  
  # Merge with Pneumonia data.
  pneumonia_data <- merged_subset$Pneumonia
  grammed_popular <- dplyr::mutate(grammed_popular, Pneumonia = pneumonia_data)
  
  # Split by Pneumonia == "p" and Pneumonia == "n".
  grammed_p <- dplyr::filter(grammed_popular, Pneumonia == "p")
  grammed_n <- dplyr::filter(grammed_popular, Pneumonia == "n")
  
  # Remove columns called "DOCUMENT_ID" and "Pnemonia" and sum up the columns.
  grammed_p_freq <- colSums(as.matrix(grammed_p[,c(-1,-ncol(grammed_p))]))
  grammed_n_freq <- colSums(as.matrix(grammed_n[,c(-1,-ncol(grammed_n))]))
  
  # Identify terms that are frequent in “p” but not “n”, and vice versa. Then, sort them in descending order.
  grammed_freq_sum <- grammed_p_freq + grammed_n_freq
  grammed_freq_diff <- abs(grammed_p_freq - grammed_n_freq)
  grammed_freq_percent <- grammed_freq_diff/grammed_freq_sum
  freq_percent_sorted <- sort(grammed_freq_percent, decreasing = TRUE)
  
  # Select words
  selected_words <- names(freq_percent_sorted[freq_percent_sorted > Threshold])
  dtm <- grammed_popular[,selected_words]
  
  # normalization - scale numerical features from 0 to 1 (not including target - species)
  normalize <- function (x) {
    return ( (x - min(x)) / (max(x) - min(x)))
  }
  dtm_normalized <- as.data.frame(lapply(dtm, normalize))
  
  # create training data set (80%) & test data set (20%)
  set.seed(10)
  num_rows <- nrow(dtm)
  num_train <- round(num_rows * 0.8)
  num_test <- num_rows - num_train
  index <- sample(1:num_rows, num_train, replace=F)
  dtm_train <- dtm_normalized[index, ]
  dtm_test <- dtm_normalized[-index, ]
  dtm_train_label <- pneumonia_data[index]
  dtm_test_label <- pneumonia_data[-index]
  
  table(dtm_train_label)
  table(dtm_test_label)
  
  # apply knn algorithm
  prediction <- class::knn(train=dtm_train, test=dtm_test, cl=dtm_train_label, k = K,prob=TRUE) # cl = class & k stands for how many nearest neighbor you want
  
  # TP, TN, FP, FN
  TP <- table(dtm_test_label, prediction)[2,2]
  TN <- table(dtm_test_label, prediction)[1,1]
  FP <- table(dtm_test_label, prediction)[1,2]
  FN <- table(dtm_test_label, prediction)[2,1]

  # Sensitivity(precision) and Specificity(recall)
  Sensitivity <- TP/(TP+FN)
  Specificity  <- TN/(TN+FP)
  
  # F1-score
  F1_score <- (TP)/(TP + 0.5 * (FP + FN))
  
  # MCC
  MCC_numerator <- TP * TN - FP * FN
  MCC_denominator <- sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  MCC <- MCC_numerator/MCC_denominator
  
  # Return the results as a tibble.
  tibble::tibble(
    K = K,
    Threshold = Threshold,
    TP = TP,
    TN = TN,
    FP = FP,
    FN = FN,
    Sensitivity = Sensitivity,
    Specificity = Specificity,
    F1_score = F1_score,
    MCC = MCC
  )
}

# All the combinations of parameters that are possible.
param_tibble <- expand.grid(K = 5:20, Threshold = seq(0.5,0.95,0.05))

# Runs each of these combinations.
out <- lapply(
  1L:nrow(param_tibble),
  function(i){
    result <- knn_fxn(param_tibble$K[i], param_tibble$Threshold[i])
    return(result)
  }
)

# lapply will return the list. bind_rows stacks the rows(list) on top of each other.
out <- dplyr::bind_rows(out)

# Export it to excel.
writexl::write_xlsx(out, "Results.xlsx")

# Graph
ggplot(out, aes(x=Specificity, y=Sensitivity, color = Threshold)) + 
  geom_point() +
  facet_wrap(~ K) + 
  theme(axis.text.x = element_text(angle = 90))

             