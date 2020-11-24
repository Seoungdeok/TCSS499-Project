###### Functions.

# Function for creating a map between radiology reports and patient ids.
# 
# data_folder: the folder that contains the subfolders "p10", "p11", etc.
# map_file: the filename of the file map that will be created.
# batch_size: the number of reports to process at a time.
create_FileMap <- function(data_folder, map_file, batch_size = 2e4){
  # Check if the file map has been created and is valid.
  valid <- FALSE
  if(file.exists(map_file)){
    if(nrow(read_FileMap(map_file)) == 227835L){
      valid <- TRUE
    }
  }
  
  # If the file map is not valid or is missing, then create it. If the file map
  # became corrupted somehow, then you should delete it and start over.
  if(!valid){
    # There are 227835 txt files distributed across about 65,000 folders. This
    # will take quite a while. Fortunately, the structure of the data folder
    # will never change so you will only need to do this once.
    #
    # The following command will create a vector where each element gives the
    # path to a report.
    files <- list.files(
      data_folder,
      recursive = TRUE
    )
    # Remove the data_folder component of the path.
    files <- gsub("files/|files\\\\", "", files)
    # Extract information.
    split_data <- strsplit(files, "/|\\\\")
    patient <- unlist(lapply(split_data, "[[", 2L))
    report <- gsub(".txt$", "", unlist(lapply(split_data, "[[", 3L)))
    
    ids <- (1L:length(files)) %/% batch_size
    num_batches <- tail(ids, 1L)
    for(i in 0L:num_batches){
      matched_ids <- which(ids == i)
      
      reports <- character(length(matched_ids))
      for(j in seq_along(reports)){
        # Note that replacing \n is necessary because \n means new line in a txt
        # file. This is reversed when the data is loaded with read_FileMap.
        reports[j] <- gsub(
          "\\n","\\\\n",
          read_report(data_folder, files[matched_ids][j])
        )
      }
      
      # Put the extracted information into a tabular data structure.
      # 
      # Note that a tibble is essentially a data.frame, except that it has some 
      # features that make it a little easier to work with.
      structured_data <- tibble::tibble(
        patient_id = patient[matched_ids],
        report_id = report[matched_ids],
        path = files[matched_ids],
        report = reports
      )
      
      # Save the data.
      if(!file.exists(map_file)){
        readr::write_csv(structured_data, map_file)
      }else{
        readr::write_csv(structured_data, map_file, append = TRUE)
      }
      print(paste0(i, " / ", num_batches, " batches."))
      flush.console()
    }
    
  }
}

# Function for reading the file map.
# 
# map_file: the filename of the file map that was created with create_FileMap.
read_FileMap <- function(map_file){
  # Load the data.
  out <- readr::read_csv(map_file, col_types = readr::cols(), trim_ws = FALSE)
  
  # Undo the replacement of \n.
  out$report <- gsub("\\\\n","\\\n",out$report)
  
  # Return the data.
  out
}

# Function for reading a report.
# 
# data_folder: the folder that contains the subfolders "p10", "p11", etc.
# path: the component of the path to a radiology report that is found in the
# path column of the file map.
read_report <- function(data_folder, path){
  full_path <- file.path(data_folder, path)
  readr::read_file(full_path)
}

