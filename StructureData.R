# Load functions.
source("fxns.R")


####### Main Setup############################################################################

# Set parameters.
data_folder <- "files" # Folder that contains the subfolders "p10", "p11", etc.
file_map_path <- "file_map.csv"
batch_size <- 2e4 # When creating the file map, how many to process at a time.

# Create the file map.
#
# This only needs to be run once. It takes a long time to run, but will save you
# loads of time later. What this command does is create a file (file_map_path)
# that has 4 columns: patient_id, report_id, path, report. The report column
# contains the text from all reports. After you have run this command once, you
# will never need to run it again. You can comment it out by putting a # symbol
# before it.
create_FileMap(data_folder, file_map_path, batch_size)
