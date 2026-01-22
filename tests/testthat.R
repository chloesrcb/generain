
library(testthat)
source("helper_load_libraries.R")
library(sf)
# Get all files in the folder "R"
functions_folder <- "./R"
files <- list.files(functions_folder, full.names = TRUE)
# # load all functions in files
invisible(lapply(files, function(f) source(f, echo = FALSE)))
# library(generain)

# test_check("generain")
test_dir("./testthat")
