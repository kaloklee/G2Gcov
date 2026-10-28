# data-raw/prepare_kb_data.R
# Script to generate and save package datasets
# Run manually via devtools::load_all() or devtools::document()

# ---- Setup ----
library(readr)      # for read_csv
library(dplyr)      # for cleaning
library(usethis)    # for use_data()

# ---- Read raw data ----
kb_raw <- read_csv("data-raw/kb_covars.csv")

# ---- Clean / Transform ----
kb_data <- kb_raw %>%
  rename_with(tolower)

# Store a plain data frame, without readr's tibble class or its "spec" and
# "problems" attributes
kb_data <- data.frame(lapply(kb_data, as.vector))

# ---- Print Columns and Data Types ----
message("--- Column Names (kb_data) ---")
print(names(kb_data))
message("\n--- Data Types (kb_data) ---")
print(str(kb_data))
message("---------------------------------\n")

# ---- Verify ----
stopifnot(!any(is.na(kb_data)))

# ---- Small sample for quick examples: all rows of the first 22 subjects ----
kb_covars_sample <- kb_data[kb_data$id %in% 1:22, ]

# ---- Save ----
usethis::use_data(kb_data, overwrite = TRUE)
usethis::use_data(kb_covars_sample, overwrite = TRUE)
