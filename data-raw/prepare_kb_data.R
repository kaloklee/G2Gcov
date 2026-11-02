# data-raw/prepare_kb_data.R
# Prepares package datasets using only base R
# Run manually to generate data/kb_data.rda and data/kb_data_sample.rda

# ---- 1. Read raw data ----
kb_data <- read.csv("data-raw/kb_covars.csv", stringsAsFactors = FALSE)

# ---- 2. Clean and standardize column names ----
names(kb_data) <- tolower(names(kb_data))

# ---- 3. Add discrete-time event indicator ----
# 1 = event occurred at terminal interval, 0 = otherwise
kb_data$status <- as.integer(kb_data$week == kb_data$num_weeks & kb_data$censor == 0)

# ---- 4. Integrity checks ----
stopifnot(!any(is.na(kb_data)))
stopifnot(all(kb_data$status %in% c(0L, 1L)))

# ---- 5. Draw 20-subject sample for quick examples ----
set.seed(42)
sample_ids <- sort(sample(unique(kb_data$id), size = 20))
kb_data_sample <- kb_data[kb_data$id %in% sample_ids, ]
rownames(kb_data_sample) <- NULL

# ---- 6. Save datasets to data/ using base R ----
if (!dir.exists("data")) {
  dir.create("data")
}

save(kb_data, file = "data/kb_data.rda", compress = "xz")
save(kb_data_sample, file = "data/kb_data_sample.rda", compress = "xz")

message("Datasets successfully created in data/ with xz compression.")
