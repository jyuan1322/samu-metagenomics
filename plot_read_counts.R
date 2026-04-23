# Load libraries
library(jsonlite)
library(ggplot2)
library(dplyr)
library(tidyr)

# raw_counts_dir <- "/data/bwh-comppath-seq/jy1008/SaMu/data/metagenomics/raw/all_merged_fastqs"
raw_counts_dir <- "/data/local/jy1008/SaMu/scripts/samu-metagenomics/read_count_fastp_jsons"

# 1. List all JSON files in the directory
json_files <- list.files(raw_counts_dir, pattern = "\\.json$", full.names = TRUE)

# 2. Function to extract reads info from one fastp JSON
extract_reads <- function(json_file) {
  data <- fromJSON(json_file)
  
  # fastp JSON has summary under data$summary$before_filtering / after_filtering
  input_reads <- data$summary$before_filtering$total_reads
  passed_reads <- data$summary$after_filtering$total_reads
  
  # Get sample name from file name (strip path and .json)
  sample_name <- tools::file_path_sans_ext(basename(json_file))
  
  data.frame(
    Sample = sample_name,
    Input_Reads = input_reads,
    Passed_QC  = passed_reads
  )
}

# 3. Apply to all JSON files and combine
df <- do.call(rbind, lapply(json_files, extract_reads))
df$Sample <- gsub("_fastp", "", df$Sample)
write.csv(df, "fastp_read_counts.csv", row.names = FALSE)

# 4. Convert to long format for ggplot
df_long <- df %>%
  pivot_longer(cols = c(InputReads, PassedReads),
               names_to = "Category", values_to = "Reads")

# 5. Plot barplot
p <- ggplot(df_long, aes(x = Sample, y = Reads, fill = Category)) +
  geom_bar(stat = "identity", position = "identity") +
  theme_minimal() +
  labs(x = "Sample", y = "Number of reads", fill = "") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("fastp_reads_barplot.pdf", plot = p, width = 20, height = 6)
