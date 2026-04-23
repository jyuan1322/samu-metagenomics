library(dplyr)
library(readr)
library(purrr)
library(tidyr)
library(ggplot2)

parent_dir <- "/data/local/jy1008/SaMu/proteomics"

files <- list.files(file.path(parent_dir, "data"),
                    pattern = "\\.csv$",
                    recursive = TRUE,
                    full.names = TRUE)

combined_df <- map_dfr(files, function(file) {
  
  #  read CSV with all columns as character and remove unnamed columns
  df <- read_csv(file, col_types = cols(.default = "c"), quote = "\"") %>%
    select(-starts_with("..."))

  # detect intensity columns dynamically
  raw_col  <- grep("^[0-9]+_.*$|^Sample", colnames(df), value = TRUE)
  iBAQ_col <- grep("^iBAQ", colnames(df), value = TRUE)
  print(paste("Detected raw count column:", raw_col))
  print(paste("Detected iBAQ column:", iBAQ_col))

  # rename columns (all_of handles dynamic names)
  df <- df %>% rename(
      RawCount       = all_of(raw_col),
      iBAQ_NormCount = all_of(iBAQ_col)
    )

  head(df)

  # Extract sample name from folder
  sample_name <- basename(file)
  print(sample_name)
  sample_name <- gsub("SaMu_sPROT1_rd1_pr1.0_rs0_", "", sample_name)
  sample_name <- gsub("_v1.csv", "", sample_name)
  print(sample_name)

  df %>%
    mutate(SampleID = sample_name)
})
combined_df$Global.Q.Value <- as.numeric(combined_df$Global.Q.Value)
combined_df$RawCount <- as.numeric(combined_df$RawCount)
combined_df$peptideCountsUnique <- as.numeric(combined_df$peptideCountsUnique)
combined_df$peptideCountsAll <- as.numeric(combined_df$peptideCountsAll)
combined_df$iBAQ_NormCount <- as.numeric(combined_df$iBAQ_NormCount)
combined_df$nbTrypticPeptides <- as.numeric(combined_df$nbTrypticPeptides)
combined_df$Relative.abundance <- as.numeric(sub("%", "", combined_df$Relative.abundance))
combined_df$contaminant <- as.logical(combined_df$contaminant)


# combined_wide <- combined_df %>%
#   select(Protein.Group, SampleID, RawCount) %>%
#   tidyr::pivot_wider(names_from = SampleID,
#                      values_from = RawCount)

# ----------------------------------------
# QC and visualization
# ----------------------------------------

# Make numeric Q-values
combined_df <- combined_df %>%
  mutate(
    Qval = as.numeric(Global.Q.Value),
    Category = ifelse(contaminant, "Contaminant", "Protein")
  )

# Define Q-value bins for non-contaminants
# bins up to 0.01, then one bin for everything above
bins <- c(seq(0, 0.01, by = 0.001), Inf)  # last bin = >0.01

combined_df <- combined_df %>%
  mutate(
    Qbin = case_when(
      Category == "Protein" ~ cut(
        Qval,
        breaks = bins,
        include.lowest = TRUE,
        labels = c(paste0(seq(0, 0.009, by = 0.001), "-", seq(0.001, 0.01, by = 0.001)), ">0.01")
      ),
      TRUE ~ "Contaminant"
    )
  )

# Count rows per category/bin
plot_data <- combined_df %>%
  group_by(Qbin) %>%
  summarise(count = n()) %>%
  ungroup()

p <- ggplot(plot_data, aes(x = Qbin, y = count, fill = Qbin)) +
  geom_col(show.legend = FALSE) +
  scale_fill_manual(values = c("Contaminant" = "red", rep("steelblue", length(unique(plot_data$Qbin)) - 1))) +
  labs(
    x = "Category / Q-value bin",
    y = "Number of rows",
    title = "Filtering Proteomics Data: Contaminants and Q-value bins"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(parent_dir, "results", "proteomics_qvalue_bins.pdf"), p, width = 8, height = 6)


# Perform the actual QC

# filter out contaminants and keep only those with Global.Q.Value < 0.01
clean_df <- combined_df %>%
  filter(!contaminant) %>% 
  filter(as.numeric(Global.Q.Value) < 0.01)
# > dim(clean_df)
# [1] 361934     17
# > dim(combined_df)
# [1] 365848     17

combined_wide <- clean_df %>%
  select(Protein.Group, Gene.names, SampleID, RawCount) %>%  # use RawCount
  pivot_wider(
    names_from = SampleID,
    values_from = RawCount
  )

# combined_wide should have Protein.Group and Protein.Names
combined_wide <- combined_wide %>%
  dplyr::rename(ID = Protein.Group,
                name = Gene.names)  # or Gene.names if you prefer

# now make unique
# combined_unique <- make_unique(
#   combined_wide,
#   "ID",
#   "name"
# )
combined_wide$ID <- paste(combined_wide$ID, combined_wide$name, sep = ";")

ecols <- grep("SaMu", colnames(combined_wide))

# ----------------------------------------
# Metadata import
# ----------------------------------------
meta_df <- read.csv("/data/local/jy1008/SaMu/metadata/SaMu_sarcopeniestatus_majorcovariates_v13_16012026.csv",
                header = TRUE, stringsAsFactors = FALSE)
meta_df$label <- paste0("SaMu", meta_df$record_id)

# Filter metadata for important covariates
meta_df <- meta_df[, c(
  "record_id",
  "label",
  "Block",
  "Full.SaMu",
  "sarc_status",
  "EWGSOP_strength",
  "EWGSOP_mass",
  "EWGSOP_performance",
  "age_def",
  "sex",
  "smke",
  "alco",
  "nutr_score",
  "bmi",
  "stvol"
)]
meta_df <- meta_df %>% filter(Full.SaMu == 1)
meta_df$smke <- ifelse(meta_df$smke %in% c("1", "2", "3"), 1,
                       ifelse(meta_df$smke == "0", 0, NA))
meta_df$alco <- ifelse(meta_df$alco %in% c("1", "2"), 0,
                       ifelse(meta_df$alco %in% c("3", "4"), 1, NA))
print("Converting sarc_status to numeric, converting non-numeric values to NA")
meta_df$sarc_status <- as.numeric(meta_df$sarc_status)
# binarize sarc_status
meta_df$sarc_status_bin <- ifelse(meta_df$sarc_status > 0, "Sarc", "NoSarc")
meta_df$sarc_status_bin <- factor(meta_df$sarc_status_bin, levels = c("NoSarc", "Sarc"))


# intensity_matrix: columns are sample IDs
# coldata: metadata table

# Keep only samples that exist in intensity_matrix
coldata <- meta_df %>%
  filter(label %in% colnames(combined_wide)[ecols])

common <- intersect(colnames(combined_wide)[ecols], coldata$label)

coldata <- coldata[match(common, coldata$label), ]



# # Optional: reorder rows to match the column order in intensity_matrix
# coldata <- coldata[match(colnames(combined_wide)[ecols], coldata$label), ]

# # check that sample names match intensity matrix
all(coldata$label %in% colnames(combined_wide)[ecols])

# import metadata
# coldata <- data.frame(
#   label = colnames(intensity_matrix),
#   condition = c("Control","Control","Treatment","Treatment"), # replace with your design
#   replicate = c(1,2,1,2)
# )


# ----------------------------------------
# Running DEP2
# ----------------------------------------
library(DEP2)

# Create SummarizedExperiment
coldata$condition <- factor(coldata$sarc_status_bin, levels = c("NoSarc", "Sarc"))
coldata$replicate <- 1

combined_wide <- combined_wide %>% as.data.frame()
rownames(combined_wide) <- combined_wide$ID
combined_wide$name <- combined_wide$ID
coldata <- coldata[, c("label", "condition", "replicate")]
coldata$condition <- as.character(coldata$condition)
# Replace NA with a placeholder like "Unknown"
coldata$condition[is.na(coldata$condition)] <- "Unknown"
coldata$replicate <- as.numeric(gsub("SaMu", "", coldata$label))
se <- make_se(
  combined_wide,
  columns = ecols,
  expdesign = coldata
)

# drop samples with unknown diagnosis
se <- se[, colData(se)$condition != "Unknown"]
# se <- se[rowData(se)$gene != "NA", ]

# Test different filtering thresholds
fractions <- seq(0.1, 0.9, by = 0.05)
thresholds <- c(10, 25, 50, 100)

param_grid <- expand.grid(
  fraction = fractions,
  thr = thresholds
)

results <- param_grid %>%
  rowwise() %>%
  mutate(
    n_rows = {
      se_tmp <- filter_se(
        se,
        thr = thr,
        fraction = fraction
      )
      nrow(se_tmp)
    }
  ) %>%
  ungroup()

p <- ggplot(results, aes(x = fraction, y = n_rows, color = factor(thr))) +
  geom_line(size = 1) +
  labs(
    x = "thr (allowed missing replicates)",
    y = "missing fraction threshold",
    color = "missing replicate threshold"
  ) +
  theme_minimal()
ggsave(file.path(parent_dir, "results", "dep2_se_filtering_thresholds.pdf"), p, width = 8, height = 6)

# Filter low-quality proteins
# se_filt <- filter_se(se,                       
#                      thr = 50,  ## the threshold of missing number in at least one condition
#                      fraction = 0.5 ## the threshold of missing occupancy in each protein
#                     )
# NOTE: the documentation for the fraction seems to be incorrect.
# It is not the fraction of missing values, but rather the fraction of valid values required.
se_filt <- filter_se(se,                       
                     thr = 50,  ## the threshold of missing number in at least one condition
                     fraction = 0.15 ## the threshold of missing occupancy in each protein
                    )
dim(se)
dim(se_filt)

pdf(file.path(parent_dir, "results", "dep2_filtering_effect.pdf"), width = 20, height = 6)
plot_frequency(se) + ggtitle("Identification overlap before filter")
plot_frequency(se_filt) + ggtitle("Identification overlap after before filter")
dev.off()





# Normalize
se_norm <- normalize_vsn(se_filt)
pdf(file.path(parent_dir, "results", "dep2_normalization_effect.pdf"), width = 6, height = 20)
plot_normalization(se_filt, se_norm)
dev.off()



# Plot mean-variance relationship: confirm there is no relationship (horizontal fit)
library(vsn)
pdf(file.path(parent_dir, "results", "meanSdPlot.pdf"), width = 8, height = 6)
meanSdPlot(assay(se_norm))
dev.off()

# Impute missing values
pdf(file.path(parent_dir, "results", "dep2_missing_values_heatmap.pdf"), width = 8, height = 6)
plot_missval(se_filt)
dev.off()

# se_imp <- DEP2::impute(se_norm) # fails
se_imp <- DEP2::impute(se_norm, fun = "MinDet")   # for minimal detection limit
# se_imp <- DEP2::impute(se_norm, fun = "MinProb", q = 0.01)

# Write the imputed, VSN-normalized data to a CSV
# Extract imputed VSN matrix
imputed_matrix <- assay(se_imp)

# transpose so rows = samples, columns = proteins
proteomics_matrix <- t(imputed_matrix)

proteomics_df <- as.data.frame(proteomics_matrix)
proteomics_df$Sample <- rownames(proteomics_df)

# move sample column first
proteomics_df <- proteomics_df[, c("Sample", setdiff(colnames(proteomics_df), "Sample"))]

write.csv(
  proteomics_df,
  file.path(parent_dir, "results", "proteomics_dep2_vsn_imputed_matrix_p15filt.csv"),
  row.names = FALSE
)

metadata_df <- as.data.frame(colData(se_imp))
write.csv(
  metadata_df,
  file.path(parent_dir, "results", "proteomics_dep2_vsn_imputed_metadata_p15filt.csv"),
  row.names = FALSE
)




# -----



# Test differential expression
colData(se_imp)$condition <- factor(colData(se_imp)$condition)
colData(se_imp)$condition <- relevel(
  colData(se_imp)$condition,
  ref = "NoSarc"
)

diff <- test_diff(
  se_imp,
  type = "manual",
  test = "Sarc_vs_NoSarc",
  fdr.type = "BH"
)

# Add significance annotations
results <- add_rejections(diff)
write.csv(as.data.frame(rowData(results)),
          file.path(parent_dir, "results", "dep2_default_results_p15filt.csv"))

# Quick volcano plot
pdf(file.path(parent_dir, "results", "dep2_default_volcano_plot_p15filt.pdf"), width = 8, height = 6)
plot_volcano(results, adjusted = F,
             add_threshold_line = "intersect", pCutoff  = 0.05, fcCutoff = 1)
dev.off()

# ----------------------------------------
# Heatmap of the data
# ----------------------------------------
library(ComplexHeatmap)
library(circlize)

# Use normalized intensities
exprs_mat <- assay(se_filt)

log_mat <- log10(exprs_mat + 1e-6)
scaled_mat <- t(scale(t(log_mat)))

# # filter low-variance proteins (for cleaner heatmap)
# # Compute row-wise variance (ignoring NAs)
# num_top_prots <- 100
# row_vars <- apply(exprs_mat, 1, var, na.rm = TRUE)
# # Remove any proteins where variance is NA
# row_vars <- row_vars[!is.na(row_vars)]
# # Get names of top 50 most variable proteins
# top50_proteins <- names(sort(row_vars, decreasing = TRUE))[1:num_top_prots]

# grab by smallest p-value from the DE results
top50_proteins <- rowData(results) %>%
  as.data.frame() %>%
  arrange(Sarc_vs_NoSarc_p.val) %>%      # sort by this column
  slice_head(n = 50) %>%
  pull(ID)           # extract this column

# 1. Get the top 50 proteins by smallest p-value
# results_df <- as.data.frame(rowData(results))
# top50_proteins <- results_df %>%
#   arrange(Sarc_vs_NoSarc_p.val) %>%      # or "qval" if that’s the column name
#   slice_head(n = 50) %>%
#   pull(name)              # or "ID" depending on your column

# Subset the expression matrix
# exprs_top50 <- exprs_mat[top50_proteins, ]
scaled_mat_top50 <- scaled_mat[top50_proteins, ]

# Get sample annotation
sample_anno <- data.frame(
  SarcStatus = colData(se_norm)$condition
)
rownames(sample_anno) <- colnames(scaled_mat_top50)

# Compute breakpoints
sorted_vals <- sort(unique(as.numeric(scaled_mat_top50)), na.last = NA)
min_val <- sorted_vals[1]              # absolute minimum (imputed)
# second_min <- sorted_vals[2]           # second smallest (observed)
max_val <- max(scaled_mat_top50, na.rm = TRUE)  # maximum

# Define color function
col_fun <- colorRamp2(
  breaks = c(min_val, max_val),
  colors = c("black", "red")
)

# Add sarc/no_sarc annotation
sample_anno <- data.frame(SarcStatus = colData(se_imp)$condition)
rownames(sample_anno) <- colnames(scaled_mat_top50)
ha <- HeatmapAnnotation(
  SarcStatus = sample_anno$SarcStatus,
  col = list(SarcStatus = c("NoSarc" = "steelblue", "Sarc" = "orange", "Unknown" = "gray"))
)


# pdf(file.path( parent_dir, "results", "proteomics_massspec_logscaled_heatmap.pdf"), width = 10, height = 12)
# Heatmap(
#   scaled_mat_top50,
#   col = col_fun,
#   show_row_names = TRUE,                # show protein names
#   row_names_side = "left",              # labels on left side
#   row_names_gp = gpar(fontsize = 10),   # adjust font size
#   name = "Intensity",
#   top_annotation = ha,
#   show_column_names = TRUE,
#   clustering_distance_rows = "euclidean",
#   clustering_distance_columns = "euclidean",
#   clustering_method_rows = "complete",
#   clustering_method_columns = "complete",
#   column_title = "Samples",
#   row_title = "Proteins"
# )
# dev.off()

col_order <- order(colnames(scaled_mat_top50))
pdf(file.path( parent_dir, "results", "proteomics_heatmap_p15filt.pdf"), width = 15, height = 12)
Heatmap(
  scaled_mat_top50,
  col = col_fun,
  na_col = "gray80",
  show_row_names = TRUE,                # show protein names
  row_names_side = "left",              # labels on left side
  row_names_gp = gpar(fontsize = 10),   # adjust font size
  name = "scaled Log10(Intensity)",
  top_annotation = ha,
  show_column_names = TRUE,
  cluster_columns = FALSE,
  column_order = col_order,
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "complete",
  column_title = "Samples",
  row_title = "Proteins"
)
dev.off()
