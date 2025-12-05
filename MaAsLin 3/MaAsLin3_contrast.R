# -------------------------
# Load libraries
# -------------------------
library(dplyr)
library(maaslin3)

# -------------------------
# Load metadata
# -------------------------
metadata_raw <- read.csv(
  "/n/holystore01/LABS/huttenhower_lab/Users/Dongyu/PET_Project/Data/HumANn4_profile/S15_metadata_consist.csv"
)

# -------------------------
# Load GO abundance data (Biological Process)
# -------------------------
go_raw <- read.delim(
  "/n/holystore01/LABS/huttenhower_lab/Users/Dongyu/PET_Project/Data/HumANn4_profile/GO/Informative_GO_NEW_1015/Unstraitified_HUMANn4_FUGAsseM/HUMAnN4_FUGAsseM_main_RM_Preva_Norm_unify_1204.tsv",
  row.names = 1,
  check.names = FALSE
)

# -------------------------
# Transpose to samples x features
# -------------------------
go_transposed <- as.data.frame(t(go_raw), check.names = FALSE)

# -------------------------
# Match metadata to abundance
# -------------------------
metadata_matched <- metadata_raw %>%
  filter(sample_id_metaphlan %in% rownames(go_transposed)) %>%
  arrange(match(sample_id_metaphlan, rownames(go_transposed)))

# Set rownames to match
rownames(metadata_matched) <- metadata_matched$sample_id_metaphlan

# Sanity checks
stopifnot(identical(rownames(metadata_matched), rownames(go_transposed)))
stopifnot("species" %in% colnames(metadata_matched))
stopifnot("study_readable" %in% colnames(metadata_matched))

# -------------------------
# Define pairwise MaAsLin3 function
# -------------------------
run_pairwise_maaslin_with_contrast <- function(species1, species2, output_dir) {
  message("\n Running MaAsLin3 for: ", species1, " vs ", species2)
  
  # Subset metadata and define factor levels
  subset_meta <- metadata_matched %>%
    filter(species %in% c(species1, species2)) %>%
    mutate(species = factor(species, levels = c(species1, species2)))
  
  if (nlevels(droplevels(subset_meta$species)) < 2) {
    warning(" Skipping comparison: ", species1, " vs ", species2,
            " — only one group present. Sample count: ", nrow(subset_meta))
    return(NULL)
  }
  
  # Subset GO abundance
  subset_abun <- go_transposed[rownames(subset_meta), ]
  
  # Create output directory
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Model formula: fixed effect (species), random effect (study)
  formula_str <- "~ species + (1 | study_readable)"
  
  # Run MaAsLin3
  fit_out <- maaslin3(
    input_data = subset_abun,
    input_metadata = subset_meta,
    output = output_dir,
    formula = formula_str,
    normalization = "TSS",      # Normalize each sample to relative abundance
    transform = "LOG",          # Apply log transformation
    standardize = FALSE,
    save_models = TRUE,
    augment = TRUE,
    max_significance = 0.1,
    cores = 2
  )
  
  # -------------------------
  # Contrast matrix setup
  # -------------------------
  coef_names <- rownames(fit_out$fit_data_abundance$coefficients)
  contrast_vals <- setNames(rep(0, length(coef_names)), coef_names)
  
  # Identify the coefficient for the target species
  contrast_term <- paste0("species", species2)
  
  if (!(contrast_term %in% coef_names)) {
    warning("⚠ Contrast coefficient not found in model: ", contrast_term)
    message("Available coefficients: ", paste(coef_names, collapse = ", "))
    return(NULL)
  }
  
  # Assign contrast weight
  contrast_vals[contrast_term] <- 1
  
  # Build contrast matrix
  contrast_matrix <- matrix(contrast_vals, nrow = 1)
  colnames(contrast_matrix) <- names(contrast_vals)
  rownames(contrast_matrix) <- "species_contrast"
  
  # Save contrast matrix
  write.table(contrast_matrix,
              file = file.path(output_dir, "contrast_matrix.txt"),
              sep = "\t", quote = FALSE, col.names = NA)
  
  # Run contrast test
  contrast_out <- maaslin_contrast_test(fit_out, contrast_matrix)
  contrast_results <- contrast_out$fit_data_abundance$results
  
  # Save full results
  write.csv(contrast_results,
            file = file.path(output_dir, "contrast_results.csv"),
            row.names = FALSE)
  
  # Filter by q-value if available
  if ("qval" %in% colnames(contrast_results)) {
    contrast_sig <- contrast_results %>% filter(qval < 0.1)
  } else {
    warning("⚠ No 'qval' column found in contrast results — skipping filtering.")
    contrast_sig <- contrast_results[0, ]
  }
  
  # Save significant results only
  write.csv(contrast_sig,
            file = file.path(output_dir, "contrast_results_qval_lt_0.1.csv"),
            row.names = FALSE)
  
  # Report count
  message("  ↳ Significant features (q < 0.1): ", nrow(contrast_sig))
}

# -------------------------
# Run pairwise comparisons
# -------------------------
output_base <- "/n/holystore01/LABS/huttenhower_lab/Users/Dongyu/PET_Project/Data/HumANn4_profile/GO/Informative_GO_NEW_1015/Unstraitified_HUMANn4_FUGAsseM/MaAsLin3/GO_Overall"

# List all host species (e.g., human, cat, dog)
species_levels <- sort(unique(na.omit(metadata_matched$species)))

# Generate all pairwise combinations
species_pairs <- combn(species_levels, 2, simplify = FALSE)

# Ensure 'human' is always the reference (first)
species_pairs <- lapply(species_pairs, function(p) {
  if ("human" %in% p) c("human", setdiff(p, "human")) else p
})

# Loop through each pair and run MaAsLin3
for (pair in species_pairs) {
  out_dir <- file.path(output_base, paste0(pair[2], "_vs_", pair[1], "_human_ref"))
  run_pairwise_maaslin_with_contrast(pair[1], pair[2], out_dir)
}
