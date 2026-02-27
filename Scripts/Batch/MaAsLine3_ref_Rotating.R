# module load gcc/12.2.0-fasrc01
# module load R/4.4.3-fasrc01
# ============================================================
# Run MaAsLin3 with Rotating Reference Levels
# human / dog / cat as reference
# ============================================================

library(dplyr)
library(maaslin3)

# Load metadata

metadata_raw <- read.csv(
  "/n/holystore01/LABS/huttenhower_lab/Users/Dongyu/PET_Project/Data/HumANn4_profile/S15_metadata_consist.csv"
)

# Load GO abundance

go_raw <- read.delim(
  "/n/holystore01/LABS/huttenhower_lab/Users/Dongyu/PET_Project/Data/HumANn4_profile/GO/Informative_GO_NEW_1015/HUMANn4_Informative_GO/Cat_Dog_Human_genefamily_main_InforGO_Preva_1016.tsv",
  row.names = 1,
  check.names = FALSE
)

go_transposed <- as.data.frame(t(go_raw), check.names = FALSE)

# Match metadata

metadata_matched <- metadata_raw %>%
  filter(sample_id_metaphlan %in% rownames(go_transposed)) %>%
  arrange(match(sample_id_metaphlan, rownames(go_transposed)))

rownames(metadata_matched) <- metadata_matched$sample_id_metaphlan

stopifnot(identical(rownames(metadata_matched), rownames(go_transposed)))


# Function to run MaAsLin with chosen reference


run_with_reference <- function(ref_species) {

  message("\n==============================")
  message("Running with reference: ", ref_species)
  message("==============================")

  metadata_ref <- metadata_matched
  
  metadata_ref$species <- factor(
    metadata_ref$species,
    levels = c(ref_species,
               setdiff(unique(metadata_ref$species), ref_species))
  )

  output_dir <- paste0(
    "/n/holystore01/LABS/huttenhower_lab/Users/Dongyu/PET_Project/Data/HumANn4_profile/GO/Informative_GO_NEW_1015/HUMANn4_Informative_GO/MaAsLin3_RotatingRef/",
    ref_species, "_reference"
  )

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  maaslin3(
    input_data = go_transposed,
    input_metadata = metadata_ref,
    output = output_dir,
    formula = "~ species + (1 | study_readable)",
    normalization = "TSS",
    transform = "LOG",
    standardize = FALSE,
    save_models = TRUE,
    augment = TRUE,
    max_significance = 0.1,
    cores = 2
  )
}


# Run for all three references

run_with_reference("human")
run_with_reference("dog")
run_with_reference("cat")
