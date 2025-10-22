# --------- Libraries ---------
library(ggplot2)
library(dplyr)
library(readr)
library(tibble)
library(vegan)
library(ape)

# --------- File Paths ---------
in_path  <- "/n/holystore01/LABS/huttenhower_lab/Users/Dongyu/PET_Project/Data/HumANn4_profile/GO/Informative_GO_NEW_1015/Unstraitified_HUMANn4_FUGAsseM/HUMAnN4_FUGAsseM_main_RM_Preva_Norm_1022.tsv"
meta_path <- "/n/holystore01/LABS/huttenhower_lab/Users/Dongyu/PET_Project/Data/HumANn4_profile/S15_metadata_consist.csv"

# --------- Step 1: Read & consolidate by GO ---------
go <- read_tsv(in_path, show_col_types = FALSE)

if (!"GO" %in% names(go)) {
  stop("Column `GO` not found in the input file. Check the header.")
}

# Updated: Exclude metadata columns like `source`
sample_cols <- setdiff(names(go), c("GO", "source"))

# Consolidate by GO (sum across duplicates)
go_summed <- go %>%
  group_by(GO) %>%
  summarise(across(all_of(sample_cols), ~ sum(.x, na.rm = TRUE)), .groups = "drop") %>%
  select(GO, all_of(sample_cols))

# --------- Step 2: Prepare GO Matrix ---------
meta <- read_csv(meta_path)

# Transpose: samples as rows, functions as columns
go_matrix <- as.data.frame(t(go_summed[,-1]))
colnames(go_matrix) <- go_summed$GO
rownames(go_matrix) <- colnames(go_summed)[-1]

# Align IDs with metadata
meta$sample_id_metaphlan <- as.character(meta$sample_id_metaphlan)
rownames(go_matrix) <- as.character(rownames(go_matrix))

meta_matched <- meta %>%
  filter(sample_id_metaphlan %in% rownames(go_matrix)) %>%
  arrange(match(sample_id_metaphlan, rownames(go_matrix)))

go_matrix <- go_matrix[meta_matched$sample_id_metaphlan, ]
stopifnot(all(rownames(go_matrix) == meta_matched$sample_id_metaphlan))

# --------- Normalize by row sums ---------
go_matrix <- sweep(go_matrix, 1, rowSums(go_matrix), FUN = "/")
go_matrix[is.na(go_matrix)] <- 0

# Drop zero-abundance samples
zero_sum_samples <- rownames(go_matrix)[rowSums(go_matrix) == 0]
if (length(zero_sum_samples) > 0) {
  message("Dropping zero-abundance samples: ", paste(zero_sum_samples, collapse = ", "))
}
go_matrix <- go_matrix[rowSums(go_matrix) > 0, ]
meta_matched <- meta_matched[meta_matched$sample_id_metaphlan %in% rownames(go_matrix), ]

# --------- Bray-Curtis Dissimilarity ---------
dist_bray <- vegdist(go_matrix, method = "bray")

# --------- Set Weights by Species ---------
species_counts <- table(meta_matched$species)
weights <- 1 / species_counts[meta_matched$species]

# --------- Weighted PCoA ---------
cmd_res <- wcmdscale(dist_bray, k = 2, eig = TRUE, w = weights)

eig_vals <- cmd_res$eig[cmd_res$eig > 0]
percent_explained <- 100 * eig_vals / sum(eig_vals)
PC1 <- round(percent_explained[1], 2)
PC2 <- round(percent_explained[2], 2)

pcoa_2PC <- tibble(
  PC1 = cmd_res$points[,1],
  PC2 = cmd_res$points[,2],
  sample_id_metaphlan = rownames(cmd_res$points)
)

pcoa_df <- left_join(pcoa_2PC, meta_matched, by = "sample_id_metaphlan")

# --------- Ensure drawing order ---------
pcoa_df$species <- factor(pcoa_df$species, levels = c("dog", "cat", "human"))
pcoa_df <- pcoa_df %>% arrange(species)

# --------- PERMANOVA ---------
set.seed(1234)  # For reproducibility
adonis_result <- adonis2(dist_bray ~ species, data = meta_matched, permutations = 999)
print(adonis_result)

adonis_df <- as.data.frame(adonis_result)
write_tsv(adonis_df, "/n/holystore01/LABS/huttenhower_lab/Users/Dongyu/PET_Project/Data/HumANn4_profile/GO/Informative_GO_NEW_1015/Unstraitified_HUMANn4_FUGAsseM/Bray_Curtis/PCoA_GO_braycurtis_PERMANOVA.tsv")

# --------- PCoA Plot ---------
species.colors <- c(
  cat = "indianred",
  dog = "#E69F00",
  human = "#56B4E0"
)

p <- ggplot(pcoa_df, aes(x = PC1, y = PC2, fill = species)) +
  geom_point(shape = 21, size = 3, alpha = 0.8, stroke = 0.2) +
  scale_fill_manual(values = species.colors) +
  labs(
    x = paste0("PCoA1 [", PC1, "%]"),
    y = paste0("PCoA2 [", PC2, "%]"),
    fill = "Host Species"
  ) +
  theme_bw() +
  theme(
    axis.line = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12)
  )

ggsave(
  filename = "/n/holystore01/LABS/huttenhower_lab/Users/Dongyu/PET_Project/Data/HumANn4_profile/GO/Informative_GO_NEW_1015/Unstraitified_HUMANn4_FUGAsseM/Bray_Curtis/PCoA_GO_braycurtis.tif",
  plot = p,
  width = 6,
  height = 4,
  dpi = 300,
  units = "in",
  device = "tiff"
)
