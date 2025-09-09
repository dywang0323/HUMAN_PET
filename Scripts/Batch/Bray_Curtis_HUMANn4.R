# braycurtis_pcoa_permanova.R

# --------- Libraries ---------
library(ggplot2)
library(dplyr)
library(readr)
library(tibble)
library(vegan)
library(ape)

# --------- Load Data ---------
go_data <- read_tsv("/n/holystore01/LABS/huttenhower_lab/Users/Dongyu/PET_Project/Data/HumANn4_profile/GO/RE_DO_0817/DepthFilter/Cat_Dog_Human_GO_rename_main_lab_RM_pre_DEPTH_4_0909.tsv")
meta <- read_csv("/n/holystore01/LABS/huttenhower_lab/Users/Dongyu/PET_Project/Data/HumANn4_profile/S15_metadata_consist.csv")

# --------- Prepare GO Matrix ---------
# Transpose GO data: samples as rows, GO terms as columns
go_matrix <- as.data.frame(t(go_data[,-1]))
colnames(go_matrix) <- go_data$GO
rownames(go_matrix) <- colnames(go_data)[-1]

# Convert sample IDs to characters
meta$sample_id_metaphlan <- as.character(meta$sample_id_metaphlan)
rownames(go_matrix) <- as.character(rownames(go_matrix))

# Filter metadata to match GO matrix
meta_matched <- meta %>%
  filter(sample_id_metaphlan %in% rownames(go_matrix)) %>%
  arrange(match(sample_id_metaphlan, rownames(go_matrix)))

# Reorder GO matrix to match metadata
go_matrix <- go_matrix[meta_matched$sample_id_metaphlan, ]

# Sanity check
stopifnot(all(rownames(go_matrix) == meta_matched$sample_id_metaphlan))

# --------- Normalize by row sums (relative abundances) ---------
go_matrix <- sweep(go_matrix, 1, rowSums(go_matrix), FUN = "/")

# --------- Replace Zeros (Optional for Aitchison) ---------
min_relab <- min(go_matrix[go_matrix > 0]) / 10
go_matrix_pseudo <- go_matrix
go_matrix_pseudo[go_matrix_pseudo == 0] <- min_relab

# --------- Bray-Curtis Dissimilarity ---------
dist_bray <- vegdist(go_matrix, method = "bray")

# --------- Set Weights by Species (dynamic) ---------
species_counts <- table(meta_matched$species)
weights <- 1 / species_counts[meta_matched$species]

# --------- Weighted PCoA ---------
cmd_res <- wcmdscale(dist_bray, k = 2, eig = TRUE, w = weights)

# Percent variance explained (use only positive eigenvalues)
eig_vals <- cmd_res$eig[cmd_res$eig > 0]
percent_explained <- 100 * eig_vals / sum(eig_vals)
PC1 <- round(percent_explained[1], 2)
PC2 <- round(percent_explained[2], 2)

# First two PCs as dataframe
pcoa_2PC <- tibble(
  PC1 = cmd_res$points[,1],
  PC2 = cmd_res$points[,2],
  sample_id_metaphlan = rownames(cmd_res$points)
)

# Merge with metadata
pcoa_df <- left_join(pcoa_2PC, meta_matched, by = "sample_id_metaphlan")

# --------- PERMANOVA (adonis2) ---------
adonis_result <- adonis2(dist_bray ~ species, data = meta_matched, permutations = 999)
print(adonis_result)

# Save PERMANOVA results to .tsv
adonis_df <- as.data.frame(adonis_result)
write_tsv(adonis_df, "/n/holystore01/LABS/huttenhower_lab/Users/Dongyu/PET_Project/Data/HumANn4_profile/GO/RE_DO_0817/DepthFilter/PCoA_GO_braycurtis_PERMANOVA.tsv")

# --------- PCoA Plot ---------
species.colors <- c(
  cat = "indianred",
  dog = "#E69F00",
  human = "#56B4E0",
  `HMP1-II` = "#56B4E0", 
  Madagascar = "#0f3b50"
)

# Choose plot style: "classic" (no top/right) or "boxed" (all borders)
plot_style <- "boxed"

p <- ggplot(pcoa_df, aes(x = PC1, y = PC2, fill = species)) +
  geom_point(shape = 21, size = 3, alpha = 0.8, stroke = 0.2) +
  scale_fill_manual(values = species.colors) +
  labs(
    x = paste0("PCoA1 [", PC1, "%]"),
    y = paste0("PCoA2 [", PC2, "%]"),
    fill = "Host Species"
  )

if (plot_style == "classic") {
  p <- p + theme_classic() +
    theme(
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12)
    )
} else if (plot_style == "boxed") {
  p <- p + theme_classic() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12)
    )
}

# --------- Save to TIFF (300 dpi, adjustable size) ---------
ggsave(
  filename = "/n/holystore01/LABS/huttenhower_lab/Users/Dongyu/PET_Project/Data/HumANn4_profile/GO/RE_DO_0817/DepthFilter/PCoA_GO_braycurtis.tif",
  plot = p,
  width = 6,   # adjust width here (inches)
  height = 4,  # adjust height here (inches)
  dpi = 300,
  units = "in",
  device = "tiff"
)
