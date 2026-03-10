# ============================================================
# Host-variable metabolic functional repertoire
# Prevalence bubble panel (publication-ready)
# ============================================================

.libPaths("/n/home07/dwangsky/R/x86_64-pc-linux-gnu-library/4.4")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(GO.db)
  library(tidyr)
})

first <- dplyr::first

# ------------------------------------------------------------
# INPUT FILES
# ------------------------------------------------------------

METADATA_FILE <- "/n/holystore01/LABS/huttenhower_lab/Users/Dongyu/PET_Project/Data/HumANn4_profile/S15_metadata_consist.csv"

GO_BP_FILE <- "/n/holystore01/LABS/huttenhower_lab/Users/Dongyu/PET_Project/Data/HumANn4_profile/GO/Informative_GO_NEW_1015/Unstraitified_HUMANn4_FUGAsseM/Split/GO_BP.tsv"

OUTPUT_FILE <- "/n/holystore01/LABS/huttenhower_lab/Users/Dongyu/PET_Project/Data/HumANn4_profile/GO/metabolic_prevalence_panel.pdf"

# ------------------------------------------------------------
# COLORS
# ------------------------------------------------------------

host_colors <- c(
  dog   = "#E69F00",
  cat   = "indianred",
  human = "#56B4E0"
)

# ------------------------------------------------------------
# METABOLIC ROOT TERMS
# ------------------------------------------------------------

metabolism_roots <- c(
  "GO:0005975","GO:0006520","GO:0006629","GO:0009117",
  "GO:0006091","GO:0051186","GO:0006807","GO:0006790",
  "GO:0019748","GO:0006082","GO:0044281"
)

metabolism_labels <- c(
  "Carbohydrate metabolism",
  "Amino acid metabolism",
  "Lipid metabolism",
  "Nucleotide metabolism",
  "Energy metabolism",
  "Cofactor metabolism",
  "Nitrogen metabolism",
  "Sulfur metabolism",
  "Secondary metabolism",
  "Organic acid metabolism",
  "Small molecule metabolism"
)

# ------------------------------------------------------------
# BUILD CATEGORY MAP
# ------------------------------------------------------------

category_map <- data.frame()

for(i in seq_along(metabolism_roots)){
  
  root <- metabolism_roots[i]
  
  terms <- c(root, unlist(GOBPOFFSPRING[[root]]))
  
  tmp <- data.frame(
    GO_ID = terms,
    category = metabolism_labels[i],
    stringsAsFactors = FALSE
  )
  
  category_map <- rbind(category_map,tmp)
  
}

category_map <- category_map %>%
  arrange(category) %>%
  distinct(GO_ID,.keep_all=TRUE)

# ------------------------------------------------------------
# LOAD DATA
# ------------------------------------------------------------

metadata <- read_csv(METADATA_FILE,show_col_types=FALSE)

go <- read_tsv(GO_BP_FILE,show_col_types=FALSE)
go <- as.data.frame(go)

colnames(go)[1] <- "feature"

# ------------------------------------------------------------
# MATCH SAMPLE IDS
# ------------------------------------------------------------

dog_cols <- intersect(colnames(go),
                      metadata$sample_id_metaphlan[metadata$species=="dog"])

cat_cols <- intersect(colnames(go),
                      metadata$sample_id_metaphlan[metadata$species=="cat"])

human_cols <- intersect(colnames(go),
                        metadata$sample_id_metaphlan[metadata$species=="human"])

# ------------------------------------------------------------
# EXTRACT GO IDS
# ------------------------------------------------------------

go$GO_ID <- str_extract(go$feature,"GO:\\d+")

go <- merge(go,category_map,by="GO_ID",all.x=TRUE)

go <- go[!is.na(go$category),]

go <- go[!str_detect(go$feature,
                     "translation|ribosome|tRNA|aminoacylation"),]

# ------------------------------------------------------------
# CLEAN LABELS
# ------------------------------------------------------------

go$feature <- gsub("GO:[0-9]+: \\[BP\\] ","",go$feature)
go$feature <- gsub("biosynthetic process","biosynthesis",go$feature)
go$feature <- gsub("catabolic process","catabolism",go$feature)
go$feature <- gsub("metabolic process","metabolism",go$feature)

go$feature <- trimws(go$feature)

# ------------------------------------------------------------
# PREVALENCE
# ------------------------------------------------------------

threshold <- 1e-6

go$dog   <- rowMeans(go[,dog_cols] > threshold)
go$cat   <- rowMeans(go[,cat_cols] > threshold)
go$human <- rowMeans(go[,human_cols] > threshold)

prev_table <- go[,c("feature","category","dog","cat","human")]

prev_table <- prev_table %>%
  group_by(feature,category) %>%
  summarise(
    dog   = mean(dog),
    cat   = mean(cat),
    human = mean(human),
    .groups="drop"
  )

# ------------------------------------------------------------
# HOST VARIABILITY
# ------------------------------------------------------------

prev_table <- prev_table %>%
  rowwise() %>%
  mutate(host_variation = sd(c(dog,cat,human))) %>%
  ungroup()

prev_table <- prev_table %>%
  group_by(category) %>%
  arrange(desc(host_variation),.by_group=TRUE) %>%
  slice_head(n=4) %>%
  ungroup()

# ------------------------------------------------------------
# LONG FORMAT
# ------------------------------------------------------------

plot_df <- pivot_longer(
  prev_table,
  cols=c("dog","cat","human"),
  names_to="Host",
  values_to="Prevalence"
)

plot_df$Host <- factor(plot_df$Host,
                       levels=c("dog","cat","human"))

# ------------------------------------------------------------
# ORDER PATHWAYS
# ------------------------------------------------------------

order_table <- prev_table %>%
  arrange(category,desc(host_variation))

ordered_features <- unique(order_table$feature)

plot_df$feature <- factor(plot_df$feature,
                          levels=rev(ordered_features))

plot_df$category <- factor(plot_df$category,
                           levels=unique(order_table$category))

# ------------------------------------------------------------
# BUBBLE PANEL
# ------------------------------------------------------------

p <- ggplot(plot_df,
            aes(x=Host,
                y=feature,
                size=Prevalence,
                color=Host)) +
  
  geom_point(alpha=0.7) +
  
  facet_grid(category~.,
             scales="free_y",
             space="free_y") +
  
  scale_color_manual(values=host_colors) +
  
  scale_size_continuous(
    range=c(5,14),
    limits=c(0,1),
    name="Prevalence"
  ) +
  
  scale_x_discrete(expand=c(0.35,0.35)) +
  
  theme_bw() +
  
  theme(
    
    strip.background=element_rect(fill="grey97",color="black"),
    
    strip.text.y=element_text(
      angle=0,
      size=12,
      face="bold"
    ),
    
    axis.text.y=element_text(size=10),
    axis.text.x=element_text(size=12,face="bold"),
    
    axis.title.y=element_text(size=12),
    axis.title.x=element_text(size=12),
    
    legend.title=element_text(size=11),
    legend.text=element_text(size=10),
    
    panel.grid.major.x=element_line(color="grey85"),
    panel.grid.major.y=element_blank(),
    panel.grid.minor=element_blank(),
    
    legend.position="right"
  ) +
  
  labs(
    x="Host",
    y="Metabolic pathway",
    title="Host-variable metabolic functional repertoire"
  )

# ------------------------------------------------------------
# SAVE
# ------------------------------------------------------------

ggsave(
  OUTPUT_FILE,
  p,
  width=9.5,
  height=16
)

message("Figure saved:")
message(OUTPUT_FILE)
