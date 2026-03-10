# ============================================================
# Forest Plot of Metabolic GO Terms from MaAsLin3
# Unique GO terms + priority-based category assignment
# ============================================================

.libPaths("/n/home07/dwangsky/R/x86_64-pc-linux-gnu-library/4.4")

suppressPackageStartupMessages({
library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(GO.db)
})

# ------------------------------------------------------------
# INPUT FILES
# ------------------------------------------------------------

DOG_HUMAN_FILE <- "/n/holystore01/LABS/huttenhower_lab/Users/Dongyu/PET_Project/Data/HumANn4_profile/GO/Informative_GO_NEW_1015/Unstraitified_HUMANn4_FUGAsseM/MaAsLin3/BP/dog_vs_human/all_results.tsv"

CAT_HUMAN_FILE <- "/n/holystore01/LABS/huttenhower_lab/Users/Dongyu/PET_Project/Data/HumANn4_profile/GO/Informative_GO_NEW_1015/Unstraitified_HUMANn4_FUGAsseM/MaAsLin3/BP/cat_vs_human/all_results.tsv"

DOG_CAT_FILE <- "/n/holystore01/LABS/huttenhower_lab/Users/Dongyu/PET_Project/Data/HumANn4_profile/GO/Informative_GO_NEW_1015/Unstraitified_HUMANn4_FUGAsseM/MaAsLin3/BP/dog_vs_cat/all_results.tsv"

OUTPUT_FILE <- "/n/holystore01/LABS/huttenhower_lab/Users/Dongyu/PET_Project/Data/HumANn4_profile/GO/metabolism_grouped_forest_plot_top20.pdf"

# ------------------------------------------------------------
# METABOLISM ROOT TERMS
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
"Secondary metabolite metabolism",
"Organic acid metabolism",
"Small molecule metabolism"
)

# ------------------------------------------------------------
# BUILD GO CATEGORY MAP
# ------------------------------------------------------------

category_map <- data.frame()

for(i in seq_along(metabolism_roots)){

root <- metabolism_roots[i]

terms <- c(root, unlist(GOBPOFFSPRING[[root]]))

category_map <- rbind(
category_map,
data.frame(
GO_ID = terms,
category = metabolism_labels[i]
)
)

}

category_map <- unique(category_map)

# ------------------------------------------------------------
# PRIORITY ORDER (NO MERGING)
# ------------------------------------------------------------

category_priority <- c(
"Carbohydrate metabolism",
"Amino acid metabolism",
"Lipid metabolism",
"Energy metabolism",
"Cofactor metabolism",
"Secondary metabolite metabolism",
"Nucleotide metabolism",
"Nitrogen metabolism",
"Sulfur metabolism",
"Organic acid metabolism",
"Small molecule metabolism"
)

category_map$category <- factor(
category_map$category,
levels = category_priority
)

# ------------------------------------------------------------
# REMOVE DAG DUPLICATES
# ------------------------------------------------------------

category_map <- category_map %>%
arrange(category) %>%
distinct(GO_ID, .keep_all = TRUE)

# ------------------------------------------------------------
# FUNCTION TO LOAD RESULTS
# ------------------------------------------------------------

process_results <- function(file,label){

df <- read_tsv(file,show_col_types=FALSE)

df <- df %>%
filter(metadata=="species") %>%
filter(model=="abundance")

df$GO_ID <- str_extract(df$feature,"GO:\\d+")

df <- df %>%
left_join(category_map,by="GO_ID")

df <- df %>% filter(!is.na(category))

df$GO_term <- str_replace(
df$feature,
"GO:[0-9]+: \\[BP\\] ",
""
)

df$comparison <- label

df$lower <- df$coef - 1.96*df$stderr
df$upper <- df$coef + 1.96*df$stderr

df

}

# ------------------------------------------------------------
# LOAD RESULTS
# ------------------------------------------------------------

dog_human <- process_results(DOG_HUMAN_FILE,"Human vs Dog")
cat_human <- process_results(CAT_HUMAN_FILE,"Human vs Cat")
dog_cat   <- process_results(DOG_CAT_FILE,"Dog vs Cat")

combined <- bind_rows(dog_human,cat_human,dog_cat)

# ------------------------------------------------------------
# REMOVE TRANSLATION PROCESSES
# ------------------------------------------------------------

combined <- combined %>%
filter(!str_detect(GO_term,"translation|ribosome|tRNA|aminoacylation"))

# ------------------------------------------------------------
# ENSURE UNIQUE GO TERMS PER COMPARISON
# ------------------------------------------------------------

combined <- combined %>%
group_by(comparison,GO_term,category) %>%
summarise(
coef=mean(coef),
stderr=mean(stderr),
lower=mean(lower),
upper=mean(upper),
.groups="drop"
)

# ------------------------------------------------------------
# SELECT TOP HUMAN PATHWAYS
# ------------------------------------------------------------

human_only <- combined %>%
filter(comparison %in% c("Human vs Dog","Human vs Cat"))

top_human <- human_only %>%
group_by(comparison) %>%
arrange(desc(abs(coef)),.by_group=TRUE) %>%
slice_head(n=20) %>%
ungroup()

selected_pathways <- unique(top_human$GO_term)

top_terms <- combined %>%
filter(GO_term %in% selected_pathways)

# ------------------------------------------------------------
# PANEL ORDER
# ------------------------------------------------------------

top_terms$comparison <- factor(
top_terms$comparison,
levels=c("Human vs Dog","Human vs Cat","Dog vs Cat")
)

# ------------------------------------------------------------
# HOST ENRICHMENT
# ------------------------------------------------------------

top_terms$host <- NA

top_terms$host[top_terms$comparison=="Human vs Dog" & top_terms$coef>0] <- "dog"
top_terms$host[top_terms$comparison=="Human vs Dog" & top_terms$coef<0] <- "human"

top_terms$host[top_terms$comparison=="Human vs Cat" & top_terms$coef>0] <- "cat"
top_terms$host[top_terms$comparison=="Human vs Cat" & top_terms$coef<0] <- "human"

top_terms$host[top_terms$comparison=="Dog vs Cat" & top_terms$coef>0] <- "cat"
top_terms$host[top_terms$comparison=="Dog vs Cat" & top_terms$coef<0] <- "dog"

# ------------------------------------------------------------
# COLORS
# ------------------------------------------------------------

host_colors <- c(
human="#56B4E0",
dog="#E69F00",
cat="indianred"
)

# ------------------------------------------------------------
# ORDER PATHWAYS
# ------------------------------------------------------------

order_table <- top_terms %>%
group_by(category,GO_term) %>%
summarise(score=max(abs(coef)),.groups="drop") %>%
arrange(category,desc(score))

ordered_terms <- unique(order_table$GO_term)

top_terms$GO_term <- factor(
top_terms$GO_term,
levels=rev(ordered_terms)
)

# ------------------------------------------------------------
# FOREST PLOT
# ------------------------------------------------------------

p <- ggplot(top_terms,
aes(x=coef,
y=GO_term,
xmin=lower,
xmax=upper,
color=host)) +

geom_vline(xintercept=0,linetype="dashed") +

geom_errorbar(height=0.2) +

geom_point(size=4) +

facet_grid(category~comparison,
scales="free_y",
space="free_y") +

scale_color_manual(values=host_colors) +

theme_bw() +

theme(
strip.text.y=element_text(angle=0,size=11,face="bold"),
strip.text.x=element_text(size=12,face="bold"),
axis.text.y=element_text(size=9),
panel.grid.major.x = element_line(color="grey85"),
panel.grid.minor = element_blank()
) +

labs(
x="Effect size (β coefficient)",
y="Metabolic pathway",
color="Enriched in",
title="Top Metabolic Functional Differences Between Hosts"
)

# ------------------------------------------------------------
# SAVE FIGURE
# ------------------------------------------------------------

ggsave(
OUTPUT_FILE,
p,
width=15,
height=11
)

message("Forest plot saved:")
message(OUTPUT_FILE)
