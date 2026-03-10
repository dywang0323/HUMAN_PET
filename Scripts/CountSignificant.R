# ============================================================
# GO Category Proportion of Significant MaAsLin3 Results
# Horizontal stacked bars (Abundance vs Prevalence)
# Largest category plotted on the LEFT
# ============================================================

.libPaths("/n/home07/dwangsky/R/x86_64-pc-linux-gnu-library/4.4")

suppressPackageStartupMessages({
library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(GO.db)
library(scales)
})

# ------------------------------------------------------------
# INPUT FILES
# ------------------------------------------------------------

DOG_HUMAN_FILE <- "/n/holystore01/LABS/huttenhower_lab/Users/Dongyu/PET_Project/Data/HumANn4_profile/GO/Informative_GO_NEW_1015/Unstraitified_HUMANn4_FUGAsseM/MaAsLin3/BP/dog_vs_human/all_results.tsv"

CAT_HUMAN_FILE <- "/n/holystore01/LABS/huttenhower_lab/Users/Dongyu/PET_Project/Data/HumANn4_profile/GO/Informative_GO_NEW_1015/Unstraitified_HUMANn4_FUGAsseM/MaAsLin3/BP/cat_vs_human/all_results.tsv"

DOG_CAT_FILE <- "/n/holystore01/LABS/huttenhower_lab/Users/Dongyu/PET_Project/Data/HumANn4_profile/GO/Informative_GO_NEW_1015/Unstraitified_HUMANn4_FUGAsseM/MaAsLin3/BP/dog_vs_cat/all_results.tsv"

OUTPUT_FILE <- "/n/holystore01/LABS/huttenhower_lab/Users/Dongyu/PET_Project/Data/HumANn4_profile/GO/go_category_proportion_panel.pdf"

# ------------------------------------------------------------
# ROOT GO TERMS
# ------------------------------------------------------------

root_terms <- c(
"GO:0008152",
"GO:0009987",
"GO:0006810",
"GO:0007165",
"GO:0006950"
)

root_labels <- c(
"Metabolism",
"Cellular processes",
"Transport",
"Signal transduction",
"Stress response"
)

# ------------------------------------------------------------
# BUILD GO CATEGORY MAP
# ------------------------------------------------------------

category_map <- data.frame()

for(i in seq_along(root_terms)){

root <- root_terms[i]

terms <- c(root, unlist(GOBPOFFSPRING[[root]]))

tmp <- data.frame(
GO_ID = terms,
category = root_labels[i],
stringsAsFactors = FALSE
)

category_map <- rbind(category_map,tmp)

}

category_map <- category_map %>%
distinct(GO_ID,.keep_all=TRUE)

# ------------------------------------------------------------
# LOAD MAASLIN RESULTS
# ------------------------------------------------------------

load_results <- function(file,label){

df <- read_tsv(file,show_col_types=FALSE)

df <- df %>%
filter(metadata=="species")

df$GO_ID <- str_extract(df$feature,"GO:\\d+")

df$significant <- df$qval_individual < 0.05

df$comparison <- label

df

}

dog_human <- load_results(DOG_HUMAN_FILE,"Human vs Dog")
cat_human <- load_results(CAT_HUMAN_FILE,"Human vs Cat")
dog_cat   <- load_results(DOG_CAT_FILE,"Dog vs Cat")

combined <- bind_rows(dog_human,cat_human,dog_cat)

# ------------------------------------------------------------
# REMOVE TRANSLATION MACHINERY
# ------------------------------------------------------------

combined <- combined %>%
filter(!str_detect(feature,
"translation|ribosome|tRNA|aminoacylation"))

# ------------------------------------------------------------
# MAP GO CATEGORY
# ------------------------------------------------------------

combined <- combined %>%
left_join(category_map,by="GO_ID") %>%
filter(!is.na(category))

combined <- combined %>%
filter(significant) %>%
distinct(comparison,model,GO_ID,category)

# ------------------------------------------------------------
# COUNT GO TERMS
# ------------------------------------------------------------

plot_data <- combined %>%
count(model,comparison,category)

# ------------------------------------------------------------
# CALCULATE PROPORTION + STACK ORDER
# ------------------------------------------------------------

plot_data <- plot_data %>%
group_by(model,comparison) %>%
mutate(
proportion = n / sum(n),
stack_order = rank(-proportion, ties.method="first")
) %>%
ungroup()

# ------------------------------------------------------------
# FACTOR ORDER
# ------------------------------------------------------------

plot_data$model <- factor(
plot_data$model,
levels=c("abundance","prevalence"),
labels=c("Abundance","Prevalence")
)

plot_data$comparison <- factor(
plot_data$comparison,
levels=c("Human vs Dog","Human vs Cat","Dog vs Cat")
)

# ------------------------------------------------------------
# COLORS
# ------------------------------------------------------------

category_colors <- c(
"Metabolism"="#4C78A8",
"Cellular processes"="#72B7B2",
"Transport"="#F58518",
"Signal transduction"="#E45756",
"Stress response"="#B279A2"
)

# ------------------------------------------------------------
# PLOT
# ------------------------------------------------------------

p <- ggplot(plot_data,
aes(
x = proportion,
y = comparison,
fill = category,
order = stack_order
)) +

geom_bar(
stat="identity",
width=0.65,
color="white",
size=0.2
) +

facet_grid(. ~ model) +

scale_fill_manual(values=category_colors) +

scale_x_continuous(
labels=scales::percent_format(),
limits=c(0,1),
breaks=c(0,0.25,0.5,0.75,1),
expand=c(0,0)
) +

scale_y_discrete(limits=rev) +

theme_bw() +

theme(
axis.text.y=element_text(size=11),
axis.text.x=element_text(size=10),
axis.title=element_text(size=11),
strip.text=element_text(size=11,face="bold"),
legend.position="right",
panel.grid.major.y=element_blank(),
panel.spacing=unit(0.2,"cm")
) +

labs(
y="",
x="Proportion of significant GO terms",
fill="GO functional category"
)

# ------------------------------------------------------------
# SAVE FIGURE
# ------------------------------------------------------------

ggsave(
OUTPUT_FILE,
p,
width=9,
height=3
)

message("Figure saved:")
message(OUTPUT_FILE)
