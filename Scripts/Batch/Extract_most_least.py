import pandas as pd
import numpy as np


# INPUT FILES

GO_TSV = "/n/holystore01/LABS/huttenhower_lab/Users/Dongyu/PET_Project/Data/HumANn4_profile/GO/Informative_GO_NEW_1015/HUMANn4_Informative_GO/Split/GO_MF_norm.tsv"
META_CSV = "/n/holystore01/LABS/huttenhower_lab/Users/Dongyu/PET_Project/Data/HumANn4_profile/S15_metadata_consist.csv"


# OUTPUT FILES

OUT_TOP = "/n/holystore01/LABS/huttenhower_lab/Users/Dongyu/PET_Project/Data/HumANn4_profile/GO/Informative_GO_NEW_1015/HUMANn4_Informative_GO/Split/Max_Min_MF/MF_top10_abundant.tsv"
OUT_BOTTOM = "/n/holystore01/LABS/huttenhower_lab/Users/Dongyu/PET_Project/Data/HumANn4_profile/GO/Informative_GO_NEW_1015/HUMANn4_Informative_GO/Split/Max_Min_MF/MF_bottom10_abundant.tsv"


# LOAD DATA

go = pd.read_csv(GO_TSV, sep="\t")
meta = pd.read_csv(META_CSV)

GO_COL = go.columns[0]


# LONG FORMAT

go_long = go.melt(
    id_vars=GO_COL,
    var_name="sample_id_metaphlan",
    value_name="abundance"
)


# MERGE METADATA

merged = go_long.merge(
    meta[["sample_id_metaphlan", "species"]],
    on="sample_id_metaphlan",
    how="inner"
)

# Keep dog / cat / human
merged = merged[merged["species"].isin(["dog", "cat", "human"])]


# MEAN ABUNDANCE PER GO PER SPECIES

mean_abundance = (
    merged
    .groupby(["species", GO_COL], as_index=False)["abundance"]
    .mean()
)


# TOP 10 MOST ABUNDANT

top10 = (
    mean_abundance
    .sort_values(["species", "abundance"], ascending=[True, False])
    .groupby("species")
    .head(10)
)

top10_terms = top10[GO_COL].unique()


# TOP 10 LEAST ABUNDANT

bottom10 = (
    mean_abundance
    .sort_values(["species", "abundance"], ascending=[True, True])
    .groupby("species")
    .head(10)
)

bottom10_terms = bottom10[GO_COL].unique()


# FILTER ORIGINAL GO TABLE

top_go = go[go[GO_COL].isin(top10_terms)]
bottom_go = go[go[GO_COL].isin(bottom10_terms)]


# SAVE FILES

top_go.to_csv(OUT_TOP, sep="\t", index=False)
bottom_go.to_csv(OUT_BOTTOM, sep="\t", index=False)

print(f"Saved top 10 abundant GO table: {OUT_TOP}")
print(f"Saved bottom 10 abundant GO table: {OUT_BOTTOM}")
print(f"Top GO terms: {top_go.shape[0]}")
print(f"Bottom GO terms: {bottom_go.shape[0]}")
