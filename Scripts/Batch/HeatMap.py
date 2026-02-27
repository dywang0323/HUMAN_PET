import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist
from matplotlib.colors import ListedColormap
from matplotlib.gridspec import GridSpec


# INPUT FILES

BP_TSV = "/n/holystore01/LABS/huttenhower_lab/Users/Dongyu/PET_Project/Data/HumANn4_profile/GO/Informative_GO_NEW_1015/HUMANn4_Informative_GO/Split/Max_Min_MF/MF_top10_abundant.tsv"
META_CSV = "/n/holystore01/LABS/huttenhower_lab/Users/Dongyu/PET_Project/Data/HumANn4_profile/S15_metadata_consist.csv"

OUTPUT_FILE = "/n/holystore01/LABS/huttenhower_lab/Users/Dongyu/PET_Project/Data/HumANn4_profile/GO/Informative_GO_NEW_1015/HUMANn4_Informative_GO/Split/Max_Min_MF/MF_top10_abundant_HeatMap.tiff"


# PARAMETERS

CELL_HEIGHT = 0.2
BASE_WIDTH = 6
DPI = 300

HEATMAP_CMAP = "cividis"
LINKAGE_METHOD = "average"
DISTANCE_METRIC = "correlation"
LOG_PSEUDOCOUNT = 1e-6

SPECIES_ORDER = ["dog", "cat", "human"]

SPECIES_COLORS = {
    "dog": "#E69F00",
    "cat": "indianred",
    "human": "#56B4E0"
}


# LOAD DATA

bp = pd.read_csv(BP_TSV, sep="\t")
data_raw = bp.set_index(bp.columns[0])
meta = pd.read_csv(META_CSV)

common_samples = [
    c for c in data_raw.columns
    if c in set(meta["sample_id_metaphlan"])
]

data_raw = data_raw[common_samples]
data_log = np.log10(data_raw + LOG_PSEUDOCOUNT)

meta_sub = meta.set_index("sample_id_metaphlan").reindex(data_log.columns)

if meta_sub["species"].isna().any():
    raise ValueError("Missing species metadata.")

# CLEAN FOR CORRELATION DISTANCE


# Remove non-finite values
data_log = data_log.replace([np.inf, -np.inf], np.nan)
data_log = data_log.dropna(axis=0, how="any")
data_log = data_log.dropna(axis=1, how="any")

# Remove zero-variance rows and columns (critical for correlation distance)
data_log = data_log.loc[data_log.std(axis=1) > 0]
data_log = data_log.loc[:, data_log.std(axis=0) > 0]

# Re-align metadata after filtering columns
meta_sub = meta_sub.loc[data_log.columns]

if data_log.shape[0] < 2 or data_log.shape[1] < 2:
    raise ValueError("Not enough data left after filtering for clustering.")


# ADAPTIVE COLOR SCALING

VMIN = np.percentile(data_log.values, 5)
VMAX = np.percentile(data_log.values, 95)

print(f"Adaptive VMIN: {VMIN:.3f}")
print(f"Adaptive VMAX: {VMAX:.3f}")

# CLUSTERING


row_linkage = linkage(
    pdist(data_log, metric=DISTANCE_METRIC),
    method=LINKAGE_METHOD
)

col_linkage = linkage(
    pdist(data_log.T, metric=DISTANCE_METRIC),
    method=LINKAGE_METHOD
)

# Get consistent leaf order
row_d = dendrogram(row_linkage, no_plot=True)
col_d = dendrogram(col_linkage, no_plot=True)

row_order = row_d["leaves"]
col_order = col_d["leaves"]

data_log_ord = data_log.iloc[row_order, col_order]
meta_ord = meta_sub.iloc[col_order]

# Clip only for visualization
data_log_ord = data_log_ord.clip(lower=VMIN, upper=VMAX)


# PREPARE LABELS


go_labels = (
    data_log_ord.index.astype(str)
    .str.replace(r"\[BP\]", "", regex=True)
    .str.strip()
)

n_rows = len(go_labels)

heatmap_height = n_rows * CELL_HEIGHT
TOP_MARGIN = 1.2
FIG_HEIGHT = heatmap_height + TOP_MARGIN

max_label_len = max(len(g) for g in go_labels)
extra_width = max_label_len * 0.06
FIG_WIDTH = BASE_WIDTH + extra_width


# FIGURE


fig = plt.figure(figsize=(FIG_WIDTH, FIG_HEIGHT), constrained_layout=True)

gs = GridSpec(
    2, 3,
    figure=fig,
    width_ratios=[1, 4, 0.4],
    height_ratios=[1, heatmap_height]
)

# --- Column dendrogram ---
ax_col = fig.add_subplot(gs[0, 1])
dendrogram(
    col_linkage,
    ax=ax_col,
    no_labels=True,
    color_threshold=0,
    above_threshold_color="black"
)
ax_col.axis("off")

# --- Row dendrogram ---
ax_row = fig.add_subplot(gs[1, 0])
dendrogram(
    row_linkage,
    ax=ax_row,
    orientation="left",
    no_labels=True,
    color_threshold=0,
    above_threshold_color="black"
)
ax_row.axis("off")

# --- Heatmap ---
ax_hm = fig.add_subplot(gs[1, 1])

im = ax_hm.imshow(
    data_log_ord,
    aspect="auto",
    cmap=HEATMAP_CMAP,
    vmin=VMIN,
    vmax=VMAX,
    interpolation="nearest"
)

ax_hm.set_xticks([])
ax_hm.set_yticks(np.arange(n_rows))
ax_hm.set_yticklabels(go_labels, fontsize=9)

ax_hm.yaxis.tick_right()
ax_hm.tick_params(axis="y", length=0, pad=6)

for label in ax_hm.get_yticklabels():
    label.set_clip_on(False)

for spine in ax_hm.spines.values():
    spine.set_visible(False)

# --- Species bar ---
species = meta_ord["species"].astype("category")
species = species.cat.set_categories(SPECIES_ORDER, ordered=True)

species_codes = np.array([
    SPECIES_ORDER.index(s) if s in SPECIES_ORDER else -1
    for s in species
]).reshape(1, -1)

species_cmap = ListedColormap([SPECIES_COLORS[s] for s in SPECIES_ORDER])

ax_species = ax_hm.inset_axes([0, 1.02, 1, 0.05])
ax_species.imshow(species_codes, aspect="auto", cmap=species_cmap)
ax_species.set_xticks([])
ax_species.set_yticks([])

# --- Colorbar ---
ax_cb = fig.add_subplot(gs[1, 2])
plt.colorbar(im, cax=ax_cb)
ax_cb.set_ylabel("log10(abundance)")

# SAVE

plt.savefig(
    OUTPUT_FILE,
    dpi=DPI,
    format="tiff",
    bbox_inches="tight",
    pil_kwargs={"compression": "tiff_lzw"}
)

plt.close()

print(f"Saved clustered heatmap with adaptive scaling to {OUTPUT_FILE}")
