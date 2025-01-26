import pandas as pd
import numpy as np
import os
from scipy.spatial.distance import pdist, squareform
from sklearn.manifold import MDS
from skbio.stats.ordination import pcoa
from skbio.stats.distance import permanova, DistanceMatrix
import matplotlib.pyplot as plt
import seaborn as sns

# --- File paths ---
abundance_file = "/ourdisk/hpc/nullspace/dywang/dont_archive/PET/Data/ECs_consist.tsv"
metadata_file = "/ourdisk/hpc/nullspace/dywang/dont_archive/PET/Data/S15_metadata_consist.csv"
output_dir = "/ourdisk/hpc/nullspace/dywang/dont_archive/PET/Results/"

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

# --- Load and preprocess abundance data ---
print("Loading abundance data...")
abundance_df = pd.read_csv(abundance_file, sep="\t")

# Transpose data to make sample IDs the first column
abundance_df = abundance_df.set_index('EC_number').transpose().reset_index()
abundance_df.rename(columns={'index': 'sample_id_metaphlan'}, inplace=True)

# Standardize sample ID formatting (remove whitespace, enforce lowercase)
abundance_df['sample_id_metaphlan'] = abundance_df['sample_id_metaphlan'].astype(str).str.strip().str.lower()

# --- Load metadata ---
print("Loading metadata...")
metadata_df = pd.read_csv(metadata_file, usecols=['sample_id_metaphlan', 'species'])
metadata_df['sample_id_metaphlan'] = metadata_df['sample_id_metaphlan'].astype(str).str.strip().str.lower()
metadata_df = metadata_df.drop_duplicates(subset=['sample_id_metaphlan'])

# Merge metadata with abundance data
abundance_df = abundance_df.merge(metadata_df, on='sample_id_metaphlan', how='inner')

# --- Prepare data for analysis ---
abundance_matrix = abundance_df.drop(columns=['sample_id_metaphlan', 'species']).apply(pd.to_numeric, errors='coerce').fillna(0)
abundance_df.set_index('sample_id_metaphlan', inplace=True)

# --- Compute Bray-Curtis dissimilarity matrix ---
print("Computing Bray-Curtis dissimilarity matrix...")
bray_curtis_dist_matrix = pdist(abundance_matrix, metric='braycurtis')
bray_curtis_dist = pd.DataFrame(
    squareform(bray_curtis_dist_matrix),
    index=abundance_df.index,
    columns=abundance_df.index
)

# Save Bray-Curtis dissimilarity matrix
bray_curtis_output_path = os.path.join(output_dir, "bray_curtis_dissimilarity_fixed.csv")
bray_curtis_dist.to_csv(bray_curtis_output_path)
print(f"Bray-Curtis dissimilarity matrix saved to '{bray_curtis_output_path}'")

# --- Align sample IDs before PERMANOVA ---
common_ids = sorted(set(bray_curtis_dist.index) & set(abundance_df.index))
print(f"Number of common sample IDs found: {len(common_ids)}")

if not common_ids:
    raise ValueError("No common sample IDs found between Bray-Curtis matrix and metadata. Please check ID formatting.")

# Subset matrices to include only common samples
bray_curtis_dist = bray_curtis_dist.loc[common_ids, common_ids]
abundance_df = abundance_df.loc[common_ids]

# --- PERMANOVA analysis ---
print("Performing PERMANOVA test...")
abundance_df = abundance_df.dropna(subset=['species'])

distance_matrix = DistanceMatrix(bray_curtis_dist.values, ids=bray_curtis_dist.index)
permanova_results = permanova(distance_matrix, abundance_df['species'], permutations=999)

# Calculate R² (explained variance)
r_squared = permanova_results['test statistic'] / (permanova_results['test statistic'] + permanova_results['number of permutations'])

# Save PERMANOVA results
permanova_output_path = os.path.join(output_dir, "permanova_results.txt")
with open(permanova_output_path, "w") as f:
    f.write(f"PERMANOVA Results:\n{permanova_results}\n")
    f.write(f"Pseudo-F: {permanova_results['test statistic']}\n")
    f.write(f"p-value: {permanova_results['p-value']}\n")
    f.write(f"R-squared: {r_squared:.4f}\n")
    f.write(f"Number of permutations: {permanova_results['number of permutations']}\n")

print(f"PERMANOVA results saved to '{permanova_output_path}'")

# --- PCoA Analysis ---
print("Performing PCoA analysis...")
pcoa_results = pcoa(bray_curtis_dist)

# Extract variance explained by axes
pcoa_var_explained = pcoa_results.proportion_explained * 100  # Convert to percentage
pcoa_df = pd.DataFrame(pcoa_results.samples[['PC1', 'PC2']])
pcoa_df['species'] = abundance_df['species'].values

# Plot PCoA with variance percentages
def plot_ordination_with_variance(df, x, y, title, color_col, var_x, var_y, output_path):
    plt.figure(figsize=(8, 8))
    sns.scatterplot(data=df, x=x, y=y, hue=color_col, palette="Dark2", s=100, edgecolor='black')
    plt.title(title, fontsize=14, weight='bold')
    plt.xlabel(f"{x} ({var_x:.2f}%)", fontsize=12, weight='bold')
    plt.ylabel(f"{y} ({var_y:.2f}%)", fontsize=12, weight='bold')
    plt.legend(title=color_col)
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.savefig(output_path, format='tiff', dpi=300, bbox_inches='tight')
    plt.show()
    print(f"PCoA plot saved to {output_path}")

pcoa_output_path = os.path.join(output_dir, "PCoA_beta_diversity.tif")
plot_ordination_with_variance(pcoa_df, 'PC1', 'PC2', "PCoA of Functional Beta Diversity", 'species', 
                              pcoa_var_explained[0], pcoa_var_explained[1], pcoa_output_path)

# --- NMDS Analysis ---
print("Performing NMDS analysis...")
nmds = MDS(n_components=2, dissimilarity='precomputed', random_state=42, metric=False, max_iter=3000, eps=1e-9)
nmds_results = nmds.fit_transform(bray_curtis_dist)
stress_value = nmds.stress_  # Extract stress value

abundance_df['NMDS1'], abundance_df['NMDS2'] = nmds_results[:, 0], nmds_results[:, 1]

# Plot NMDS with stress value
def plot_nmds_with_stress(df, x, y, title, color_col, stress, output_path):
    plt.figure(figsize=(8, 8))
    sns.scatterplot(data=df, x=x, y=y, hue=color_col, palette="Dark2", s=100, edgecolor='black')
    plt.title(f"{title}\nStress: {stress:.4f}", fontsize=14, weight='bold')
    plt.xlabel(x, fontsize=12, weight='bold')
    plt.ylabel(y, fontsize=12, weight='bold')
    plt.legend(title=color_col)
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.savefig(output_path, format='tiff', dpi=300, bbox_inches='tight')
    plt.show()
    print(f"NMDS plot saved to {output_path}")

nmds_output_path = os.path.join(output_dir, "NMDS_beta_diversity.tif")
plot_nmds_with_stress(abundance_df, 'NMDS1', 'NMDS2', "NMDS of Functional Beta Diversity", 'species', 
                      stress_value, nmds_output_path)

print("Analysis completed successfully.")
