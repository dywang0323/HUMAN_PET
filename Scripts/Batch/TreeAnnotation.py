#!/usr/bin/env python3

import os
import bz2
import pandas as pd
from Bio import Phylo


# INPUT FILES

TREE_FILE = "/n/holystore01/LABS/huttenhower_lab/Users/Dongyu/PET_Project/Data/HumANn4_profile/GO/Informative_GO_NEW_1015/HUMANn4_Informative_GO/PhylogeneticTree/GO_SGB_subtree.nwk"

ABUNDANCE_FILE = "/n/holystore01/LABS/huttenhower_lab/Users/Dongyu/PET_Project/Data/HumANn4_profile/GO/Informative_GO_NEW_1015/HUMANn4_Informative_GO/PhylogeneticTree/GO_SGB_abundance.tsv"

METADATA_FILE = "/n/holystore01/LABS/huttenhower_lab/Users/Dongyu/PET_Project/Data/HumANn4_profile/S15_metadata_consist.csv"

SGB_TAXONOMY_FILE = "/n/holystore01/LABS/huttenhower_lab/Users/Dongyu/PET_Project/Data/MetaPhlAn4/mpa_vJun23_CHOCOPhlAnSGB_202403_species.txt.bz2"

OUTPUT_DIR = "/n/holystore01/LABS/huttenhower_lab/Users/Dongyu/PET_Project/Data/HumANn4_profile/GO/Informative_GO_NEW_1015/HUMANn4_Informative_GO/PhylogeneticTree/GO_Annotation/"

os.makedirs(OUTPUT_DIR, exist_ok=True)

# HOST COLORS

HOST_COLORS = {
    "dog": "#E69F00",
    "cat": "#CD5C5C",
    "human": "#56B4E0"
}


# LOAD TREE (NUMERIC SGB TIPS)


print("Loading numeric subtree...")
tree = Phylo.read(TREE_FILE, "newick")

tree_ids = [t.name.strip() for t in tree.get_terminals()]
tree_id_set = set(tree_ids)

print(f"Tree leaves: {len(tree_ids)}")


# LOAD ABUNDANCE MATRIX


print("Loading abundance matrix...")
abund = pd.read_csv(ABUNDANCE_FILE, sep="\t")

# Ensure correct index
abund = abund.set_index("SGB_ID")

# Convert SGB#### → digits only (to match tree tips)
abund.index = abund.index.str.replace("SGB", "", regex=False)

abund_id_set = set(abund.index)

common_ids = sorted(tree_id_set.intersection(abund_id_set))
abund = abund.loc[common_ids]

print(f"Matched abundance rows: {len(abund)}")


# LOAD METADATA


print("Loading metadata...")
meta = pd.read_csv(METADATA_FILE)
meta.columns = meta.columns.str.strip()
meta["species"] = meta["species"].str.lower().str.strip()

sample_to_host = dict(zip(meta["sample_id_metaphlan"], meta["species"]))

host_samples = {h: [] for h in HOST_COLORS}

for sample in abund.columns:
    if sample in sample_to_host:
        host = sample_to_host[sample]
        if host in host_samples:
            host_samples[host].append(sample)

for host in host_samples:
    print(f"{host}: {len(host_samples[host])} samples")


# COMPUTE HOST MEAN ABUNDANCE


host_means = {}

for host, samples in host_samples.items():
    if samples:
        host_means[host] = abund[samples].mean(axis=1)


# LOAD SGB → SPECIES MAPPING


print("Loading SGB taxonomy mapping...")
sgb_tax = {}

with bz2.open(SGB_TAXONOMY_FILE, "rt") as f:
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) < 2:
            continue

        sgb_id = parts[0].replace("SGB", "")
        taxonomy = parts[1]

        # Extract species name
        if "s__" in taxonomy:
            species = taxonomy.split("s__")[-1]
        else:
            species = taxonomy

        # Make labels clean for visualization
        species = species.replace("_", " ")

        sgb_tax[sgb_id] = species


# iTOL SPECIES LABEL DATASET


label_file = os.path.join(OUTPUT_DIR, "iTOL_species_labels.txt")

with open(label_file, "w") as f:
    f.write("DATASET_TEXT\n")
    f.write("SEPARATOR TAB\n")
    f.write("DATASET_LABEL\tSpecies\n")
    f.write("COLOR\t#000000\n")
    f.write("DATA\n")

    for leaf in tree_ids:
        species = sgb_tax.get(leaf, "Unknown")
        f.write(f"{leaf}\t{species}\n")

print("Species label dataset written.")


# HOST BAR DATASETS

for host in HOST_COLORS:

    bar_file = os.path.join(OUTPUT_DIR, f"iTOL_{host}_mean_bar.txt")

    with open(bar_file, "w") as f:

        f.write("DATASET_SIMPLEBAR\n")
        f.write("SEPARATOR TAB\n")
        f.write(f"DATASET_LABEL\t{host}_mean_abundance\n")
        f.write(f"COLOR\t{HOST_COLORS[host]}\n")
        f.write("WIDTH\t50\n")
        f.write("MARGIN\t5\n")
        f.write("SHOW_INTERNAL\t0\n")
        f.write("DATA\n")

        if host in host_means:
            for leaf in tree_ids:
                value = host_means[host].get(leaf, 0)
                f.write(f"{leaf}\t{value}\n")
        else:
            for leaf in tree_ids:
                f.write(f"{leaf}\t0\n")

print("Host bar datasets written.")


# VALIDATION REPORT


print("\nValidation summary:")
print(f"Tree leaves: {len(tree_id_set)}")
print(f"Abundance IDs: {len(abund_id_set)}")
print(f"Matched IDs: {len(common_ids)}")

print("\nUpload this tree to iTOL:")
print(TREE_FILE)

print("\nUpload these annotation files:")
print(label_file)

for host in HOST_COLORS:
    print(os.path.join(OUTPUT_DIR, f"iTOL_{host}_mean_bar.txt"))

print("\nAnnotation generation completed successfully.")
