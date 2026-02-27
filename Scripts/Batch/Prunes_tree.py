# to run this script:
# python /Prunes_tree.py \
# --humann /HUMAnN4_Species_Stratified_preva.tsv \
# --tree /mpa_vJun23_CHOCOPhlAnSGB_202403.nwk \
# --out-subtree /GO_SGB_subtree.nwk \
# --out-matrix /GO_SGB_abundance.tsv

#!/usr/bin/env python3

import argparse
import re
from skbio import TreeNode


SGB_RE = re.compile(r"SGB(\d+)")

# Aggregate HUMAnN abundances per SGB

def aggregate_by_sgb(humann_file):
    sgb_matrix = {}

    with open(humann_file, "r", encoding="utf-8", errors="replace") as fin:
        header = fin.readline().strip().split("\t")
        header_samples = header[1:]

        for line in fin:
            if "|" not in line:
                continue

            parts = line.rstrip("\n").split("\t")
            first_col = parts[0]

            match = SGB_RE.search(first_col)
            if not match:
                continue

            sgb_digits = match.group(1)

            numeric_vals = [
                float(x) if x not in ("", "NA") else 0.0
                for x in parts[1:]
            ]

            if sgb_digits not in sgb_matrix:
                sgb_matrix[sgb_digits] = [0.0] * len(numeric_vals)

            row = sgb_matrix[sgb_digits]
            for i, v in enumerate(numeric_vals):
                row[i] += v

    return sgb_matrix, header_samples

# Extract subtree using numeric SGB tree

def generate_subtree(tree_file, sgb_matrix,
                     header_samples,
                     out_tree, out_matrix):

    print("Loading numeric SGB tree...")
    tree = TreeNode.read(tree_file)

    # Tree tips are pure digits
    tree_tip_ids = {tip.name for tip in tree.tips()}

    humann_sgbs = set(sgb_matrix.keys())
    overlap = sorted(humann_sgbs.intersection(tree_tip_ids))

    print(f"SGBs in HUMAnN: {len(humann_sgbs)}")
    print(f"SGB tips in tree: {len(tree_tip_ids)}")
    print(f"Overlap: {len(overlap)}")

    if not overlap:
        raise ValueError("No matching SGB IDs found between HUMAnN and tree.")

    # Extract subtree
    print("Extracting subtree using shear()...")
    subtree = tree.shear(set(overlap))
    subtree.write(out_tree)

    # Write abundance matrix
    with open(out_matrix, "w", encoding="utf-8") as f:
        f.write("SGB_ID\t" + "\t".join(header_samples) + "\n")

        for sgb in overlap:
            values = sgb_matrix[sgb]
            row = [f"SGB{sgb}"] + [f"{v:.6f}" for v in values]
            f.write("\t".join(row) + "\n")

    print("Subtree and matrix successfully written.")


# Main

def main():
    parser = argparse.ArgumentParser(
        description="Fast subtree extraction using numeric CHOCOPhlAn SGB tree."
    )

    parser.add_argument("--humann", required=True)
    parser.add_argument("--tree", required=True)
    parser.add_argument("--out-subtree", default="subtree_SGB.nwk")
    parser.add_argument("--out-matrix", default="SGB_abundance_matrix.tsv")

    args = parser.parse_args()

    print("Aggregating HUMAnN abundances...")
    sgb_matrix, header_samples = aggregate_by_sgb(args.humann)

    generate_subtree(
        args.tree,
        sgb_matrix,
        header_samples,
        args.out_subtree,
        args.out_matrix
    )

    print("Done.")


if __name__ == "__main__":
    main()
