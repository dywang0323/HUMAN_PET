# This script is used to filter out the records their prevalence is lower than 10% 

import pandas as pd

# Load the TSV file
file_path = "/ourdisk/hpc/nullspace/dywang/dont_archive/PET/Data/Merge/ecs_Human.tsv"  # Update this with your actual file path
df = pd.read_csv(file_path, sep="\t")

# Ensure the first column is not an index-like column
df.reset_index(drop=True, inplace=True)

# Calculate the total number of columns (excluding the first column)
total_columns = df.shape[1] - 1

# Compute prevalence: percentage of non-zero values in each row
non_zero_counts = (df.iloc[:, 1:] != 0).sum(axis=1)
prevalence = (non_zero_counts / total_columns) * 100

# Apply filtering conditions:
# 1. Prevalence > 10%
# 2. At least 10% of the columns must have nonzero values
filtered_df = df[(prevalence > 10) & (non_zero_counts >= 0.1 * total_columns)]

# Drop the first column if it's just an index-like series
filtered_df = filtered_df.iloc[:, 1:]

# Save the filtered data to a new TSV file
filtered_file_path = "/ourdisk/hpc/nullspace/dywang/dont_archive/PET/Data/Merge/Prevalence_filter/Reaction/Human_10_ECs.tsv"
filtered_df.to_csv(filtered_file_path, sep="\t", index=False)

# Display the result path
print(f"Filtered data saved to: {filtered_file_path}")
