# This script is used to filter out the records their prevalence is lower than 10% 

import pandas as pd

# Load the TSV file
file_path = "/ourdisk/hpc/nullspace/dywang/dont_archive/PET/Data/Merge/ecs_Human.tsv" 
df = pd.read_csv(file_path, sep="\t")

# Identify the first column (assumed to be pathway or ID column)
id_column = df.columns[0]

# Calculate the total number of numerical data columns (excluding the first column)
total_columns = df.shape[1] - 1

# Compute prevalence: percentage of non-zero values in each row (excluding the ID column)
non_zero_counts = (df.iloc[:, 1:] != 0).sum(axis=1)
prevalence = (non_zero_counts / total_columns) * 100

# Apply filtering condition: Prevalence > 10%
filtered_df = df[prevalence > 10]

# Ensure the first column (ID) is retained
filtered_file_path = "/ourdisk/hpc/nullspace/dywang/dont_archive/PET/Data/Merge/Prevalence_filter/Reaction/Human_10_ECs.tsv" 
filtered_df.to_csv(filtered_file_path, sep="\t", index=False)

# Display success message
print(f"Filtered data saved to: {filtered_file_path}")
