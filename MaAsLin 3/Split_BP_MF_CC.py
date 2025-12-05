import pandas as pd
import os

# Load the TSV file
input_file = '/n/holystore01/LABS/huttenhower_lab/Users/Dongyu/PET_Project/Data/HumANn4_profile/GO/Informative_GO_NEW_1015/Unstraitified_HUMANn4_FUGAsseM/MaAsLin3/GO_Overall/cat_vs_human_human_ref/all_results.tsv'
df = pd.read_csv(input_file, sep='\t')

# Identify the GO term column
go_column = df.columns[0]

# Create masks for each GO category
bp_mask = df[go_column].str.contains(r'\[BP\]')
mf_mask = df[go_column].str.contains(r'\[MF\]')
cc_mask = df[go_column].str.contains(r'\[CC\]')

# Split the dataframe accordingly
df_bp = df[bp_mask].copy()
df_mf = df[mf_mask].copy()
df_cc = df[cc_mask].copy()

# Rename the first column
df_bp.rename(columns={go_column: 'GO_BP'}, inplace=True)
df_mf.rename(columns={go_column: 'GO_MF'}, inplace=True)
df_cc.rename(columns={go_column: 'GO_CC'}, inplace=True)

# Define output directory
output_dir = '/n/holystore01/LABS/huttenhower_lab/Users/Dongyu/PET_Project/Data/HumANn4_profile/GO/Informative_GO_NEW_1015/Unstraitified_HUMANn4_FUGAsseM/MaAsLin3/GO_Overall/cat_vs_human_human_ref/Split'
os.makedirs(output_dir, exist_ok=True)

# Save the files to the output directory
df_bp.to_csv(os.path.join(output_dir, 'GO_BP.tsv'), sep='\t', index=False)
df_mf.to_csv(os.path.join(output_dir, 'GO_MF.tsv'), sep='\t', index=False)
df_cc.to_csv(os.path.join(output_dir, 'GO_CC.tsv'), sep='\t', index=False)

print("Files saved in the 'GO_output' directory.")
