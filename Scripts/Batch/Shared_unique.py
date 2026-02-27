import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

# CONFIGURATION

threshold = 1e-6
pseudocount = 1e-6

metadata_file = '/n/holystore01/LABS/huttenhower_lab/Users/Dongyu/PET_Project/Data/HumANn4_profile/S15_metadata_consist.csv'
abundance_file = '/n/holystore01/LABS/huttenhower_lab/Users/Dongyu/PET_Project/Data/HumANn4_profile/GO/Informative_GO_NEW_1015/HUMANn4_Informative_GO/Split/GO_MF.tsv'

output_dir = '/n/holystore01/LABS/huttenhower_lab/Users/Dongyu/PET_Project/Data/HumANn4_profile/GO/Informative_GO_NEW_1015/HUMANn4_Informative_GO/Shared_unique/'
os.makedirs(output_dir, exist_ok=True)

x_limits = (-6, 0)
dot_size = 8
axis_linewidth = 1.5

group_colors = {
    'Shared': '#ffffcc',
    'DogCat_Shared': '#d0ffd0',
    'Unique': '#e6ccff'
}

species_palette = {
    "cat": "indianred",
    "dog": "#E69F00",
    "human": "#56B4E0"
}


# LOAD DATA

metadata_df = pd.read_csv(metadata_file)
main_df = pd.read_csv(abundance_file, sep='\t')

feature_col = main_df.columns[0]
sample_columns = main_df.columns[1:]
df = main_df.set_index(feature_col)


# MAP SPECIES SAMPLES

dog_samples = metadata_df.loc[metadata_df['species']=='dog','sample_id_metaphlan']
cat_samples = metadata_df.loc[metadata_df['species']=='cat','sample_id_metaphlan']
human_samples = metadata_df.loc[metadata_df['species']=='human','sample_id_metaphlan']

dog_cols = [c for c in sample_columns if c in dog_samples.values]
cat_cols = [c for c in sample_columns if c in cat_samples.values]
human_cols = [c for c in sample_columns if c in human_samples.values]


# CALCULATE MEANS 

df['dog_mean'] = df[dog_cols].mean(axis=1)
df['cat_mean'] = df[cat_cols].mean(axis=1)
df['human_mean'] = df[human_cols].mean(axis=1)

df['mean_abundance'] = df[sample_columns].mean(axis=1)
df = df[df['mean_abundance'] > threshold]


# SELECT FEATURES 

# Shared
shared = df[
    (df['dog_mean'] > threshold) &
    (df['cat_mean'] > threshold) &
    (df['human_mean'] > threshold)
].copy()

shared['avg_mean'] = shared[['dog_mean','cat_mean','human_mean']].mean(axis=1)
top_shared = shared.sort_values('avg_mean', ascending=False).head(5)

# DogCat only
dogcat = df[
    (df['dog_mean'] > threshold) &
    (df['cat_mean'] > threshold) &
    (df['human_mean'] <= threshold)
].copy()

dogcat['avg_mean'] = dogcat[['dog_mean','cat_mean']].mean(axis=1)
top_dogcat = dogcat.sort_values('avg_mean', ascending=False).head(5)

# Unique per species
dog_unique = df[
    (df['dog_mean'] > threshold) &
    (df['cat_mean'] <= threshold) &
    (df['human_mean'] <= threshold)
].sort_values('dog_mean', ascending=False).head(5)

cat_unique = df[
    (df['cat_mean'] > threshold) &
    (df['dog_mean'] <= threshold) &
    (df['human_mean'] <= threshold)
].sort_values('cat_mean', ascending=False).head(5)

human_unique = df[
    (df['human_mean'] > threshold) &
    (df['dog_mean'] <= threshold) &
    (df['cat_mean'] <= threshold)
].sort_values('human_mean', ascending=False).head(5)


# PREPARE FOR PLOTTING


def melt_group(data, species_cols, group_name):
    melted = data[species_cols].reset_index().melt(
        id_vars=[feature_col],
        value_vars=species_cols,
        var_name='Species',
        value_name='Abundance'
    )
    melted['Species'] = melted['Species'].str.replace('_mean','',regex=False)
    melted['x_value'] = np.log10(melted['Abundance'] + pseudocount)
    melted['GroupType'] = group_name
    return melted

shared_plot = melt_group(top_shared,
                         ['dog_mean','cat_mean','human_mean'],
                         'Shared')

dogcat_plot = melt_group(top_dogcat,
                         ['dog_mean','cat_mean'],
                         'DogCat_Shared')

# Unique plotting
unique_frames = []

if not dog_unique.empty:
    tmp = dog_unique[['dog_mean']].rename(columns={'dog_mean':'Abundance'})
    tmp['Species'] = 'dog'
    tmp['GroupType'] = 'Unique'
    unique_frames.append(tmp)

if not cat_unique.empty:
    tmp = cat_unique[['cat_mean']].rename(columns={'cat_mean':'Abundance'})
    tmp['Species'] = 'cat'
    tmp['GroupType'] = 'Unique'
    unique_frames.append(tmp)

if not human_unique.empty:
    tmp = human_unique[['human_mean']].rename(columns={'human_mean':'Abundance'})
    tmp['Species'] = 'human'
    tmp['GroupType'] = 'Unique'
    unique_frames.append(tmp)

if unique_frames:
    unique_plot = pd.concat(unique_frames).reset_index()
    unique_plot['x_value'] = np.log10(unique_plot['Abundance'] + pseudocount)
else:
    unique_plot = pd.DataFrame()

combined_df = pd.concat([shared_plot, dogcat_plot, unique_plot], ignore_index=True)


# SORT WITHIN EACH GROUP

shared_order = top_shared.index.tolist()
dogcat_order = top_dogcat.index.tolist()
dog_unique_order = dog_unique.index.tolist()
cat_unique_order = cat_unique.index.tolist()
human_unique_order = human_unique.index.tolist()

feature_order = (
    shared_order +
    dogcat_order +
    dog_unique_order +
    cat_unique_order +
    human_unique_order
)

combined_df[feature_col] = pd.Categorical(
    combined_df[feature_col],
    categories=feature_order,
    ordered=True
)


# PLOT

per_row_height = 0.55
fig_height = per_row_height * len(feature_order)

max_label_len = max(len(str(x)) for x in feature_order)
left_margin_inch = 6 + 0.08 * max(0, max_label_len - 40)
right_margin_inch = 2
plot_width_inch = 3

fig_width = left_margin_inch + plot_width_inch + right_margin_inch
left_frac = left_margin_inch / fig_width
right_frac = 1 - (right_margin_inch / fig_width)

fig, ax = plt.subplots(figsize=(fig_width, fig_height))

# Background shading
for group in ['Shared','DogCat_Shared','Unique']:
    group_feats = combined_df[
        combined_df['GroupType']==group
    ][feature_col].unique()
    for feat in group_feats:
        pos = feature_order.index(feat)
        ax.axhspan(pos-0.5, pos+0.5,
                   color=group_colors[group],
                   alpha=0.3)

sns.stripplot(
    data=combined_df,
    y=feature_col,
    x='x_value',
    hue='Species',
    dodge=True,
    size=dot_size,
    palette=species_palette,
    orient='h',
    ax=ax
)

ax.set_xlim(x_limits)
ax.set_xlabel('log10(Mean Abundance)')
ax.set_ylabel('')

for spine in ax.spines.values():
    spine.set_linewidth(axis_linewidth)

ax.legend(title='Species',
          bbox_to_anchor=(1.02,1),
          loc='upper left',
          frameon=False)

plt.subplots_adjust(
    left=left_frac,
    right=right_frac,
    top=0.98,
    bottom=0.05
)

plt.savefig(
    os.path.join(output_dir,'Integrated_Shared_Unique_final.tiff'),
    dpi=300,
    format='tiff'
)

plt.close()

print("Final sorted shared & unique figure saved.")
