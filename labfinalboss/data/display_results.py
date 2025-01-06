#%%
import os

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sympy import false


def unpack_txt_file(file_path):
    # Read the data from the txt file
    df = pd.read_csv(file_path)

    # Assuming df is your dataframe
    grouped = df.groupby(df.FileName)

    # Create a list of dataframes
    list_of_dfs = [group for _, group in grouped]

    for i, dframe in enumerate(list_of_dfs):
        # Calculate the mean for each column and append as a new row
        mean_row = dframe.mean(numeric_only=True).to_frame().T
        mean_row['caseName'] = 'mean'
        mean_row['FileName'] = dframe.FileName.iloc[0]
        dframe = pd.concat([dframe, mean_row], ignore_index=True)
        list_of_dfs[i] = dframe

    return pd.concat(list_of_dfs, ignore_index=True)

# File paths
root = r"C:\Users\gimes\OneDrive\MAIA\3_UdG\classes\MIRA\lab\labfinalboss\data"
files = [
    rf'{root}\TRE_dataPar0003_lungseg3_1400clip_all_params.txt',
    rf'{root}\Param0007.MI.Coarse.Bspline_tuned.csv',
    rf'{root}\TRE_dataPar0008.csv',
    rf'{root}\TRE_dataPar0015.csv',
    rf'{root}\TRE_dataPar0016.csv'
]

# Read data from each file into a dataframe
dfs = []
for file in files:
    if file.endswith('.txt'):
        df = unpack_txt_file(file)
    else:
        df = pd.read_csv(file)
        df['FileName'] = file.split("/")[-1]  # Add a column to indicate the source file

        # Calculate the mean for each column and append as a new row
        mean_row = df.mean(numeric_only=True).to_frame().T
        mean_row['caseName'] = 'mean'
        mean_row['FileName'] = file.split("/")[-1]
        df = pd.concat([df, mean_row], ignore_index=True)

    dfs.append(df)

# Concatenate the dataframes
combined_df = pd.concat(dfs, ignore_index=True)
combined_df

#%%
# Filter out rows containing "-ug" in the FileName column
filtered_df = combined_df[combined_df['FileName'].str.contains('-fg')]
# Remove these rows from the original DataFrame
remaining_df = combined_df[~combined_df['FileName'].str.contains('-fg')]
# Append the filtered DataFrame to the start of the remaining DataFrame
final_df = pd.concat([filtered_df, remaining_df], ignore_index=True)
# Display the final DataFrame
print(final_df)

combined_df = final_df
combined_df['FileName'] = combined_df['FileName'].apply(lambda x: os.path.splitext(os.path.split(x)[-1])[0])

#%%
output_file_path = rf'{root}\combined_result_df.csv'
combined_df.to_csv(output_file_path, index=False)

#%%
# Define a custom color palette
unique_case_names = combined_df['caseName'].unique()
palette = sns.color_palette("colorblind", len(unique_case_names))
palette_dict = {case_name: palette[i] for i, case_name in enumerate(unique_case_names)}
palette_dict['mean'] = 'black'
# Plot each column separately
columns_to_plot = ['voxel', 'mm', 'NCC']
for column in columns_to_plot:
    plt.figure(figsize=(14, 10))
    sns.scatterplot(
        x='FileName',
        y=column,
        hue='caseName',
        data=combined_df,
        palette=palette_dict,
        legend='full'
    )
    # Add line plot for mean caseName
    sns.lineplot(
        x='FileName',
        y=column,
        hue='caseName',
        data=combined_df[combined_df['caseName'] == 'mean'],
        palette={'mean': palette_dict['mean']},
        legend=False,
        linestyle='--'
    )
    plt.title(f'Performance Comparison by File ({column})')
    plt.xlabel('Param Name')
    plt.xticks(rotation=90)
    plt.ylabel(column)
    plt.legend(title='Case Name', loc='upper left')
    plt.tight_layout()
    plt.show()

#%%
# Filter the DataFrame to include only the desired cases
cases_to_plot = ['copd1', 'copd2', 'copd3', 'copd4']
fs = (14, 10)
filtered_df = combined_df[combined_df['caseName'].isin(cases_to_plot)]

# Define a custom color palette
unique_case_names = combined_df['caseName'].unique()
palette = sns.color_palette("deep", len(unique_case_names))
palette_dict = {case_name: palette[i] for i, case_name in enumerate(unique_case_names)}
palette_dict['mean'] = 'black'

# Plot each column separately
columns_to_plot = ['voxel', 'mm', 'NCC']
for column in columns_to_plot:
    plt.figure(figsize=fs)
    boxplot = sns.boxplot(
        x='FileName',
        y=column,
        data=filtered_df,
        hue='FileName',
        palette="deep"
        # palette=palette_dict
    )

    plt.title(f'Boxplot of {column} by Parameter Name')
    plt.xlabel('Parameter Name')
    plt.xticks(rotation=90)
    plt.ylabel(column)
    # plt.legend(title='Case Name', loc='upper left')
    plt.tight_layout()
    plt.show()

#%%
for column in columns_to_plot:
    plt.figure(figsize=fs)
    sns.scatterplot(
        x='FileName',
        y=column,
        data=filtered_df,
        palette="deep",
        hue='caseName',
        legend='full',
        s = 100,  # Increase the size of the dots
    )
    sns.lineplot(
        x='FileName',
        y=column,
        data=combined_df[combined_df['caseName'] == 'mean'],
        palette=palette_dict,
        hue='caseName',
        legend=False,
        linestyle='--'
    )

    plt.title(f'Boxplot of {column} by Parameter')
    plt.xlabel('Parameter Name')
    plt.xticks(rotation=90)
    plt.ylabel(column)
    plt.legend(title='Case Name', loc='upper left')
    plt.tight_layout()
    plt.show()