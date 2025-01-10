# %%
import os

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from mpmath import timing


def get_timing_filepaths(directory):
    filepaths = []
    for root, _, files in os.walk(directory):
        for file in files:
            if "timing_data" in file:
                filepaths.append(os.path.join(root, file))
    return filepaths

def get_lungseg3_tre_filepaths(directory):
    filepaths = []
    for root, _, files in os.walk(directory):
        for file in files:
            if "lungseg3" in file and "TRE" in file:
                filepaths.append(os.path.join(root, file))
    return filepaths

def get_par003_tre_filepaths(directory):
    filepaths = []
    for root, _, files in os.walk(directory):
        for file in files:
            if "dataPar0003" in file and "TRE" in file:
                filepaths.append(os.path.join(root, file))
    return filepaths

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

# %% FILE PATHS
root = r"C:\Users\gimes\OneDrive\MAIA\3_UdG\classes\MIRA\lab\labfinalboss\data"
files = [
    rf'{root}\TRE_dataPar0003_lungseg3_1400clip_all_params.txt',
    # rf'{root}\Param0007.MI.Coarse.Bspline_tuned.csv',
    # rf'{root}\Param0008.csv',
    # rf'{root}\Param0011.csv',
    # rf'{root}\Param0015.csv',
    # rf'{root}\Param0016.csv'
]

root = r"C:\Users\gimes\Desktop\MIRA IMAGES"
files = [
    rf'{root}\TRE_dataPar0003_lungseg3_1400clip_all_params.txt',
    rf'{root}\Param0007.MI.Coarse.Bspline_tuned.csv',
    rf'{root}\Param0008.csv',
    rf'{root}\Param0011.csv',
    rf'{root}\Param0015.csv',
    rf'{root}\Param0016.csv'
]

results_directory = r'C:\Users\gimes\Desktop\MIRA IMAGES\pars'

files = get_lungseg3_tre_filepaths(results_directory)
for path in files:
    print(path)

# Read data from each file into a dataframe
dfs = []
for file in files:
    if file.endswith('.txt'):
        df = unpack_txt_file(file)
    else:
        df = pd.read_csv(file)
        filename = os.path.splitext(os.path.split(file)[-1])[-2].__str__()
        filename = filename.replace("TRE_data", "")
        filename = filename.replace("_lungseg3", "")
        print(filename)
        df['FileName'] =  filename # Add a column to indicate the source file

        # Calculate the mean for each column and append as a new row
        mean_row = df.mean(numeric_only=True).to_frame().T
        mean_row['caseName'] = 'mean'
        mean_row['FileName'] = filename
        df = pd.concat([df, mean_row], ignore_index=True)

    dfs.append(df)

# Concatenate the dataframes
combined_df = pd.concat(dfs, ignore_index=True)
combined_df

# %% CELL VALUE MANIPULATIONS
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

#%% SAVING RESULTS
output_file_path = rf'{root}\combined_result_df.csv'
combined_df.to_csv(output_file_path, index=False)

# %%  PLOTTING SETUP
# Define a custom color palette
unique_case_names = combined_df['caseName'].unique()
palette = sns.color_palette("colorblind", len(unique_case_names))
palette_dict = {case_name: palette[i] for i, case_name in enumerate(unique_case_names)}
palette_dict['mean'] = 'black'
# Plot each column separately
columns_to_plot = ['voxel', 'mm', 'NCC', 'NGC']

# %% BOXPLOT FOR EACH PARAMETERS SEGMENTAION RESULTS
# Filter the DataFrame to include only the desired cases
cases_to_plot = ['copd1', 'copd2', 'copd3', 'copd4']
fs = (10, 6)
filtered_df = combined_df[combined_df['caseName'].isin(cases_to_plot)]

# Define a custom color palette
unique_case_names = combined_df['caseName'].unique()
palette = sns.color_palette("deep", len(unique_case_names))
palette_dict = {case_name: palette[i] for i, case_name in enumerate(unique_case_names)}
palette_dict['mean'] = 'black'

# Plot each column separately
columns_to_plot = ['voxel', 'mm', 'NCC', 'NGC']
for column in columns_to_plot:
    plt.figure(figsize=fs)
    boxplot = sns.boxplot(
        x='FileName',
        y=column,
        data=filtered_df,
        hue='FileName',
        palette="deep",
        # palette=palette_dict
    )

    plt.title(f'TRE ({column}) by Parameter')
    plt.xticks(rotation=90)
    plt.xlabel("")

    plt.ylabel(column)
    # plt.legend(title='Case Name', loc='upper left')
    plt.tight_layout()
    plt.show()

#%% SCATTERPLOT FOR EACH PARAMETERS SEGMENTAION RESULTS
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

    plt.title(f'TRE ({column}) by Parameter')
    plt.xlabel("")
    plt.xticks(rotation=90)
    plt.ylabel(column)
    plt.legend(title='Case Name', loc='upper left')
    plt.tight_layout()
    plt.show()


#%% EXTRACTING DFs FOR EACH METRIC AND PARAMETERS
columns_to_pivot = ['NCC', 'NGC', 'mm', 'voxel']
pivoted_dfs = {}

for column in columns_to_pivot:
    pivoted_df = combined_df.pivot(index='caseName', columns='FileName', values=column)
    pivoted_dfs[column] = pivoted_df

# Access the separate DataFrames
ncc_df = pivoted_dfs['NCC']
ngc_df = pivoted_dfs['NGC']
mm_df = pivoted_dfs['mm']
voxel_df = pivoted_dfs['voxel']


# Example: print the NCC DataFrame
print(ncc_df)

ncc_df.to_csv(rf'{root}\result_ncc_df.csv', index=True)
ngc_df.to_csv(rf'{root}\result_ngc_result_df.csv', index=True)
mm_df.to_csv(rf'{root}\result_mm_df.csv', index=True)
voxel_df.to_csv(rf'{root}\result_voxel_df.csv', index=True)


# %% LUNG SEGMENTATION COMPARISON
par003_filepaths = get_par003_tre_filepaths(fr"{root}\pars")
print(par003_filepaths)

# Read data from each file into a dataframe
dfs = []
for file in par003_filepaths:
    df = pd.read_csv(file)
    filename = os.path.splitext(os.path.split(file)[-1])[-2].__str__()
    filename = filename.replace("TRE_data", "")
    print(filename)
    df['FileName'] =  filename # Add a column to indicate the source file

    # Calculate the mean for each column and append as a new row
    mean_row = df.mean(numeric_only=True).to_frame().T
    mean_row['caseName'] = 'mean'
    mean_row['FileName'] = filename
    df = pd.concat([df, mean_row], ignore_index=True)

    dfs.append(df)
dfs = pd.concat(dfs, ignore_index=True)

# print(dfs)
column = "mm"
mm_df_lungseg = dfs.pivot(index='caseName', columns='FileName', values=column)
print(mm_df_lungseg)
mm_df_lungseg.to_csv(fr"{root}\results_Par0003_lungmasks.csv", index=True)


# %% PARAMETER COMPARISON
column = "NCC"
mm_df = combined_df.pivot(index='caseName', columns='FileName', values=column)
print(mm_df)
mm_df.to_csv(fr"{root}\result_{column}_by_cases.csv", index=True)


# %% TIMING DATA LOADING
timing_filepaths = get_timing_filepaths(results_directory)

dfs = []
for filename in timing_filepaths:
    df = pd.read_csv(filename)
    filename = os.path.splitext(os.path.split(filename)[-1])[-2].replace("timing_data","")
    df["FileName"] = filename

    mean_row = df.mean(numeric_only=True).to_frame().T
    mean_row["FileName"] = filename
    mean_row["caseName"] = "mean"
    df = pd.concat([df, mean_row], ignore_index=True)


    # print(mean_row)
    # print(mean)
    print(df)

    dfs.append(df)

timing_df = pd.concat(dfs)
print(timing_df)

# %% TIMINGS FOR PARAMETERS

par003_times = timing_df[timing_df['FileName'].str.contains('Par0003')]

# Print the save times
print(par003_times)

# Assuming par003_times is your DataFrame containing the timing data for Par0003
# Pivot the DataFrame to have distinct filenames as columns and cases as rows
pivoted_df = par003_times.pivot(index='caseName', columns='FileName', values='saveTime')

# Calculate the mean for each column and append as a new row
mean_row = pivoted_df.mean().to_frame().T
# mean_row.index = ['mean']
# pivoted_df = pd.concat([pivoted_df, mean_row])

# Print the resulting DataFrame
print(pivoted_df)
pivoted_df.save(rf'{root}\result_par003_times.csv', index=True)

# %%
lungseg3_times = timing_df[timing_df['FileName'].str.contains('lungseg3')]

# Print the save times
print(lungseg3_times)

# Assuming par003_times is your DataFrame containing the timing data for Par0003
# Pivot the DataFrame to have distinct filenames as columns and cases as rows
pivoted_df = lungseg3_times.pivot(index='caseName', columns='FileName', values='registrationTime')

# Calculate the mean for each column and append as a new row
mean_row = pivoted_df.mean().to_frame().T
# mean_row.index = ['mean']
# pivoted_df = pd.concat([pivoted_df, mean_row])

# Print the resulting DataFrame
print(pivoted_df)
pivoted_df.save(rf'{root}\result_lungmask3_times.csv', index=True)
