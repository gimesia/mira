% Paths
elastixBasePath = "C:\Users\gimes\OneDrive\MAIA\3_UdG\classes\MIRA\lab\lab2\";
basePath = "C:\Users\gimes\OneDrive\MAIA\3_UdG\classes\MIRA\lab\labfinalboss\";
dataPath = "C:\Users\gimes\OneDrive\MAIA\3_UdG\classes\MIRA\lab\labfinalboss\data\";

rawImgFiles = dir(fullfile(dataPath, "**", "*BHCT.img")); % only raw
rawImgFiles = {rawImgFiles.name};
rawImgFiles = strcat(strtok(rawImgFiles, '_'), '\', rawImgFiles);

pointFiles = dir(fullfile(dataPath, "**", "*xyz_r1.txt")); % only points txt
pointFiles= {pointFiles.name};  % Get all filenames in a cell array
pointFiles = strcat(strtok(pointFiles, '_'), '\', pointFiles);

cases = unique(strtok(pointFiles, "\"));

% Define the data
Label = {'copd1', 'copd2', 'copd3', 'copd4', 'copd5', 'copd6', 'copd7', 'copd8', 'copd9', 'copd10'};
ImageDims = {[512, 512, 121], [512, 512, 102], [512, 512, 126], [512, 512, 126], [512, 512, 131], ...
             [512, 512, 119], [512, 512, 112], [512, 512, 115], [512, 512, 116], [512, 512, 135]};
Spacing = {[0.625, 0.625, 2.5], [0.645, 0.645, 2.5], [0.652, 0.652, 2.5], [0.59, 0.59, 2.5], ...
           [0.647, 0.647, 2.5], [0.633, 0.633, 2.5], [0.625, 0.625, 2.5], [0.586, 0.586, 2.5], ...
           [0.664, 0.664, 2.5], [0.742, 0.742, 2.5]};
NumFeatures = [773, 612, 1172, 786, 1029, 633, 575, 791, 447, 480];
Displacement = {'25.90 (11.57)', '21.77 (6.46)', '12.29 (6.39)', '30.90 (13.49)', '30.90 (14.05)', ...
                '28.32 (9.20)', '21.66 (7.66)', '25.57 (13.61)', '14.84 (10.01)', '22.48 (10.64)'};
NumRepeats = repmat({'150/3'}, 10, 1);
Observers = {'0.65 (0.73)', '1.06 (1.51)', '0.58 (0.87)', '0.71 (0.96)', '0.65 (0.87)', ...
             '1.06 (2.38)', '0.65 (0.78)', '0.96 (3.07)', '1.04 (2.54)', '0.87 (1.65)'};

% Create the table
data = table(Label', ImageDims', Spacing', NumFeatures', Displacement', NumRepeats, Observers', ...
                 'VariableNames', {'Label', 'ImageDims', 'Spacing', 'NumFeatures', 'Displacement', 'NumRepeats', 'Observers'});
