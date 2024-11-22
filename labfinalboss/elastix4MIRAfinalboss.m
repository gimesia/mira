function read_raw_binary_to_nifti(filename, shape, spacing)
    % Read binary data
    fid = fopen(filename, 'rb');
    data = fread(fid, inf, 'int16');
    fclose(fid);

    % Reshape data
    data = reshape(data, shape);

    % Create NIfTI image
    img = make_nii(data, -spacing);

    % Save the NIfTI image
    save_nii(img, strrep(filename, '.img', '_fromMATLAB.nii.gz'));
end

% Paths
elastixBasePath = "C:\Users\gimes\OneDrive\MAIA\3_UdG\classes\MIRA\lab\lab2\";
basePath = "C:\Users\gimes\OneDrive\MAIA\3_UdG\classes\MIRA\lab\labfinalboss\";
dataPath = "C:\Users\gimes\OneDrive\MAIA\3_UdG\classes\MIRA\lab\labfinalboss\data\";

rawImgFiles = dir(fullfile(dataPath, "**", "*.img")); % only raw
rawImgFiles = {rawImgFiles.name}

for i = rawImgFiles
    
end

imgFiles = dir(fullfile(dataPath, "**", "*.nii.gz")); % only nifti
imgFiles = {imgFiles.name};  % Get all filenames in a cell array
imgFiles = strcat(strtok(imgFiles, '_'), '\', imgFiles);

pointFiles = dir(fullfile(dataPath, "**", "*.txt")); % only txt
pointFiles= {pointFiles.name};  % Get all filenames in a cell array
pointFiles = strcat(strtok(pointFiles, '_'), '\', pointFiles);

cases = unique(strtok(pointFiles, "\"));

% Path to affine parameter file
% parameterFileAffine = elastixBasePath + "elastix_model_zoo\models\Par0010\Par0010affine.txt";
% Path to B-spline parameter file
parameterFileBSpline = "C:\Users\gimes\OneDrive\MAIA\3_UdG\classes\" + ...
    "MIRA\lab\lab2\elastix_model_zoo\models\Par0016\" + ...
    "Par0016.multibsplines.lung.sliding.txt";
    
% Paths to elastix and transformix executables
elastixPath = elastixBasePath + "elastix-5.0.0-win64\elastix.exe";
transformixPath = elastixBasePath + "elastix-5.0.0-win64\transformix.exe";


for caseID = 1:length(cases)
    caseName = cases{caseID};
    disp("Processing Case [" + caseName + "]")
    idxs = cellfun(@(x) contains(x, caseName), imgFiles);
    idxs = find(idxs == 1);
    caseImgs = imgFiles(idxs);
    casePoints = pointFiles(idxs);

    inhaleImgId = cellfun(@(x) contains(x, "iBH"), caseImgs);
    exhaleImgId = cellfun(@(x) contains(x, "eBH"), caseImgs);
    inhaleImg= fullfile(dataPath, caseImgs{inhaleImgId == 1});
    exhaleImg = fullfile(dataPath, caseImgs{exhaleImgId == 1});

    inhalePointsId = cellfun(@(x) contains(x, "iBH"), casePoints);
    exhalePointsId = cellfun(@(x) contains(x, "eBH"), casePoints);
    inhalePoints = fullfile(dataPath, casePoints{inhalePointsId == 1});
    exhalePoints = fullfile(dataPath, casePoints{exhalePointsId == 1});

    outputDirectory = basePath + "data\" + caseName + "\out";
    if ~exist(outputDirectory, 'dir')
        mkdir(outputDirectory);
    end


    disp("Registration with elastix ([exhale] -> [inhale]): " + caseName);
    elastixCommand = sprintf('"%s" -f "%s" -m "%s" -p "%s" -out "%s"', ...
                      elastixPath, exhaleImg, inhaleImg, parameterFileBSpline, outputDirectory);
    [status, result] = system(elastixCommand);
    if status == 0
        disp("Registration with elastix completed.")
    else
        disp("There has been an error, see log!")
    end
end



% Moving
movingImage = basePath + "result-atlas2\atlas_intensity.nii.gz";


for caseId = 1:length(shortFileNames)
    realCurrentMovingImage = string(shortFileNames(caseId));

    baseOutputDirectory = basePath + "output2\" + realCurrentMovingImage + "\";
    disp("Processing Image No. " + realCurrentMovingImage);
    % Fixed
    fixedImage = basePath + "test-set\testing-images\" + realCurrentMovingImage + ".nii.gz";
    fixedMask = basePath + "test-set\testing-mask\" + realCurrentMovingImage + "_1C.nii.gz";
    fixedLabels = basePath + "test-set\testing-labels\" + realCurrentMovingImage + "_3C.nii.gz";
    
    
    % Create output directory for affine registration
    outputDirectoryAffine = baseOutputDirectory + "registration";
    if ~exist(outputDirectoryAffine, 'dir')
        mkdir(outputDirectoryAffine);
    end
    
  
    % Start affine registration process
    disp('RUNNING REGISTRATION PIPELINE FOR CASE: ' + realCurrentMovingImage)
    elastixCommand = sprintf('"%s" -f "%s" -m "%s" -fMask "%s" -p "%s" -p "%s" -out "%s"', ...
                      elastixPath, fixedImage, movingImage, fixedMask, parameterFileAffine, parameterFileBSpline, outputDirectoryAffine);
    [status, result] = system(elastixCommand);
    
    if status == 0
        % Affine registration succeeded, renaming output file
        elastixResultFile = fullfile(outputDirectoryAffine, "result.0.nii.gz");
        newElastixResultFileName = fullfile(baseOutputDirectory, realCurrentMovingImage + "_registered_image_" + "_affine.nii.gz");
        if isfile(elastixResultFile)
            movefile(elastixResultFile, newElastixResultFileName);
        end
        affineTransformationParamFile = fullfile(outputDirectoryAffine, 'TransformParameters.0.txt');
        
        % Affine registration succeeded, renaming output file
        elastixResultFile = fullfile(outputDirectoryAffine, "result.1.nii.gz");
        newElastixResultFileName = fullfile(baseOutputDirectory, realCurrentMovingImage + "_registered_image_" + "_bspline.nii.gz");
        
        if isfile(elastixResultFile)
            movefile(elastixResultFile, newElastixResultFileName);
        end
        transformationParamFile = fullfile(outputDirectoryAffine, 'TransformParameters.1.txt')
        
        disp('Elastix registration completed successfully.');
    else
        disp('Error running Elastix');
        disp(result);
        return; % Exit if affine registration fails
    end
    
    
    % Modifying transformation parameter files
    % List of transformation parameter files to modify
    transformParamFiles = {affineTransformationParamFile, transformationParamFile};
    for i = 1:length(transformParamFiles)
        currentFile = transformParamFiles{i};
        
        % Check if the parameter file exists
        if isfile(currentFile)
            % Read the file content
            fileText = fileread(currentFile);
    
            % Replace the interpolator with "FinalNearestNeighborInterpolator"
            fileText = regexprep(fileText, ...
                '\(ResampleInterpolator\s+"FinalBSplineInterpolator"\)', ...
                '(ResampleInterpolator "FinalNearestNeighborInterpolator")');
            
            % Replace ResultImagePixelType with "float"
            fileText = regexprep(fileText, ...
                '\(ResultImagePixelType\s+".*?"\)', ...
                '(ResultImagePixelType "float")');
            
            % Write the modified content back to the file
            fileID = fopen(currentFile, 'w');
            fprintf(fileID, '%s', fileText);
            fclose(fileID);
        else
            disp(['Transformation parameter file not found: ', currentFile]);
        end
    end
    
    % Setup output directories for Transformix
    outputTransformixLabels = fullfile(baseOutputDirectory, "transformix");
    if ~exist(outputTransformixLabels, 'dir')
        mkdir(outputTransformixLabels);
    end
    
    for i = 0:3
        movingLabel = basePath + "\result-atlas2\atlas_probability_label_" + string(num2str(i)) + ".nii.gz"
        % Run Transformix for mask transformation
        transformixMaskCommand = sprintf('"%s" -in "%s" -tp "%s" -out "%s"', ...
                             transformixPath, movingLabel, transformationParamFile, outputTransformixLabels);
        [statusMask, resultMask] = system(transformixMaskCommand);
    
        % Rename transformed mask file
        transformedMaskFile = fullfile(outputTransformixLabels, "result.nii.gz");
        newMaskFileName = fullfile(baseOutputDirectory, realCurrentMovingImage + "_transformed_atlas_label_" + string(num2str(i)) + ".nii.gz");
        if isfile(transformedMaskFile)
            movefile(transformedMaskFile, newMaskFileName);
        end
        if statusMask == 0
            disp("Probability-map for label " + string(i) + " transformed successfully")
        else
            disp("Probability-map for label " + string(i) + " transformation failed")
        end
    end
end