metaData

% Path to parameter file
parameter = "Par0011"

parameterFileAffine = ...
    "C:\Users\gimes\OneDrive\MAIA\3_UdG\classes\MIRA\lab\lab2\elastix_model_zoo\models\Par0011\Parameters.Par0011.affine.txt"
parameterFileBSpline = ...
    "C:\Users\gimes\OneDrive\MAIA\3_UdG\classes\MIRA\lab\lab2\" + ...
    "elastix_model_zoo\models\Par0011\Parameters.Par0011.bspline1_s.txt"
parameterFileBSpline2 = ...
    "C:\Users\gimes\OneDrive\MAIA\3_UdG\classes\MIRA\lab\lab2\" + ...
    "elastix_model_zoo\models\Par0011\Parameters.Par0011.bspline2_s.txt"


REGISTER = 1;
PREPROCESS = 1;

% Paths to elastix and transformix executables
elastixPath = elastixBasePath + "elastix-5.0.0-win64\elastix.exe";
transformixPath = elastixBasePath + "elastix-5.0.0-win64\transformix.exe";

% To track timing for each case
timingData = struct();
TREdata = struct();

for caseID = 1:length(cases)
    caseName = string(cases{caseID});
    disp("Processing Case [" + caseName + "]");
    idxs = cellfun(@(x) contains(x, caseName), rawImgFiles);
    idxs = find(idxs == 1);
    caseImgs = rawImgFiles(idxs);
    casePoints = pointFiles(idxs);
    index = find(strcmp(data.Label, caseName));
    imageDims = data.ImageDims{index};
    imageSpacing = data.Spacing{index};

% ---------------------------------------------------------------------
% READING IN RAW IMAGES
% ---------------------------------------------------------------------
    % Measure time for reading data
    tic;
    % Moving Img
    inhaleImgId = cellfun(@(x) contains(x, "iBH"), caseImgs);
    inhaleImg = fullfile(dataPath, caseImgs{inhaleImgId == 1});
    inhaleImgNifti= read_raw_to_nifti(inhaleImg, imageDims, imageSpacing);

    % Fixed Img
    exhaleImgId = cellfun(@(x) contains(x, "eBH"), caseImgs);
    exhaleImg = fullfile(dataPath, caseImgs{exhaleImgId == 1});
    exhaleImgNifti = read_raw_to_nifti(exhaleImg, imageDims, imageSpacing);

    % Fixed Points
    inhalePointsId = cellfun(@(x) contains(x, "iBH"), casePoints);
    inhalePointsPath = fullfile(dataPath, casePoints{inhalePointsId == 1});
    inhalePoints = readmatrix(inhalePointsPath);

    % Moving Points
    exhalePointsId = cellfun(@(x) contains(x, "eBH"), casePoints);
    exhalePointsPath = fullfile(dataPath, casePoints{exhalePointsId == 1});
    exhalePoints = readmatrix(exhalePointsPath);

    readTime = toc;
    disp("Reading in completed in " + readTime + " seconds.");

    outputDirectory = basePath + "data\" + caseName + "\out";
    if ~exist(outputDirectory, 'dir')
        mkdir(outputDirectory);
    end


    % Measure time for preprocessing data
    tic;
    if PREPROCESS
% ---------------------------------------------------------------------
% PREPROCESSING
% ---------------------------------------------------------------------
        inhaleImgPreprocNifti = preprocess(inhaleImgNifti);
        exhaleImgPreprocNifti = preprocess(exhaleImgNifti);

% ---------------------------------------------------------------------
% LUNG SEGMENTATION
% ---------------------------------------------------------------------
        inhaleLungMaskNifti = inhaleImgPreprocNifti;
        inhaleLungMaskNifti.img = segmentLungs3(inhaleLungMaskNifti.img);
        exhaleLungMaskNifti = exhaleImgPreprocNifti;
        exhaleLungMaskNifti.img = segmentLungs3(exhaleLungMaskNifti.img);
    
% ---------------------------------------------------------------------
% KEYPOINT MASK GENERATION
% ---------------------------------------------------------------------
        inhalePointsMaskNifti = create_keypoints_mask_nifti(inhalePoints, inhaleImgNifti);
        exhalePointsMaskNifti = create_keypoints_mask_nifti(exhalePoints, exhaleImgNifti);
    end

% ---------------------------------------------------------------------
% INTERMEDIARY IMAGES SAVING
% ---------------------------------------------------------------------
    [~, inhaleImgFileName, ~] = fileparts(inhaleImg);
    [~, exhaleImgFileName, ~] = fileparts(exhaleImg);

    % Inhale
    inhaleImgNiftiPath = char(fullfile(outputDirectory, ...
        inhaleImgFileName + ".nii"));
    inhaleImgPreprocNiftiPath = char(fullfile(outputDirectory, ...
        inhaleImgFileName + "_preproc.nii"));
    inhalePointsMaskNiftiPath = char(fullfile(outputDirectory, ...
        inhaleImgFileName + '_keypoint_mask.nii'));
    inhaleLungMaskNiftiPath = char(fullfile(outputDirectory, ...
        inhaleImgFileName + '_mask.nii'));

    % Exhale
    exhaleImgNiftiPath = char(fullfile(outputDirectory, ...
        exhaleImgFileName + ".nii"));
    exhaleImgPreprocNiftiPath = char(fullfile(outputDirectory, ...
        exhaleImgFileName + "_preproc.nii"));
    exhalePointsMaskNiftiPath = char(fullfile(outputDirectory, ...
        exhaleImgFileName + '_keypoint_mask.nii'));
    exhaleLungMaskNiftiPath = char(fullfile(outputDirectory, ...
        exhaleImgFileName + '_mask.nii'));

    if PREPROCESS
        save_nii(inhaleImgNifti, inhaleImgNiftiPath);
        save_nii(inhaleImgPreprocNifti, inhaleImgPreprocNiftiPath);
        save_nii(inhalePointsMaskNifti, inhalePointsMaskNiftiPath);
        save_nii(inhaleLungMaskNifti, inhaleLungMaskNiftiPath);
            
        save_nii(exhaleImgNifti, exhaleImgNiftiPath);
        save_nii(exhaleImgPreprocNifti, exhaleImgPreprocNiftiPath);
        save_nii(exhalePointsMaskNifti, exhalePointsMaskNiftiPath);
        save_nii(exhaleLungMaskNifti, exhaleLungMaskNiftiPath);
    end

    saveTime = toc;
    disp("Preprocessing completed in " + saveTime + " seconds.");

% ---------------------------------------------------------------------
% REGISTRATION
% ---------------------------------------------------------------------
    caseOutputDirectory = fullfile(outputDirectory, parameter)
    if ~exist(caseOutputDirectory, 'dir')
        mkdir(caseOutputDirectory);
    end

    if REGISTER
        % Measure time for registration
        disp("Starting registration with elastix ([exhale] -> [inhale]): " + caseName);
        tic;  % Start timing
        elastixCommand = sprintf('"%s" -m "%s" -mMask "%s" -f "%s" -fMask "%s" -p "%s" -p "%s" -p "%s" -out "%s"', ...
                          elastixPath, ...
                          exhaleImgNiftiPath, exhaleLungMaskNiftiPath, ...
                          inhaleImgNiftiPath, inhaleLungMaskNiftiPath, ...
                          parameterFileAffine, parameterFileBSpline, parameterFileBSpline2, ...
                          caseOutputDirectory);

        [status, result] = system(elastixCommand);
    else
        status = 0;
    end

    registrationTime = toc;  % End timing

    if status == 0
        % Rename the result image
        resultImagePath = fullfile(caseOutputDirectory, "result." + parameter + ".nii");
        if REGISTER
            movefile(fullfile(caseOutputDirectory, "result.2.nii"), resultImagePath);
        end

        % Rename the TransformParameters.0.txt file
        originalTransformParametersPath = fullfile( ...
            caseOutputDirectory, "TransformParameters.2.txt");
        newTransformParametersPath = fullfile( ...
            caseOutputDirectory, "TransformParameters." + parameter + ".txt");
        if REGISTER
            movefile(originalTransformParametersPath, newTransformParametersPath);
        end

        disp("Registration with elastix completed in " + registrationTime + " seconds.");


% -----------------------------------------------------------------
% TRANSFORMIX
% -----------------------------------------------------------------
        tic; % start timer

% -----------------------------------------------------------------
% KEYPOINT TRANSFORM
% -----------------------------------------------------------------
        
        inPtsPath2 = fullfile(caseOutputDirectory, "inhalePts.txt");
        exPtsPath2 = fullfile(caseOutputDirectory, "exhalePts.txt");
        formatPointsFile(inhalePointsPath, inPtsPath2);
        formatPointsFile(exhalePointsPath, exPtsPath2);

        disp("Starting transformation with transformix: " + caseName);
        transformixCommand = sprintf('"%s" -def "%s" -tp "%s" -out "%s"', ...
                      transformixPath, inPtsPath2, ...
                      newTransformParametersPath, caseOutputDirectory);
        [statusKeypoint, resultKeypoint] = system(transformixCommand);

% -----------------------------------------------------------------
% LUNG MASK TRANSFORM
% -----------------------------------------------------------------
        % Prepare the inhale lung mask transformation
        inhaleLungMaskTransformParamPath = fullfile(caseOutputDirectory, "TransformParameters_mask." + parameter + ".txt");
        copyfile(newTransformParametersPath, inhaleLungMaskTransformParamPath);
        
        % Update the interpolator to NN for mask transformation
        updateTransformixParameterFile(inhaleLungMaskTransformParamPath, 'ResampleInterpolator', 'FinalNearestNeighborInterpolator');
        updateTransformixParameterFile(inhaleLungMaskTransformParamPath, 'FinalBSplineInterpolationOrder', '0');
        updateTransformixParameterFile(inhaleLungMaskTransformParamPath, 'WriteResultImage', 'true');

        
        % Transform the inhale lung mask
        disp("Starting transformation of inhale lung mask with transformix: " + caseName);
        transformixLungMaskCommand = sprintf('"%s" -in "%s" -tp "%s" -out "%s"', ...
                          transformixPath, inhaleLungMaskNiftiPath, ...
                          inhaleLungMaskTransformParamPath, caseOutputDirectory);
        [statusLungMask, resultLungMask] = system(transformixLungMaskCommand);
        
        if statusLungMask == 0
            transformedLungMaskPath = char(...
                fullfile(caseOutputDirectory, "transformed_inhale_lung_mask." + parameter + ".nii"));
            movefile(fullfile(caseOutputDirectory, "result.nii"), transformedLungMaskPath);
            disp("Inhale lung mask transformation completed: " + transformedLungMaskPath);
        else
            disp("Error during inhale lung mask transformation, check transformix log!");
        end


        transformationTime = toc;
        disp("Transformation with transformix completed in " ...
            + transformationTime + " seconds.");


% -----------------------------------------------------------------
% CALCULATE RESULTS
% -----------------------------------------------------------------
        % Masked registered image
        registeredImgNifti = load_nii(char(resultImagePath));
        transformedMaskNifti = load_nii(char(transformedLungMaskPath));
        registeredImgMaskedImg = registeredImgNifti.img .* transformedMaskNifti.img;
        
        % Masked ground truth image
        inhaleImgMaskedImg = inhaleImgNifti.img .* inhaleLungMaskNifti.img;

        % Compute Metrics
        % diceCoefficient = computeDSC(binarizedMovingMask, binarizedFixedMask);
        % ncc = computeNCC(inhaleImgMaskedImg, registeredImgMaskedImg);
        % ngc = computeNGC(inhaleImgMaskedImg, registeredImgMaskedImg);

        ncc = computeNCC( ...
            inhaleImgNifti.img .* inhaleLungMaskNifti.img, ...
            registeredImgNifti.img .* inhaleLungMaskNifti.img);
        ngc = computeNGC( ...
            inhaleImgNifti.img .* inhaleLungMaskNifti.img, ...
            registeredImgNifti.img .* inhaleLungMaskNifti.img);

        if statusKeypoint == 0
            % Rename the transformation result
            transformedResultPath = char(...
                fullfile(caseOutputDirectory, "result_transform."+ parameter + ".txt"));
            movefile(fullfile(caseOutputDirectory, "outputpoints.txt"), transformedResultPath);


            fixedCoordinates = exhalePoints;
            transformedCoordinates = extractOutputPoints(transformedResultPath);

            % Ensure the point sets have the same number of points
            if size(fixedCoordinates) ~= size(transformedCoordinates)
                error('The number of points in the two files do not match.');
            end

            % Calculate the Euclidean distances
            distances = fixedCoordinates - transformedCoordinates;
            distancesVoxel = sqrt(sum(distances.^2, 2));
            distancesMilimeter = sqrt(sum((distances.*imageSpacing).^2, 2));

            % Save the distances into a text file
            distanceFilePath = char( ...
                fullfile(caseOutputDirectory, 'euclidean_distances.txt'));
            writematrix(distancesVoxel, distanceFilePath);

            disp("Average distance: " + mean(distancesVoxel) + " voxels;")
            disp("Average distance: " + mean(distancesMilimeter) + " mm;")
            disp("NCC: " + ncc)
            disp("NGC: " + ngc)
        else
            disp("There has been an error, see transformix log!");
        end
    else
        disp("There has been an error, see elastix log!");
    end

% -----------------------------------------------------------------
% STORE RESULTS
% -----------------------------------------------------------------
    timingData.(caseName).readTime = readTime;
    timingData.(caseName).saveTime = saveTime;
    timingData.(caseName).registrationTime = registrationTime;
    timingData.(caseName).transformationTime = transformationTime;
    timingData.(caseName).total = readTime + saveTime + registrationTime + transformationTime;
    TREdata.(caseName).voxel = mean(distancesVoxel);
    TREdata.(caseName).mm = mean(distancesMilimeter);
    TREdata.(caseName).NCC = ncc;
    TREdata.(caseName).NGC = ngc;

% ---------------------------------------------------------------------
% DELETE UNNECESSARY RESULT FILES
% ---------------------------------------------------------------------
    % List of files to retain
    retainFiles = {
        resultImagePath, ...
        newTransformParametersPath, ...
        transformedResultPath, ...
        transformedLungMaskPath, ...
        distanceFilePath, ...
        inPtsPath2, ...
        exPtsPath2
    };

    % Add logs to retain list
    elastixLogPath = fullfile(caseOutputDirectory, "elastix.log");
    transformixLogPath = fullfile(caseOutputDirectory, "transformix.log");
    retainFiles = [retainFiles, elastixLogPath, transformixLogPath];
    
    % Get all files in the output directory
    allFiles = dir(caseOutputDirectory);
    allFilePaths = fullfile({allFiles.folder}, {allFiles.name});
    
    % Files to delete: those not in retainFiles
    filesToDelete = setdiff(allFilePaths, retainFiles);
    
    % Delete unwanted files
    for fileIdx = 1:length(filesToDelete)
        delete(filesToDelete{fileIdx});
    end

    disp("---------------------------------------------------------------")
end

% ---------------------------------------------------------------------
% SAVE RESULTS
% ---------------------------------------------------------------------
% Initialize tables to store the unpacked data
timingTable = table();
TRETable = table();

% Unpack timingData
caseNames = fieldnames(timingData);
for i = 1:numel(caseNames)
    caseName = caseNames{i};
    caseTimingData = timingData.(caseName);
    caseTimingData.caseName = caseName; % Add case name as a field
    timingTable = [timingTable; struct2table(caseTimingData)];
end

% Unpack TREdata
caseNames = fieldnames(TREdata);
for i = 1:numel(caseNames)
    caseName = caseNames{i};
    caseTREData = TREdata.(caseName);
    caseTREData.caseName = caseName; % Add case name as a field
    TRETable = [TRETable; struct2table(caseTREData)];
end

% Save the tables to CSV files
timingCSVPath = fullfile(dataPath, 'timing_data' + parameter + '.csv');
TRECSVPath = fullfile(dataPath, 'TRE_data' + parameter + '.csv');
writetable(timingTable, timingCSVPath);
writetable(TRETable, TRECSVPath);

% Display summary of times for all cases
disp("Timing data for all cases:");
disp(timingTable);
disp("TRE data for all cases:");
disp(TRETable);
