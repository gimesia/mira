
metaData

% Path to parameter file
parameter = "Par0008"

parameterFileAffine = ...
    "C:\Users\gimes\OneDrive\MAIA\3_UdG\classes\MIRA\lab\lab2\" + ...
    "elastix_model_zoo\models\Par0008\Parameters.Par0008.affine.txt"
parameterFileBSpline = ...
    "C:\Users\gimes\OneDrive\MAIA\3_UdG\classes\MIRA\lab\lab2\" + ...
    "elastix_model_zoo\models\Par0008\Parameters.Par0008.elastic.txt"

REGISTER = 1
PREPROCESS = 0

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
        % Create preprocessed NIfTI images
        inhaleImgPreprocNifti = preprocess(inhaleImgNifti);
        exhaleImgPreprocNifti = preprocess(exhaleImgNifti);
    
        % Create lung masks using the NIfTI template
        inhaleLungMaskNifti = inhaleImgPreprocNifti;
        inhaleLungMaskNifti.img = segmentLungs3(inhaleLungMaskNifti.img);
        exhaleLungMaskNifti = exhaleImgPreprocNifti;
        exhaleLungMaskNifti.img = segmentLungs3(exhaleLungMaskNifti.img);
    
        % Create point masks using the NIfTI template
        inhalePointsMaskNifti = create_keypoints_mask_nifti(inhalePoints, inhaleImgNifti);
        exhalePointsMaskNifti = create_keypoints_mask_nifti(exhalePoints, exhaleImgNifti);
    end
    
    % Saving
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

    caseOutputDirectory = fullfile(outputDirectory, parameter)
    if ~exist(caseOutputDirectory, 'dir')
        mkdir(caseOutputDirectory);
    end

    if REGISTER
        % Measure time for registration
        disp("Starting registration with elastix ([exhale] -> [inhale]): " + caseName);
        tic;  % Start timing
        elastixCommand = sprintf('"%s" -m "%s" -mMask "%s" -f "%s" -fMask "%s" -p "%s" -p "%s" -out "%s"', ...
                          elastixPath, ...
                          exhaleImgNiftiPath, exhaleLungMaskNiftiPath, ...
                          inhaleImgNiftiPath, inhaleLungMaskNiftiPath, ...
                          parameterFileAffine, parameterFileBSpline, ...
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
            movefile(fullfile(caseOutputDirectory, "result.1.nii"), resultImagePath);
        end

        registeredImg = load_nii(char(resultImagePath));
        % Compute Metrics
        % diceCoefficient = computeDSC(binarizedMovingMask, binarizedFixedMask);
        ncc = computeNCC(inhaleImgNifti.img, registeredImg.img);
        ngc = computeNGC(inhaleImgNifti.img, registeredImg.img);


        % Rename the TransformParameters.0.txt file
        originalTransformParametersPath = fullfile( ...
            caseOutputDirectory, "TransformParameters.1.txt");
        newTransformParametersPath = fullfile( ...
            caseOutputDirectory, "TransformParameters." + parameter + ".txt");
        if REGISTER
            movefile(originalTransformParametersPath, newTransformParametersPath);
        end

        disp("Registration with elastix completed in " + registrationTime + " seconds.");

        inPtsPath2 = fullfile(caseOutputDirectory, "inhalePts.txt");
        exPtsPath2 = fullfile(caseOutputDirectory, "exhalePts.txt");

        formatPointsFile(inhalePointsPath, inPtsPath2);
        formatPointsFile(exhalePointsPath, exPtsPath2);

        tic;
        % Use transformix to transform the inhalePoints
        disp("Starting transformation with transformix: " + caseName);
        transformixCommand = sprintf('"%s" -def "%s" -tp "%s" -out "%s"', ...
                      transformixPath, inPtsPath2, ...
                      newTransformParametersPath, caseOutputDirectory);
        [status1, result1] = system(transformixCommand);

                transformixCommand = sprintf('"%s" -def "%s" -tp "%s" -out "%s"', ...
                      transformixPath, inPtsPath2, ...
                      newTransformParametersPath, caseOutputDirectory);
        [status1, result1] = system(transformixCommand);

        transformationTime = toc;
        disp("Transformation with transformix completed in " ...
            + transformationTime + " seconds.");

        if status1 == 0
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
        else
            disp("There has been an error, see transformix log!");
        end
    else
        disp("There has been an error, see elastix log!");
    end

    % Store the times
    timingData.(caseName).readTime = readTime;
    timingData.(caseName).saveTime = saveTime;
    timingData.(caseName).registrationTime = registrationTime;
    timingData.(caseName).transformationTime = transformationTime;
    timingData.(caseName).total = readTime + saveTime + registrationTime + transformationTime;
    TREdata.(caseName).voxel = mean(distancesVoxel);
    TREdata.(caseName).mm = mean(distancesMilimeter);
    TREdata.(caseName).NCC = ncc;
    TREdata.(caseName).NGC = ngc;
    disp("---------------------------------------------------------------")
end


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
