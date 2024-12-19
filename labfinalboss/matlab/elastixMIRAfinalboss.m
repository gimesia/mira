function nii = read_raw_to_nifti(filename, shape, spacing)

    % Read binary data
    fid = fopen(filename, 'rb');
    data = fread(fid, inf, 'int16','ieee-le'); % little-endian unsigned 16bit integer
    fclose(fid);

    data = reshape(data, shape);

    % Create NIfTI image
    nii = make_nii(data, spacing);

    nii.hdr.hist.sform_code = 1;  % Use a valid sform
    nii.hdr.hist.srow_x = [-spacing(1) 0 0 0];  % spacing along x-axis
    nii.hdr.hist.srow_y = [0 -spacing(2) 0 0];  % spacing along y-axis
    nii.hdr.hist.srow_z = [0 0 -spacing(3) 0];  % spacing along z-axis
end

function binarized_img_3d = binarize_lung_3d(img_3d, threshold_val)
    if ndims(img_3d) ~= 3
        error('Input image must be 3D');
    end

    binarized_img_3d = zeros(size(img_3d));
    for i = 1:size(img_3d, 3)
        img = img_3d(:,:,i);

        % Apply Gaussian blur
        blurred_img = imgaussfilt(img, 5);

        % Threshold the image
        binary_img = imbinarize(blurred_img, threshold_val);
        binary_img = imcomplement(binary_img);

        % Remove objects touching the image border
        binary_img = imclearborder(binary_img);

        % Perform binary closing
        closed_img = imclose(binary_img, strel('disk', 5));

        % Label connected components
        labeled_img = bwlabel(closed_img);
        props = regionprops(labeled_img, 'Area');

        % Filter out small objects
        idx = find([props.Area] > 1500); % Adjust the threshold as needed
        filtered_img = ismember(labeled_img, idx);

        % Remove non-central objects (replace with your specific implementation)
        % ...

        % Fill remaining holes
        filtered_img = imfill(filtered_img, 'holes');

        binarized_img_3d(:,:,i) = filtered_img;
    end
end

function mask = create_keypoints_mask(points, shape)
    mask = zeros(shape);
    for i = 1:size(points, 1)
        x = round(points(i, 1));
        y = round(points(i, 2));
        z = round(points(i, 3));
        if x > 0 && x <= shape(1) && y > 0 && y <= shape(2) && z > 0 && z <= shape(3)
            mask(x, y, z) = 1;
        end
    end
end

function maskNifti = create_keypoints_mask_nifti(points, templateNifti)
    % Extract image dimensions from the template NIfTI
    dims = templateNifti.hdr.dime.dim(2:4);  % Dimensions of the image   ű
    spacing = templateNifti.hdr.dime.pixdim(2:4);
    origin = templateNifti.hdr.hist.originator(2:4);

    % Create an empty mask with the same dimensions as the template
    mask = create_keypoints_mask(points, dims);

    % Dilate (optionally)
    se = strel("cube", 3);
    mask = imdilate(mask, se);

    maskNifti = make_nii(mask, spacing, origin);
    maskNifti.hdr.hist.sform_code = templateNifti.hdr.hist.sform_code;  % Use a valid sform
    maskNifti.hdr.hist.srow_x = templateNifti.hdr.hist.srow_x;
    maskNifti.hdr.hist.srow_y = templateNifti.hdr.hist.srow_y;
    maskNifti.hdr.hist.srow_z = templateNifti.hdr.hist.srow_z;
end

function nii = preprocess(niftiImg)
    img = int32(niftiImg.img);

    % img_filtered = imgaussfilt3(img, 0.5);
    img_filtered = medfilt3(img, [5 5 5]);

    mx = 2000; % max(img(:));
    clipped = rescale(img_filtered, 0, mx, 'InputMin', 0, 'InputMax', mx);

    normalizedImg = rescale(single(clipped), 0, 1, 'InputMin', 0, ...
        'InputMax', 1400);

    % Plot histograms for visualization before and after normalization
    % figure;
    % subplot(1, 3, 1);
    % histogram(clipped, 50);
    % title('Before Normalization');
    %
    % subplot(1, 3, 2);
    % histogram(normalizedImg(:), 50);
    % title('After Normalization');
    %
    % subplot(1, 3, 3);
    % histogram(histeq(normalizedImg(:)), 50);
    % title('After equalization');


    % Update the NIfTI structure with the normalized image
    nii = niftiImg;
    nii.img = (normalizedImg);

    % Automatically update the header fields related to intensity
    nii.hdr.dime.cal_min = min(normalizedImg(:));
    nii.hdr.dime.cal_max = max(normalizedImg(:));
end

function mask = segmentLungs3(img_3d)
    mask = zeros(size(img_3d));

    for i=1:size(mask, 3)
        blurred_img = imgaussfilt(img_3d(:,:,i), 3);

        % Threshold the image
        binary_img = imbinarize(blurred_img, 0.5);
        binary_img = imcomplement(binary_img);

        if ~any(binary_img(:))
            % Skip this iteration
            continue;
        end

        % Remove objects touching the image border
        binary_img = imclearborder(binary_img);

        % Perform binary closing
        closed_img = imclose(binary_img, strel('disk', 5));

        closed_img = imfill(closed_img, 'holes');

        mask(:,:,i) = closed_img;
    end

    for i=1:size(mask,1)
        slice = mask(i,:,:);

        if ~any(slice(:))
            % Skip this iteration
            continue;
        end

        reshaped = reshape(slice, [size(slice,2),size(slice,3)]);

        closed_img = imclose(reshaped, strel('disk', 5));

        filled = imfill(closed_img, 'holes');

        % imshow(filled)

        filled = reshape(filled,size(slice));
        mask(i,:,:) = filled;
    end

    for i=1:size(mask,2)
        slice = mask(:,i,:);

        if ~any(slice(:))
            % Skip this iteration
            continue;
        end

        reshaped = reshape(slice, [size(slice,1),size(slice,3)]);

        closed_img = imclose(reshaped, strel('disk', 5));

        filled = imfill(closed_img, 'holes');

        % imshow(filled)

        filled = reshape(filled,size(slice));
        mask(:,i,:) = filled;
    end

end

function formatPointsFile(inputFile, outputFile)
    % Read the data from the input file
    data = readmatrix(inputFile);

    % Ensure the file contains exactly 300 rows
    numRows = size(data, 1);
    if numRows ~= 300
        error('The input file must contain exactly 300 rows of data.');
    end

    % Open the output file for writing
    fid = fopen(outputFile, 'w');
    if fid == -1
        error('Could not open the output file for writing.');
    end

    % Write the required headers
    fprintf(fid, 'index\n');
    fprintf(fid, '%d\n', numRows);

    % Write the data rows
    for i = 1:numRows
        fprintf(fid, '%d\t%d\t%d\n', data(i, :));
    end

    % Close the file
    fclose(fid);

    % fprintf('Formatted file saved as: %s\n', outputFile);
end

function outputMatrix = extractOutputPoints(filename)
    % Open the file for reading
    fid = fopen(filename, 'r');
    if fid == -1
        error('Could not open the file: %s', filename);
    end

    % Initialize an empty array to hold the extracted points
    outputMatrix = [];

    % Read the file line by line
    while ~feof(fid)
        line = fgetl(fid);
        % Check if the line contains "OutputPoint"
        if contains(line, 'OutputIndexFixed =')
            % Extract the numbers after "OutputPoint = ["
            tokens = regexp(line, 'OutputIndexFixed = \[([-\d.e ]+)\]', 'tokens');
            if ~isempty(tokens)
                % Convert the string of numbers into a row vector
                point = str2num(tokens{1}{1}); %#ok<ST2NM>
                outputMatrix = [outputMatrix; point]; %#ok<AGROW>
            end
        end
    end

    % Close the file
    fclose(fid);

    % Ensure the output matrix has 300 rows
    if size(outputMatrix, 1) ~= 300
        error('The number of extracted points (%d) does not match the expected count (300).', size(outputMatrix, 1));
    end
end

metaData
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

% Path to parameter file
parameter = "Par0016"
parameterFileBSpline = ...
    "C:\Users\gimes\OneDrive\MAIA\3_UdG\classes\MIRA\lab\lab2\"+...
    "elastix_model_zoo\models\Par0016\Par0016.multibsplines.lung.sliding.txt";
REGISTER = 1

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

    save_nii(inhaleImgNifti, inhaleImgNiftiPath);
    save_nii(inhaleImgPreprocNifti, inhaleImgPreprocNiftiPath);
    save_nii(inhalePointsMaskNifti, inhalePointsMaskNiftiPath);
    save_nii(inhaleLungMaskNifti, inhaleLungMaskNiftiPath);

    % Exhale
    exhaleImgNiftiPath = char(fullfile(outputDirectory, ...
        exhaleImgFileName + ".nii"));
    exhaleImgPreprocNiftiPath = char(fullfile(outputDirectory, ...
        exhaleImgFileName + "_preproc.nii"));
    exhalePointsMaskNiftiPath = char(fullfile(outputDirectory, ...
        exhaleImgFileName + '_keypoint_mask.nii'));
    exhaleLungMaskNiftiPath = char(fullfile(outputDirectory, ...
        exhaleImgFileName + '_mask.nii'));

    save_nii(exhaleImgNifti, exhaleImgNiftiPath);
    save_nii(exhaleImgPreprocNifti, exhaleImgPreprocNiftiPath);
    save_nii(exhalePointsMaskNifti, exhalePointsMaskNiftiPath);
    save_nii(exhaleLungMaskNifti, exhaleLungMaskNiftiPath);

    saveTime = toc;
    disp("Preprocessing completed in " + saveTime + " seconds.");

    caseOutputDirectory = fullfile(outputDirectory, parameter);
    if ~exist(caseOutputDirectory, 'dir')
        mkdir(caseOutputDirectory);
    end

    if REGISTER

        % Measure time for registration
        disp("Starting registration with elastix ([exhale] -> [inhale]): " + caseName);
        tic;  % Start timing
        elastixCommand = sprintf('"%s" -m "%s" -mMask "%s" -f "%s" -fMask "%s" -p "%s" -out "%s"', ...
                          elastixPath, ...
                          exhaleImgNiftiPath, exhaleLungMaskNiftiPath, ...
                          inhaleImgNiftiPath, inhaleLungMaskNiftiPath, ...
                          parameterFileBSpline, caseOutputDirectory);

        [status, result] = system(elastixCommand);
    else
        status = 0;
    end

    registrationTime = toc;  % End timing

    if status == 0
        % Rename the result image
        resultImagePath = fullfile(caseOutputDirectory, "result." + parameter + ".nii");
        if REGISTER
            movefile(fullfile(caseOutputDirectory, "result.0.nii"), resultImagePath);
        end

        % Rename the TransformParameters.0.txt file
        originalTransformParametersPath = fullfile( ...
            caseOutputDirectory, "TransformParameters.0.txt");
        newTransformParametersPath = fullfile( ...
            caseOutputDirectory, "TransformParameters." + parameter + ".txt");
        if REGISTER
            movefile(originalTransformParametersPath, newTransformParametersPath);
        end


        disp("Registration with elastix completed in " + registrationTime + " seconds.");

        inPtsPath2 = fullfile(caseOutputDirectory, "inhalePts.txt");
        exPtsPath2 = fullfile(caseOutputDirectory, "exhalePts.txt");

        % transformParamsPath2 = fullfile(caseOutputDirectory, "TransformParameters2.test.txt");
        % updateTransformixParameterFile(newTransformParametersPath, transformParamsPath2)

        formatPointsFile(inhalePointsPath, inPtsPath2);
        formatPointsFile(exhalePointsPath, exPtsPath2);

        tic;
        % Use transformix to transform the inhalePoints
        disp("Starting transformation with transformix: " + caseName);
        transformixCommand = sprintf('"%s" -def "%s" -tp "%s" -out "%s"', ...
                      transformixPath, inPtsPath2, ...
                      newTransformParametersPath,caseOutputDirectory);
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
    TREdata.(caseName).voxel = mean(distancesVoxel);
    TREdata.(caseName).mm = mean(distancesMilimeter);
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
timingCSVPath = fullfile(dataPath, 'timing_dataPar0016.csv');
TRECSVPath = fullfile(dataPath, 'TRE_dataPar0016.csv');
writetable(timingTable, timingCSVPath);
writetable(TRETable, TRECSVPath);

% Display summary of times for all cases
disp("Timing data for all cases:");
disp(timingTable);
disp("TRE data for all cases:");
disp(TRETable);
