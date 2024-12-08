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

    % Save the NIfTI image
    dstPath = strrep(filename, '.img', '.nii');
    if ~ischar(dstPath) % && ~isstring(dstPath)
        dstPath = char(dstPath); % Convert to char if needed
    end
end

function mask = create_binary_mask(points, shape)
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

function maskNifti = create_binary_mask_nifti(points, templateNifti)
    % Extract image dimensions from the template NIfTI
    dims = templateNifti.hdr.dime.dim(2:4);  % Dimensions of the image

    % Create an empty mask with the same dimensions as the template
    mask = create_binary_mask(points, dims);

    % Dilate (optionally)
    se = strel("cube", 3);
    mask = imdilate(mask, se);

    maskNifti = make_nii(mask, [1,1,1], [1,1,1]);
end

function showMidSlice(niftiImg)
    midSlice = round(size(niftiImg.img,3)/2);
    imshow(niftiImg.img(:,:,midSlice));
end

function hist(V)
    figure
    histogram(V)
    xlim([min(V,[],"all") max(V,[],"all")])
    ylim([0 2e6])
    xlabel("Intensity")
    ylabel("Number of Voxels")
end

function nii = preprocess(niftiImg)
    nii = niftiImg;
    nii.img(nii.img<0) = 0;
    min_val = min(nii.img(:));
    max_val = max(nii.img(:));
    normalizedImg = (nii.img - min_val) / (max_val - min_val);

    hist(niftiImg.img);
    title("before");
    hist(nii.img)
    title("clipping");
    hist(normalizedImg)
    title("after");
    % Update the NIfTI structure with the normalized image
    nii.img = normalizedImg;

    % Update the header information
    nii.hdr.dime.cal_min = 0;
    nii.hdr.dime.cal_max = 1;
end

metaData
% Paths
elastixBasePath = "C:\Users\gimes\OneDrive\MAIA\3_UdG\classes\MIRA\lab\lab2\";
basePath = "C:\Users\gimes\OneDrive\MAIA\3_UdG\classes\MIRA\lab\labfinalboss\";
dataPath = "C:\Users\gimes\OneDrive\MAIA\3_UdG\classes\MIRA\lab\labfinalboss\data\";

rawImgFiles = dir(fullfile(dataPath, "**", "*BHCT.img")); % only raw
rawImgFiles = {rawImgFiles.name};
rawImgFiles = strcat(strtok(rawImgFiles, '_'), '\', rawImgFiles);

imgFiles = dir(fullfile(dataPath, "**", "*.nii.gz")); % only nifti
imgFiles = {imgFiles.name};  % Get all filenames in a cell array
imgFiles = strcat(strtok(imgFiles, '_'), '\', imgFiles);

pointFiles = dir(fullfile(dataPath, "**", "*xyz_r1.txt")); % only points txt
pointFiles= {pointFiles.name};  % Get all filenames in a cell array
pointFiles = strcat(strtok(pointFiles, '_'), '\', pointFiles);

cases = unique(strtok(pointFiles, "\"));

% Path to B-spline parameter file
parameter = "Par0016"
parameterFileBSpline = "C:\Users\gimes\OneDrive\MAIA\3_UdG\classes\MIRA\lab\lab2\elastix_model_zoo\models\Par0016\Par0016.multibsplines.lung.sliding.txt"

% Paths to elastix and transformix executables
elastixPath = elastixBasePath + "elastix-5.0.0-win64\elastix.exe";
transformixPath = elastixBasePath + "elastix-5.0.0-win64\transformix.exe";

% To track timing for each case
registrationTimes = struct();

for caseID = 1:1%length(cases)
    caseName = string(cases{caseID});
    disp("Processing Case [" + caseName + "]");
    idxs = cellfun(@(x) contains(x, caseName), rawImgFiles);
    idxs = find(idxs == 1);
    caseImgs = rawImgFiles(idxs);
    casePoints = pointFiles(idxs);
    index = find(strcmp(data.Label, caseName));
    imageDims = data.ImageDims{index};
    imageSpacing = data.Spacing{index};

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
    inhalePoints = fullfile(dataPath, casePoints{inhalePointsId == 1});
    inhalePointsData = readmatrix(inhalePoints);

    % Moving Points
    exhalePointsId = cellfun(@(x) contains(x, "eBH"), casePoints);
    exhalePoints = fullfile(dataPath, casePoints{exhalePointsId == 1});
    exhalePointsData = readmatrix(exhalePoints);

    outputDirectory = basePath + "data\" + caseName + "\out";
    if ~exist(outputDirectory, 'dir')
        mkdir(outputDirectory);
    end

    % Create lung masks using the NIfTI template
    % inhaleLungMaskNifti = uNetLungMask(inhaleImgNifti);
    % exhaleLungMaskNifti = uNetLungMask(exhaleImgNifti);




    % Create point masks using the NIfTI template
    inhalePointsMaskNifti = create_binary_mask_nifti(inhalePointsData, inhaleImgNifti);
    exhalePointsMaskNifti = create_binary_mask_nifti(exhalePointsData, exhaleImgNifti);



    % Create preprocessed NIfTI images
    inhaleImgPreprocNifti = preprocess(inhaleImgNifti);
    exhaleImgPreprocNifti = preprocess(exhaleImgNifti);

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
    % inhaleLungMaskNiftiPath = char(fullfile(outputDirectory, ...
    %     inhaleImgFileName + '_mask.nii'));

    save_nii(inhaleImgNifti, inhaleImgNiftiPath);
    save_nii(inhaleImgPreprocNifti, inhaleImgPreprocNiftiPath);
    save_nii(inhalePointsMaskNifti, inhalePointsMaskNiftiPath);
    % save_nii(inhaleLungMaskNifti, inhaleLungMaskNiftiPath);

    % Exhale
    exhaleImgNiftiPath = char(fullfile(outputDirectory, ...
        exhaleImgFileName + ".nii"));
    exhaleImgPreprocNiftiPath = char(fullfile(outputDirectory, ...
        exhaleImgFileName + "_preproc.nii"));
    exhalePointsMaskNiftiPath = char(fullfile(outputDirectory, ...
        exhaleImgFileName + '_keypoint_mask.nii'));
    % exhaleLungMaskNiftiPath = char(fullfile(outputDirectory, ...
         % exhaleImgFileName + '_mask.nii'));

    save_nii(exhaleImgNifti, exhaleImgNiftiPath);
    save_nii(exhaleImgPreprocNifti, exhaleImgPreprocNiftiPath);
    save_nii(exhalePointsMaskNifti, exhalePointsMaskNiftiPath);
    % save_nii(exhaleLungMaskNifti, exhaleLungMaskNiftiPath);

    caseOutputDirectory = fullfile(outputDirectory, parameter);
    if ~exist(caseOutputDirectory, 'dir')
        mkdir(caseOutputDirectory);
    end

    % Measure time for registration
    disp("Starting registration with elastix ([exhale] -> [inhale]): " + caseName);
    tic;  % Start timing
    elastixCommand = sprintf('"%s" -f "%s" -m "%s" -p "%s"  -labels "%s" -out "%s"', ...
                      elastixPath, exhaleImgPreprocNiftiPath, inhaleImgPreprocNiftiPath, ...
                      parameterFileBSpline, ...
                      inhalePointsMaskNiftiPath, caseOutputDirectory);
    [status, result] = system(elastixCommand);
    elapsedTime = toc;  % End timing

    % Store the time
    registrationTimes.(caseName) = elapsedTime;

    if status == 0
        disp("Registration with elastix completed in " + elapsedTime + " seconds.");
    else
        disp("There has been an error, see log!");
    end
end

% Display summary of registration times
disp("Registration times for all cases:");
disp(registrationTimes);