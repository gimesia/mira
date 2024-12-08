function nii = read_raw_to_nifti(filename, shape, spacing)
    % Read binary data
    fid = fopen(filename, 'rb');
    data = fread(fid, inf, 'int16','ieee-le'); % little-endian unsigned 16bit integer
    fclose(fid);

    size(data);
    shape(1)*shape(2)*shape(3);

    data = reshape(data, shape);

    % Create NIfTI image
    nii = make_nii(data, spacing, [0,0,0], 4);

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
    se = strel("cube", 5);
    %mask = imdilate(mask, se);

    maskNifti = make_nii(mask, [1,1,1], [1,1,1]);
end


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


for caseID = 1:1%length(cases)
    caseName = string(cases{caseID});
    disp("Processing Case [" + caseName + "]");
    idxs = cellfun(@(x) contains(x, caseName), rawImgFiles);
    idxs = find(idxs == 1);
    caseImgs = rawImgFiles(idxs);
    casePoints = pointFiles(idxs);
    caseMetaData = metaData(metaData.Label == caseName, :);
    imageDims = caseMetaData.ImageDims;
    imageSpacing = caseMetaData.Spacing;

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

    % Create binary masks
    % Create binary masks using the NIfTI template
    inhalePointsMaskNifti = create_binary_mask_nifti(inhalePointsData, inhaleImgNifti);
    exhalePointsMaskNifti = create_binary_mask_nifti(exhalePointsData, exhaleImgNifti);

    % Saving
    [~, inhaleImgFileName, ~] = fileparts(inhaleImg);
    [~, exhaleImgFileName, ~] = fileparts(exhaleImg);


    % Inhale
    inhaleImgNiftiPath = char(fullfile(outputDirectory, ...
        inhaleImgFileName + ".nii"))
    inhalePointsMaskNiftiPath = char(fullfile(outputDirectory, ...
        inhaleImgFileName + '_keypoint_mask.nii'))
    save_nii(inhaleImgNifti, inhaleImgNiftiPath);
    save_nii(inhalePointsMaskNifti, inhalePointsMaskNiftiPath);


    % Exhale
    exhaleImgNiftiPath = char(fullfile(outputDirectory, ...
        exhaleImgFileName + ".nii"))
    exhalePointsMaskNiftiPath = char(fullfile(outputDirectory, ...
        exhaleImgFileName + '_keypoint_mask.nii'))
    save_nii(exhaleImgNifti, exhaleImgNiftiPath);
    save_nii(exhalePointsMaskNifti, exhalePointsMaskNiftiPath);

    disp("Registration with elastix ([exhale] -> [inhale]): " + caseName);
    elastixCommand = sprintf('"%s" -f "%s" -m "%s" -p "%s" -labels "%s" -out "%s"', ...
                      elastixPath, exhaleImgNiftiPath, inhaleImgNiftiPath, ...
                      parameterFileBSpline, inhalePointsMaskNiftiPath, ...
                      outputDirectory);
    [status, result] = system(elastixCommand);
    if status == 0
        disp("Registration with elastix completed.")
    else
        disp("There has been an error, see log!")
    end
end
