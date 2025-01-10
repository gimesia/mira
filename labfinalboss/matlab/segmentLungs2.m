function binaryImg = segmentLungs2(img_3d)
    downscaled = preprocessLungCT(img_3d);
    seg = lungSeg(downscaled, batchSize=8);
    binaryImg = postprocessBinaryLung(seg);
    % rescaled = rescale(binaryImg, 0, 2000);
    % binaryImg = activecontour(rescaled, binaryImg);
end

function downscaled = preprocessLungCT(img_3d)    
    % Orientation correction for pretrained UNet
    rotatedImg = imrotate3(single(img_3d), 90, [0 0 1]);

    % Downscaling for faster processing
    downscaled = zeros([256 256 size(rotatedImg, 3)]);
    for i=1:size(rotatedImg,3)
        % Resize the first two dimensions of V to size imSize
        downscaled(:,:,i) = imresize(rotatedImg(:,:,i), [256 256]);
    end
end

function upscaled = postprocessBinaryLung(img_3d)
    % Upscaling for faster processing
    upscaled = zeros([512 512 size(img_3d, 3)]);
    for i=1:size(img_3d,3)
        % Resize the first two dimensions
        resized = imresize(img_3d(:,:,i), [512 512]);
        slice = resized;
        slice(resized>0) = 1;
        slice = imfill(slice, 'holes');
        upscaled(:,:,i) = slice;
    end

    % Orientation correction for pretrained UNet
    upscaled = imrotate3(upscaled,-90,[0 0 1]);
    upscaled(upscaled ~= 0) = 1;

    % upscaled = imdilate(upscaled, strel('sphere',5));
end

function labelOut = lungSeg(inputVol,NameValueArgs)
% Perform semantic segementation of LungCT
%
% Input:
%   inputVol: 2-D slice or 3-D volume with x and y
%           dimensions 256x256 of single class (normalized values)
%
% Output:
%   labelOut: 2-D slice or 3-D volume (same size as input) of labels
%     Note: For Left-Right segmentation, labels with range from 0 (background) to 2

    arguments
        inputVol single 
        NameValueArgs.batchSize = 16;
        NameValueArgs.modeFilt = [9 9];
    end
    
    params = loadLungCTParams();
    labelOut = zeros(size(inputVol));

    imHeight = size(inputVol,1);
    imWidth = size(inputVol,2);
    numSlices = size(inputVol,3);

    % Reshape the Volume to SSCB
    inputVol = reshape(inputVol,[imHeight, imWidth, 1, numSlices]);
    inputVol = dlarray(inputVol,"SSCB"); % For batch processing (speed improvement)
    
    
    % Leverage GPU for speed if available
    if canUseGPU
        inputVol = gpuArray(inputVol);
    end
    
    for i = 1:NameValueArgs.batchSize:numSlices
        iSlice = sliceArray(i,NameValueArgs.batchSize,numSlices);
        outputVol = lungmaskFcn_R231(inputVol(:,:,:,iSlice),params.R231);
        labelOut(:,:,iSlice) = postprocessLungCT(outputVol,NameValueArgs.modeFilt);
    end

end

function params = loadLungCTParams()
    load('uNetLungmaskParams_R231.mat','params');
    parameters.R231=params;
    params = parameters;
end

function iSlice = sliceArray(i,batchSize,nSlices)
    iSlice = i:i+batchSize-1;
    if iSlice(end)>nSlices
        iSlice=i:nSlices;
    end
end

function labelOut = postprocessLungCT(outputVol,modeFilt)
% Postprocess to create labels matrix:
% Parse the labels from the initial label output to a single label matrix
% (2D for slice, 3D for volume). Also, a modefilt is used to smooth out the
% edges of the labels.
    outputVol = extractdata(outputVol);
    if canUseGPU
        outputVol = gather(outputVol);
    end
    [~, labelOutPrefilt] = max(outputVol,[],3);
    labelOutSqueezed = squeeze(labelOutPrefilt);
    for idx = 1:size(labelOutSqueezed,3)
        labelOut(:,:,idx) = modefilt(labelOutSqueezed(:,:,idx),modeFilt)-1;
    end
end
