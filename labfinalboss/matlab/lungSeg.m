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

%% LungSeg Utility Functions
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