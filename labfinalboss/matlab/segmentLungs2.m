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

