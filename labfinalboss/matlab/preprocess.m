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
