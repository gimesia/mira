function binarized_img_3d = segmentLungs1(img_3d, threshold_val)
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

        % Remove non-central objects
        % ...

        % Fill remaining holes
        filtered_img = imfill(filtered_img, 'holes');

        binarized_img_3d(:,:,i) = filtered_img;
    end
end