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