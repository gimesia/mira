
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
