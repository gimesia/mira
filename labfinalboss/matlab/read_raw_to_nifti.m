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