#%%
print(walrus := "I am walrus")
walrus

#%%

import numpy as np
import nibabel as nib
import os
import time
from skimage.morphology import ball, binary_erosion, binary_dilation
from skimage.measure import label, regionprops
from scipy.ndimage import binary_fill_holes, gaussian_filter
import subprocess

import pandas as pd

data = {
    "Label": ["copd1", "copd2", "copd3", "copd4", "copd5", "copd6", "copd7", "copd8", "copd9", "copd10"],
    "ImageDims": [(512, 512, 121), (512, 512, 102), (512, 512, 126), (512, 512, 126), (512, 512, 131),
                   (512, 512, 119), (512, 512, 112), (512, 512, 115), (512, 512, 116), (512, 512, 135)],
    "Spacing": [(0.625, 0.625, 2.5), (0.645, 0.645, 2.5), (0.652, 0.652, 2.5), (0.59, 0.59, 2.5),
                    (0.647, 0.647, 2.5), (0.633, 0.633, 2.5), (0.625, 0.625, 2.5), (0.586, 0.586, 2.5),
                    (0.664, 0.664, 2.5), (0.742, 0.742, 2.5)],
    "#Features": [773, 612, 1172, 786, 1029, 633, 575, 791, 447, 480],
    "Displacement(mm)": ["25.90 (11.57)", "21.77 (6.46)", "12.29 (6.39)", "30.90 (13.49)", "30.90 (14.05)",
                          "28.32 (9.20)", "21.66 (7.66)", "25.57 (13.61)", "14.84 (10.01)", "22.48 (10.64)"],
    "#Repeats": ["150/3"] * 10,
    "Observers(mm)": ["0.65 (0.73)", "1.06 (1.51)", "0.58 (0.87)", "0.71 (0.96)", "0.65 (0.87)",
                       "1.06 (2.38)", "0.65 (0.78)", "0.96 (3.07)", "1.04 (2.54)", "0.87 (1.65)"]
}
metaData = pd.DataFrame(data)
metaData.set_index("Label", inplace=True)
print(metaData)

def read_raw_to_nifti(filename, shape, spacing):
    # Read binary data
    data = np.fromfile(filename, dtype=np.int16)
    data = data.reshape(shape)

    # Create NIfTI image
    affine = np.diag([-spacing[0], -spacing[1], -spacing[2], 1])
    nii = nib.Nifti1Image(data, affine)

    # Save the NIfTI image
    dst_path = filename.replace('.img', '.nii')
    nib.save(nii, dst_path)
    return nii


def create_binary_mask(points, shape):
    mask = np.zeros(shape, dtype=np.uint8)
    for point in points:
        x, y, z = map(int, np.round(point))
        if 0 <= x < shape[0] and 0 <= y < shape[1] and 0 <= z < shape[2]:
            mask[x, y, z] = 1
    return mask


def create_binary_mask_nifti(points, template_nifti):
    dims = template_nifti.shape
    mask = create_binary_mask(points, dims)
    mask_nifti = nib.Nifti1Image(mask, template_nifti.affine)
    return mask_nifti


def segment_lungs(image):
    return image


# Paths
elastix_base_path = fr"C:\Users\gimes\OneDrive\MAIA\3_UdG\classes\MIRA\lab\lab2"
base_path = fr"C:\Users\gimes\OneDrive\MAIA\3_UdG\classes\MIRA\lab\labfinalboss"
data_path = fr"C:\Users\gimes\OneDrive\MAIA\3_UdG\classes\MIRA\lab\labfinalboss\data"

raw_img_files = [os.path.join(dp, f) for dp, dn, filenames in os.walk(data_path) for f in filenames if
                 f.endswith('BHCT.img')]
img_files = [os.path.join(dp, f) for dp, dn, filenames in os.walk(data_path) for f in filenames if
             f.endswith('.nii.gz')]
point_files = [os.path.join(dp, f) for dp, dn, filenames in os.walk(data_path) for f in filenames if f.endswith('.txt')]

cases = list(set(os.path.basename(f).split('_')[0] for f in point_files if ("Info" not in f) and ("Transform" not in f)) )
case = sorted(cases)
# Path to B-spline parameter file
parameter_file_bspline = os.path.join(elastix_base_path, "elastix_model_zoo", "models", "Par0016",
                                      "Par0016.multibsplines.lung.sliding.txt")

# Paths to elastix and transformix executables
elastix_path = os.path.join(elastix_base_path, "elastix-5.0.0-win64", "elastix.exe")
transformix_path = os.path.join(elastix_base_path, "elastix-5.0.0-win64", "transformix.exe")

# To track timing for each case
registration_times = {}

for case_name in cases:
    print(f"Processing Case [{case_name}]")
    case_imgs = [f for f in raw_img_files if case_name in f]
    case_points = [f for f in point_files if case_name in f]

    print(cases)
    # Assuming metaData is a dictionary or similar structure
    case_meta_data = metaData.loc[case_name]
    image_dims = case_meta_data['ImageDims']
    image_spacing = case_meta_data['Spacing']

    # Moving Img
    inhale_img = next(f for f in case_imgs if "iBH" in f)
    inhale_img_nifti = read_raw_to_nifti(inhale_img, image_dims, image_spacing)
    print(inhale_img_nifti)

    # Fixed Img
    exhale_img = next(f for f in case_imgs if "eBH" in f)
    exhale_img_nifti = read_raw_to_nifti(exhale_img, image_dims, image_spacing)

    # Fixed Points
    inhale_points = next(f for f in case_points if "iBH" in f)
    inhale_points_data = np.loadtxt(inhale_points)

    # Moving Points
    exhale_points = next(f for f in case_points if "eBH" in f)
    exhale_points_data = np.loadtxt(exhale_points)

    output_directory = os.path.join(base_path, "data", case_name, "out")
    os.makedirs(output_directory, exist_ok=True)

    # Create lung masks using the NIfTI template
    inhale_lung_mask = segment_lungs(inhale_img_nifti.get_fdata())
    exhale_lung_mask = segment_lungs(exhale_img_nifti.get_fdata())
    inhale_lung_mask_nifti = nib.Nifti1Image(inhale_lung_mask, inhale_img_nifti.affine)
    exhale_lung_mask_nifti = nib.Nifti1Image(exhale_lung_mask, exhale_img_nifti.affine)

    # Create point masks using the NIfTI template
    inhale_points_mask_nifti = create_binary_mask_nifti(inhale_points_data, inhale_img_nifti)
    exhale_points_mask_nifti = create_binary_mask_nifti(exhale_points_data, exhale_img_nifti)

    # Saving
    inhale_img_file_name = os.path.splitext(os.path.basename(inhale_img))[0]
    exhale_img_file_name = os.path.splitext(os.path.basename(exhale_img))[0]

    # Inhale
    inhale_img_nifti_path = os.path.join(output_directory, f"{inhale_img_file_name}.nii")
    inhale_points_mask_nifti_path = os.path.join(output_directory, f"{inhale_img_file_name}_keypoint_mask.nii")
    inhale_lung_mask_nifti_path = os.path.join(output_directory, f"{inhale_img_file_name}_mask.nii")

    nib.save(inhale_img_nifti, inhale_img_nifti_path)
    nib.save(inhale_points_mask_nifti, inhale_points_mask_nifti_path)
    nib.save(inhale_lung_mask_nifti, inhale_lung_mask_nifti_path)

    # Exhale
    exhale_img_nifti_path = os.path.join(output_directory, f"{exhale_img_file_name}.nii")
    exhale_points_mask_nifti_path = os.path.join(output_directory, f"{exhale_img_file_name}_keypoint_mask.nii")
    exhale_lung_mask_nifti_path = os.path.join(output_directory, f"{exhale_img_file_name}_mask.nii")

    nib.save(exhale_img_nifti, exhale_img_nifti_path)
    nib.save(exhale_points_mask_nifti, exhale_points_mask_nifti_path)
    nib.save(exhale_lung_mask_nifti, exhale_lung_mask_nifti_path)

    # Measure time for registration
    print(f"Starting registration with elastix ([exhale] -> [inhale]): {case_name}")
    start_time = time.time()
    elastix_command = [
        elastix_path,
        '-f', exhale_img_nifti_path,
        '-m', inhale_img_nifti_path,
        '-p', parameter_file_bspline,
        '-labels', inhale_points_mask_nifti_path,
        '-out', output_directory
    ]
    result = subprocess.run(elastix_command, capture_output=True, text=True)
    elapsed_time = time.time() - start_time

    # Store the time
    registration_times[case_name] = elapsed_time

    if result.returncode == 0:
        print(f"Registration with elastix completed in {elapsed_time} seconds.")
    else:
        print("There has been an error, see log!")

# Display summary of registration times
print("Registration times for all cases:")
print(registration_times)