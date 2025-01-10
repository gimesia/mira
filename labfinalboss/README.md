# Lab Final Boss

This project contains scripts and functions for processing and analyzing lung CT images. The main functionalities include preprocessing, segmentation, registration, and transformation of lung CT images.

- `matlab/`
  - `elastixMIRAfinalboss<ParameterName>.m`: MATLAB scripts for preprocessing, registration, and transformation using elastix and transformix.
  - `segmentLungs<N>.m`: Functions for lung segmentation.
- `data/`: Directory containing input and output data in raw and NIfTI format.

## Requirements

- MATLAB
- elastix
- transformix
- Pretrained UNet model parameters (`uNetLungmaskParams_R231.mat`)

## Usage

1. Copy the contents of the `matlab/` directory to your MATLAB working directory.
2. Open the script file for the desired parameter.
3. Run the script in MATLAB to process the lung CT images in the given directory.
4. (Optional) Modify the script to change the input and output directory or other parameters.

## Pipeline
### Preprocessing and Registration

1. Preprocesses the inhale and exhale CT images.
2. Saves the preprocessed images and masks.
3. Registers the exhale image to the inhale image using elastix.
4. Transforms the inhale points using transformix.
5. Calculates and saves the Euclidean distances between the fixed and transformed points.

The intermediate results are saved in the configured output directory.

### Data

The project expects the input data to be in raw format. The output data, including preprocessed images, masks, and registration results, are saved in NIfTI format.

### Output

The output of the scripts includes:
- Preprocessed images and masks.
- Registration results and transformation parameters.
- Euclidean distances between fixed and transformed points.
- Timing data and TRE (Target Registration Error) data.
