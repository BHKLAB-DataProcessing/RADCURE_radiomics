# `helper_functions.py`

`loadDicomSITK`: This function reads a DICOM series and returns a SimpleITK Image.

`loadNRRDSITK`: This reads an NRRD file and returns it as a SimpleITK Image.

`loadNiftiSITK`: Similar to the above, but for NIFTI files.

`loadRTSTRUCTSITK`: Loads RTSTRUCT into a SimpleITK Image.

`saveImageAsNRRD`: Saves a SimpleITK Image as an NRRD file.

`loadSegmentation`: Loads segmentation from a binary file.

`makeBinROIMask`: Generates a binary mask from an ROI label.

`flattenImage`: Flattens a 3D image into 2D arrays.

`alignImages`: Aligns two images in terms of origin, direction, and spacing.

`displayImageSlice`: Displays a slice of an image.

`padSEGtoMatchCT`: Pads a segmentation image to match the dimensions of a CT image.

`getAcquisitionFileList`: Retrieves a list of files from an image acquisition.

`getROIVoxelLabel`: Retrieves the label value that identifies segmentation voxels.

`shuffleImage`: Shuffles all pixel values in an image.
