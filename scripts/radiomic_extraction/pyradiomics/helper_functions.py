import pandas as pd
from radiomics import featureextractor, getFeatureClasses
import radiomics
import os
import SimpleITK as sitk
import numpy as np
import matplotlib.pyplot as plt
from imgtools.ops import StructureSetToSegmentation
from imgtools.io import read_dicom_auto
from dicom_parser import Series
import pydicom

def loadDicomSITK(imgPath: str):
    """Read DICOM series as SimpleITK Image.

    Parameters
    ----------
    img_path
        Path to directory containing the DICOM series to load.

    Returns
    -------
    The loaded image. 

    """
    # Set up the reader for the DICOM series
    reader = sitk.ImageSeriesReader()
    dicomNames = reader.GetGDCMSeriesFileNames(imgPath)
    reader.SetFileNames(dicomNames)
    return reader.Execute()


def loadNRRDSITK(imgPath: str):
    """Read in NRRD as SimpleITK Image.

    Parameters
    ----------
    img_path
        Path to the NRRD file to load.

    Returns
    -------
    The loaded image. 

    """
    reader = sitk.ImageFileReader()
    dicomNames = reader.SetImageIO("NrrdImageIO")
    reader.SetFileName(imgPath)
    return reader.Execute()


def loadNiftiSITK(imgPath: str):
    """Read in NIFTI as SimpleITK Image.

    Parameters
    ----------
    img_path
        Path to the NIFTI file to load.

    Returns
    -------
    The loaded image. 

    """
    reader = sitk.ImageFileReader()
    reader.SetImageIO("NiftiImageIO")
    reader.SetFileName(imgPath)
    return reader.Execute()


def loadRTSTRUCTSITK(rtstructPath: str,
                     originalImageDirPath: str,
                     roiNames:str = None):
    """ Load RTSTRUCT into SimpleITK Image.

    Parameters
    ----------
    rtstructPath
        Path to the DICOM file containing the RTSTRUCT
    originalImageDirPath
        Path to the directory containing the DICOMS for the original image the segmentation 
        was created from. This is required to load the RTSTRUCT.
    roiNames
        Identifier for which region(s) of interest to load from the total segmentation file
    
    Returns
    -------
    The loaded RTSTRUCT image as a SimpleITK image object.
    The segmentation label is set to 1.
    """

    # Set up segmentation loader
    makeMask = StructureSetToSegmentation(roi_names=roiNames)
    
    # Read in the CT and segmentation DICOMs into SITK Images
    ctImage = read_dicom_auto(originalImageDirPath)
    segImage = read_dicom_auto(rtstructPath)

    try:
        # Get the individual ROI masks
        segMasks = makeMask(segImage, ctImage.image, existing_roi_indices={"background":0}, ignore_missing_regex=False)
    except ValueError:
        return {}
        
    # Get list of ROIs present in this rtstruct
    roiNames = segMasks.raw_roi_names
    # Initialize dictionary to store ROI names and images
    roiStructs = {}
    # Get each roi and its label and store in dictionary
    for roi in roiNames:
        # Get the mask for this ROI
        roiMask = segMasks.get_label(name=roi)
        # Store the ROI name and image
        roiStructs[roi] = roiMask
    
    return roiStructs


def saveImageAsNRRD(image: sitk.Image,
                fileName: str,
                outputDirPath:str):
    """ Save out sitk.Image as .nrrd

    Parameters
    ----------
    image
        sitk Image to save out
    fileName
        Name to save file as
    outputDirPath
        Path to directory to save image to

    Returns
    -------
    Full path to the saved image file
    """
    fullFilename = os.path.join(outputDirPath, (fileName + ".nrrd"))

    writer = sitk.ImageFileWriter()
    writer.SetFileName(fullFilename)
    writer.Execute(image)

    return fullFilename


def loadSegmentation(imgPath: str,
                     modality: str,
                     originalImageDirPath: str = None,
                     roiNames: str = None):
    ''' Function to load a segmentation with the correct function.
    
    Parameters
    ----------
    imgPath
        Path to the segmentation file to load
    modality
        Type of image that imgPath points to to load. If RTSTRUCT, must set originalImageDirPath
    originalImageDirPath
        Path to the directory containing the DICOMS for the original image the segmentation 
        was created from. 
    roiNames
        Identifier for which region(s) of interest to load from the total segmentation file
        
    Returns
    -------
    A dictionary of each of the ROIs and their name in the segmentation image as sitk.Image objects.
    '''
    #TODO: add in check for the file suffix
    #TODO: once med-imagetools SEG loading works, make this function return the CT as well as the segmentation
    #TODO: manage multiple segmentations in SEG, NIFTI, and NRRD

    if modality in ['SEG', 'seg']:
        # Loading SEG requires directory containing file, not the actual file path
        imgFolder, _ = os.path.split(imgPath)
        segHeader = pydicom.dcmread(imgPath, stop_before_pixels=True)
        roiName = segHeader.SegmentSequence[0].SegmentLabel
        return {roiName: loadDicomSITK(imgFolder)}
    
    elif modality in ['RTSTRUCT', 'rtstruct']:
        if originalImageDirPath == None:
            raise ValueError("Missing path to original image segmentation was taken from. RTSTRUCT loader requires original image.")
        else:
            return loadRTSTRUCTSITK(imgPath, originalImageDirPath, roiNames)
        
    elif modality in ['NIFTI', 'nifti']:
        return {"ROI": loadNiftiSITK(imgPath)}
    
    elif modality in ['NRRD', 'nrrd']:
        return {"ROI": loadNRRDSITK(imgPath)}
    
    else:
        raise ValueError('This segmentation modality is not supported. Must be one of RTSTRUCT, SEG, NIFTI, or NRRD')


def makeBinROIMask(image: sitk.Image,
                pixelType = sitk.sitkUInt32) -> sitk.Image:
    """Make a binarized mask from a segmentation image.

    Parameters
    ----------
    image
        Segmentation sitk.Image with background as 0s.
    
    pixelType
        Pixel type to cast binarized mask to.
    
    Returns
    -------
    A binarized map of the ROI as a sitk.Image.
    """

    imageArray = sitk.GetArrayFromImage(image)
    
    # Expecting a 3D image - segmentation can have a 4th dimension with just 1, so flatten to remove this
    if imageArray.ndim == 4:
        imageArray = np.squeeze(imageArray)
    
    binImageArray = (imageArray != 0).astype(float)

    maskedImage = sitk.GetImageFromArray(binImageArray)
    # maskedImage = sitk.Cast(maskedImage, sitk.sitkUInt32)
    
    return maskedImage


def flattenImage(image: sitk.Image) -> sitk.Image:
    """Remove axes of image with length one. (ex. shape is [1, 100, 256, 256])

    Parameters
    ----------
    image
        sitk.Image to flatten.
    
    Returns
    -------
    A sitk.Image with axes of length one removed.
    """
    imageArr = sitk.GetArrayFromImage(image)

    imageArr = np.squeeze(imageArr)

    return sitk.GetImageFromArray(imageArr)


def alignImages(originImage: sitk.Image, movingImage: sitk.Image):
    """Align movingImage to the originImage so origin and direction match

    Parameters
    ----------
    originImage
        sitk.Image to use to set direction and origin for the moving image

    movingImage
        sitk.Image to align to originImage
    
    Returns
    -------
    movingImage now aligned to originImage
    """
    movingImage.SetDirection(originImage.GetDirection())
    movingImage.SetOrigin(originImage.GetOrigin())
    movingImage.SetSpacing(originImage.GetSpacing())

    return movingImage


def displayImageSlice(imgArray, sliceIdx, cmap=plt.cm.Greys_r, dispMin = None, dispMax = None):
    # Function to display a 2D slice from a 3D image
    # By default, displays slice in greyscale with min and max range set to min and max value in the slice
    # imgArray - 3D ndarray object, the full array you'd like to display a slice of, must have slices as first dimension
    # sliceIdx - int, slice index from img_array you'd like to display
    # cmap - color map to use for plot, see https://matplotlib.org/stable/tutorials/colors/colormaps.html for options
    # dispMin - Value to use as min for cmap in display
    # dispMax - Value to use as max for cmap in display

    if dispMin == None:
        dispMin = imgArray.min()
    if dispMax == None:
        dispMax = imgArray.max()

    plt.imshow(imgArray[sliceIdx,:,:], cmap=cmap, vmin=dispMin, vmax=dispMax)
    plt.axis('off')


def padSEGtoMatchCT(ctFolderPath:str,
                    segFilePath:str,
                    ctImage:sitk.Image = None,
                    alignedSegImage:sitk.Image = None) -> sitk.Image:
    ''' Function to take a segmentation that doesn't have the same slice count as the base CT, maps it to the corresponding
        CT slices, and pads it with slices containing 0s so it maps properly onto the original image.

    Parameters
    ----------
    ctFolderPath
        Path to DICOM series folder containing all CT image files. Must be a directory.
    
    segFilePath
        Path to the DICOM SEG file that corresponds with CT in ctFolderPath that has incorrect slice count.

    ctImage
        Optional argument, CT image to align the padded segmentation image to. If None is passed, will be loaded in from ctFolderPath.
    
    alignedSegImage
        Optional argument, if image has already been loaded it can be passed in to be adjusted.
        Assumes that flattenImage and alignImages has already been run.
        If not passed, will use segFilePath to load the image.
    
    Returns
    -------
    Padded segmentation as a sitk.Image object. Will have the same dimensions as the CT.
    '''

    # Load the CT image to align the segmentation to if not passed as argument
    if ctImage == None:
        ctImage = loadDicomSITK(ctFolderPath)

    # Load in the segmentation image if not passed as argument
    if alignedSegImage == None:
        segImage = loadSegmentation(segFilePath, modality="SEG")
        # Segmentation contains extra axis, flatten to 3D by removing it
        segImage = flattenImage(segImage)
        # Segmentation has different origin, align it to the CT for proper feature extraction
        alignedSegImage = alignImages(ctImage, segImage)
    
    # Load in header information for the CT and SEG files
    ctSeries = Series(ctFolderPath)
    segWithHeader = pydicom.dcmread(segFilePath, stop_before_pixels=True)

    # Get the first and last reference ID for the slices of the CT that are in the SEG file
    lastSliceRef = segWithHeader.ReferencedSeriesSequence[0].ReferencedInstanceSequence[0].ReferencedSOPInstanceUID
    firstSliceRef = segWithHeader.ReferencedSeriesSequence[0].ReferencedInstanceSequence[-1].ReferencedSOPInstanceUID

    # Get the index of the reference IDs in the CT image
    firstSliceIdx = ctSeries['SOPInstanceUID'].index(firstSliceRef)
    lastSliceIdx = ctSeries['SOPInstanceUID'].index(lastSliceRef)

    # Convert the segmentation image to an array and pad with 0s so segmentation mask is in the correct indices
    arrSeg = sitk.GetArrayFromImage(alignedSegImage)
    padArrSeg = np.pad(arrSeg, (((firstSliceIdx, (ctSeries.data.shape[-1]-lastSliceIdx-1)), (0,0), (0,0))), 'constant', constant_values=(0))

    # Convert back to Image object
    paddedSegImage = sitk.GetImageFromArray(padArrSeg)
    paddedSegImage = alignImages(ctImage, paddedSegImage)

    return paddedSegImage


def getAcquisitionFileList(acquisitionInfo:pd.Series, 
                           summaryJSON:dict, 
                           imageDirPath:str):
    ''' Function to get all files for a single image acquisition (CT/MRI) to load the image.

    Parameters
    ----------
    acquisitionInfo
        pandas Series, row from the med-imagetools summary CSV containing info about the acquisition of interest

    summaryJSON
        dictionary, loaded JSON from med-imagetools output this acquisition is contained in

    imageDirPath
        Path to DICOM series folder containing all image files for the acquisition. Must be a directory. Will be appended
        to the acquisition file list entries
    
    Returns
    -------
    List of paths for files associated with a single image acquisition. Files have path to image directory added to front of entries.
    '''

    # Get info to access the acquisition file list from the summary JSON dictionary
    patID         = acquisitionInfo["patient_ID"]
    study_CT      = acquisitionInfo["study_CT"]
    series_CT     = acquisitionInfo["series_CT"]
    subseries_CT  = str(acquisitionInfo["subseries_CT"])
    instanceCount = acquisitionInfo["instances_CT"]

    # Get list of files for this acquisition
    acquisitionFileDict = summaryJSON[patID][study_CT][series_CT][subseries_CT]['instances']

    # Double check the right acquisition has been loaded
    if len(acquisitionFileDict) != instanceCount:
        raise ValueError("Acquisition matching failed.")

    # Combine the file list with the absolute path
    acquisitionFileList = [os.path.join(imageDirPath, imageFile) for imageFile in acquisitionFileDict.values()]
    # Have to reverse sort the list and convert tuple to match expected input of simpleITK loader function
    acquisitionFileList = sorted(tuple(acquisitionFileList), reverse=True)
    return acquisitionFileList


def getROIVoxelLabel(segImage:sitk.Image):
    ''' A function to find the non-zero value that identifies segmentation voxels in a loaded RTSTRUCT or SEG file.
    
    Parameters
    ----------
    segImage
        sitk.Image, a loaded segmentation image, should be binary with segmentation voxels as a non-zero value
    
    Returns
    -------
    labelValue
        int, the label value for the segmentation voxels

    '''
    # Convert segmentation image to a numpy array
    arrSeg = sitk.GetArrayFromImage(segImage)
    # Get all values that aren't 0 - these will identify the ROI
    roiVoxels = arrSeg[arrSeg != 0]
    # Confirm that all of these are the same value
    if np.all(roiVoxels == roiVoxels[0]):
        labelValue = roiVoxels[0]
        return labelValue
    else:
        raise ValueError("Multiple label values present in this segmentation. Must all be the same.")


def shuffleImage(imageToShuffle: sitk.Image):
    ''' Function to shuffle all pixel values in a sitk Image (developed for 3D, should work on 3D as well)

    Parameters
    ----------
    imageToShuffle
        sitk.Image, image to shuffle the pixels in
    
    Returns
    -------
        sitk.Image, image with all pixel values randomly shuffled with same dimensions as input image
    '''
    # Convert the image to an array
    arrImage = sitk.GetArrayFromImage(imageToShuffle)
    
    # Get array dimensions to reshape back to 
    imgDimensions = arrImage.shape

    # Flatten the 3D array to 1D so values can be shuffled
    flatArrImage = arrImage.flatten()

    # Shuffle the flat array
    np.random.shuffle(flatArrImage)

    # Reshape the array back into the original image dimensions
    shuffled3DArrImage = np.reshape(flatArrImage, imgDimensions)

    # Convert back to sitk Image
    shuffledImage = sitk.GetImageFromArray(shuffled3DArrImage)
    
    # Set the origin/direction/spacing from original image to shuffled image
    alignedShuffledImage = alignImages(imageToShuffle, shuffledImage)

    return alignedShuffledImage
