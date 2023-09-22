import pandas as pd
import glob

from pyradiomics.helper_functions import * 

def findCTAndSegmentation(dicomFileListPath: str,
                          segType: str = "SEG",
                          outputFilePath: str = None):
    """From full list of image files, extract CT and corresponding SEG files and create new table.
    One row of the table contains both the CT and segmentation data for one patient.
    This function currently assumes there is one segmentation for each patient.

    Parameters
    ----------
    dicomFileListPath
        Path to csv containing list of DICOM directories in the dataset. 
        Expecting output from med-imagetools autopipeline .imgtools_[dataset]
    segType
        Type of file segmentation is in. Can be SEG or RTSTRUCT.
    outputFilePath
        Optional file path to save the dataframe to as a csv.
    
    Returns
    -------
    pandas dataframe containing the CT and corresponding segmentation data for each patient

    Note: All subseries of CT will be kept in the list in this function
    """
    # Check that segmentation file type is acceptable
    if segType != "RTSTRUCT" and segType != "SEG":
        raise ValueError("Incorrect segmentation file type. Must be RTSTRUCT or SEG.")

    # Load in complete list of patient dicom directories of all modalities (output from med-imagetools crawl)
    fullDicomList = pd.read_csv(dicomFileListPath, index_col=0)

    # Extract all CT rows
    allCTRows = fullDicomList.loc[fullDicomList['modality'] == "CT"]

    # Extract all SEG rows
    allSEGRows = fullDicomList.loc[fullDicomList['modality'] == segType]

    # Merge the CT and SEG dataframes based on the CT ID (referenced in the SEG rows)
    # Uses only SEG keys, so no extra CTs are kept
    # If multiple CTs have the same ID, they are both included in this table
    samplesWSeg = allCTRows.merge(allSEGRows, how='right', 
                                  left_on=['series', 'patient_ID'], 
                                  right_on=['reference_ct', 'patient_ID'], 
                                  suffixes=('_CT','_seg'))

    # Sort dataframe by ascending patient ID value
    samplesWSeg.sort_values(by='patient_ID', inplace=True)

    # Save out the combined list
    if outputFilePath != None:
        # Make directory if it doesn't exist
        if not os.path.exists(os.path.dirname(outputFilePath)):
            os.mkdir(os.path.dirname(outputFilePath))
        # Save out feature set
        samplesWSeg.to_csv(outputFilePath, index=False)

    return samplesWSeg


def matchCTandSegmentation(ctFileListPath: str,
                           segFileListPath: str,
                           segType: str = 'NIFTI',
                           outputFilePath: str = None):
    """Function takes two separate csv files containing paths to the CT and segmentation files and matches
    them and outputs a merged table. One row of the table contains the shared information of the CT and
    segmentation for one patient.

    Parameters
    ----------
    ctFileListPath
        Path to csv containing list of DICOM directories for the dataset.
    segFileListPath
        Path to csv containing list of file names of the segmentations. Must be NIFTI, RTSTRUCT, or SEG files.
    segType
        Type of file segmentation is in. Can be NIFTI, RTSTRUCT, or SEG.
    outputFilePath
        Optional file path to save the dataframe to as a csv.
    
    Returns
    -------
    pandas dataframe containing the CT and corresponding segmentation data for each patient
    
    """
    
    if segType != 'NIFTI' and segType != "RTSTRUCT" and segType != "SEG":
        raise ValueError("Incorrect segmentation file type. Must be NIFTI, RTSTRUCT, or SEG.")
    
    # Load in the CT list and remove any dicoms that aren't CT
    fullDicomList = pd.read_csv(ctFileListPath, index_col=0)
    ctDicomList = fullDicomList.loc[fullDicomList['modality'] == "CT"]
    # Load in the segmentation list
    segFileList = pd.read_csv(segFileListPath)

    samplesWSeg = ctDicomList.merge(segFileList, how='right', 
                                    left_on=['patient_ID', 'instances'], 
                                    right_on=['patient_ID', 'dimension'], 
                                    suffixes=('_CT', '_seg'))

    # Save out the combined list
    if outputFilePath != None:
        # Make directory if it doesn't exist
        if not os.path.exists(os.path.dirname(outputFilePath)):
            os.mkdir(os.path.dirname(outputFilePath))
        # Save out feature set
        samplesWSeg.to_csv(outputFilePath, index=False)

    return samplesWSeg


def makeNIFTIList(dirPath: str,
                  outputFilePath: str = None
                  ):
    ''' Function to create dataframe containing list of NIFTI files and some useful attributes.
    Can be used to associate segmentations with the original CT in matchCTandSegmentation.
    Function assumes the names of the files are the sample IDs.

    Parameters
    ----------
    dirPath
        Path to the directory containing the NIFTI files to create the dataframe for. 
    outputFilePath
        Optional file path to save the dataframe to as a csv.

    Returns
    -------
    Dataframe containing the list of NIFTI file sample IDs, the modality, the number of slices, and the full file path.
    '''

    # Get list of nifti files in the dirPath, recursively checking directories
    fileList = glob.glob(os.path.join(dirPath, '**', '*.nii.gz'), recursive=True)

    # Initialize dataframe to store data for each NIFTI file
    pdFileInfo = pd.DataFrame(columns = ['patient_ID', 'modality', 'dimension', 'folder'])

    for niftiFile in fileList:
        # Load the nifti file to get dimensions
        image = loadNiftiSITK(niftiFile)

        # Get just the sampleID for the output table by removing the nifti file suffix
        _, fileNameOnly = os.path.split(niftiFile)
        sampleID = fileNameOnly.removesuffix('.nii.gz')

        # Get relative file path from the image directory to mimic the med-imagetools output
        # Get full path to nifti file
        relImgPath = niftiFile.removeprefix(dirPath + '/')

        # Combine sample ID, the modality, the number of slices in the image, and the full path 
        # to add to output dataframe
        tempRow = pd.DataFrame([sampleID, "NIFTI", image.GetSize()[-1], relImgPath]).transpose()
        tempRow.columns = ['patient_ID', 'modality', 'dimension', 'folder']
        pdFileInfo = pd.concat([pdFileInfo, tempRow])
    
    # Reset the index of the info dataframe by dropping the existing one that just has 0s
    pdFileInfo.reset_index(drop=True, inplace=True)

    # Save out NIFTI file descriptions
    if outputFilePath != None:
        # Make directory if it doesn't exist
        if not os.path.exists(os.path.dirname(outputFilePath)):
            os.mkdir(os.path.dirname(outputFilePath))
        # Save out feature set
        pdFileInfo.to_csv(outputFilePath, index=False)

    return pdFileInfo