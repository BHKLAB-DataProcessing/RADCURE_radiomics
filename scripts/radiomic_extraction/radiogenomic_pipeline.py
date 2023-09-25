# Script to run the radiomic feature extraction portion of radiogenomic_pipeline
import argparse
import os
import pandas as pd
from shutil import copy
import yaml
from scalpl import Cut

from pyradiomics.radiomic_feature_extraction_parallel import ctsegRadiomicFeatureExtractionParallel
from pyradiomics.negative_control_radiomic_feature_extraction_parallel import ctsegNegativeControlRadiomicFeatureExtractionParallel
from utils.find_segmentations import * 
from utils.quality_control import QualityControl

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("config", help="Path to file with settings for radiogenomic pipeline")
    parser.add_argument("--update", help="Flag to force rerun all steps of pipeline.",
                        action="store_true")
    args = parser.parse_args()

    CONFIGFILE = args.config
    # Flag to rerun everything even if it's already been run
    UPDATE = args.update

    # HERE FOR DEBUGGING, REMOVE WHEN DONE
    # CONFIGFILE = os.path.join(os.getcwd() + "/scripts/radiomic_extraction/RADCURE_config.yaml")
    # UPDATE = True

    # Load in run settings from config 
    with open(CONFIGFILE, 'r') as f:
        configDict = yaml.load(f, Loader=yaml.FullLoader)

    # Making config dictionary accessible with . operators
    config = Cut(configDict)

    ### OUTPUT SETUP ###
    # Make output directory for the experiment
    if config['file_paths.output_dir'] == None:
        outputDir = os.path.join(config['file_paths.top_dir'], "radiogenomic_output/", config["meta.experiment"])
    else:  
        outputDir = os.path.join(os.getcwd(), config['file_paths.output_dir'], config["meta.experiment"])

    if not os.path.exists(outputDir):
        os.makedirs(outputDir)

    # Copy the config file to the output dir so settings can be viewed
    copy(CONFIGFILE, outputDir)

    # Check for .imgtools dir and dicoms file in the image_dir
    imageDirPath, ctDirName = os.path.split(config['file_paths.image_dir'])    
    dicomSummaryFile = os.path.join(imageDirPath + '/.imgtools/imgtools_' + ctDirName + '.csv')
    dicomSummaryJSON = os.path.join(imageDirPath + '/.imgtools/imgtools_' + ctDirName + '.json')
    if not os.path.exists(dicomSummaryFile):
        raise FileNotFoundError("Output for med-imagetools not found for this image set. Check the image_dir argument or run med-imagetools.")
    
    ## QC Checks
    if config['quality_checks'] == 'True':
        # placeholder parameters for quality control checks
        qc_params = {
            'modality': '',
            'patient_ID': '',
            'kernel_list': list(),
            'slice_thickness_range': list(),
            'scan_length_range': list(),
            'spacing_range': list()
        }
        # Quality control takes two arguments: qc params and path to dicom images
        qcChecksDataFrame = QualityControl(qc_params=qc_params, imgPath=imageDirPath).run_checks()

    # Check if summary file listing the corresponding images and segmentations has been made yet
    # If not, make it
    summaryFilePath = os.path.join(outputDir, config['meta.experiment'] + "_seg_and_ct_dicom_list.csv")
    if not os.path.exists(summaryFilePath) or UPDATE:
        if config['radiomic_extraction.segmentation_modality'] == 'NIFTI':
            # NIFTI files won't be caught by med-imagetools, need to make a separate list of them
            niftiOut = os.path.join(outputDir, "segmentation_file_summary.csv")
            # Output file contains sample IDs, number of slices in image, and full path to NIFTI file
            makeNIFTIList(dirPath = os.path.join(config['file_paths.segmentation_dir']),
                                        outputFilePath = niftiOut)
            # Find the intersection of the segmentation and CT lists so correct files get loaded later
            matchCTandSegmentation(ctFileListPath=dicomSummaryFile,
                                   segFileListPath=niftiOut,
                                   segType=config['radiomic_extraction.segmentation_modality'],
                                   outputFilePath=summaryFilePath)

        else: 
            # Segmentation is SEG or RTSTRUCT
            # Find the intersection of the segmentation and CT lists so correct files get loaded later
            findCTAndSegmentation(dicomSummaryFile,
                                  segType = config['radiomic_extraction.segmentation_modality'],
                                  outputFilePath = summaryFilePath)


    ### RADIOMIC FEATURE EXTRACTION ###
    # Set up input variables for radiomic feature extraction
    # summaryFilePath = os.path.join(config['file_paths.top_dir'], config['file_paths.img_summary_file'])
    pyradiomicsParamFilePath = os.path.join(os.getcwd(), config['radiomic_extraction.pyrad_param_file'])
    idColumnName = config['radiomic_extraction.id_column_label']
    segmentationDirPath = config['file_paths.segmentation_dir']
    segmentationLabel = config['radiomic_extraction.segmentation_label']
    radOutputFilePath = os.path.join(outputDir, "features/" + config['meta.experiment'] + "_radiomic_features.csv")
    roiNames = config['radiomic_extraction.roi_names']
    parallel = config['radiomic_extraction.parallel']

    # Check if radiomic features are already present and user has not requested them to be updated
    if os.path.exists(radOutputFilePath) and not UPDATE:
        print("Radiomic features have already been extracted. Loading from existing spreadsheet.")
        radiomicFeatures = pd.read_csv(radOutputFilePath)

    else:
        radiomicFeatures = ctsegRadiomicFeatureExtractionParallel(summaryFilePath, dicomSummaryJSON, pyradiomicsParamFilePath, idColumnName,
                                                                  imageDirPath, segmentationDirPath, segmentationLabel, roiNames, radOutputFilePath, parallel)
    
    #### NEGATIVE CONTROL RADIOMIC FEATURE SELECTION ####
    if config['radiomic_extraction.negative_control']:
        ncRadOutputFilePath = os.path.join(outputDir, "features/" + config['meta.experiment'] + "_negative_control_radiomic_features.csv")
        if os.path.exists(ncRadOutputFilePath) and not UPDATE:
            print("Negative control features have already been extracted. Loading from existing spreadsheet.")
            negControlRadFeatures = pd.read_csv(ncRadOutputFilePath)
        else:
            negControlRadFeatures = ctsegNegativeControlRadiomicFeatureExtractionParallel(
                                                    summaryFilePath, dicomSummaryJSON, pyradiomicsParamFilePath, idColumnName,
                                                    imageDirPath, segmentationDirPath, segmentationLabel, roiNames, ncRadOutputFilePath, parallel)


if __name__ == "__main__":
    main()