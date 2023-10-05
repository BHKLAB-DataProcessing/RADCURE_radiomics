# Script to run the radiomic feature extraction portion of radiogenomic_pipeline
import argparse
import os
import sys
import pandas as pd
from shutil import copy
import yaml
from scalpl import Cut
import json

from pyradiomics.radiomic_feature_extraction_parallel import ctsegRadiomicFeatureExtractionParallel
from pyradiomics.negative_control_radiomic_feature_extraction_parallel import ctsegNegativeControlRadiomicFeatureExtractionParallel
from utils.find_segmentations import * 
from utils.quality_control import QualityControl

def main():
    config = Cut(parse_args())    
   
    # print config as a prettified dictionary
    # print(json.dumps(config.data, indent=4))
    UPDATE =  config['update']
    
    OUTPUTDIR = config['file_paths.output_dir']
    outputDir = OUTPUTDIR

    if not os.path.exists(outputDir):
        print("Setting up output directory: ", outputDir)
        os.makedirs(outputDir)

    # Copy the config file to the output dir so settings can be viewed
    # copy(CONFIGFILE, outputDir)
    # INSTEAD: write the config dictionary to a yaml file in the output directory
    with open(os.path.join(outputDir, "config.yaml"), 'w') as f:
        yaml.dump(config.data, f)
    
    
    # Check for .imgtools dir and dicoms file in the image_dir
    print("Loading med-imagetools outputs.")
    imageDirPath, ctDirName = os.path.split(config['file_paths.image_dir'])    

    dicomSummaryFile = os.path.join(imageDirPath + '/.imgtools/imgtools_' + ctDirName + '.csv')
    dicomSummaryJSON = os.path.join(imageDirPath + '/.imgtools/imgtools_' + ctDirName + '.json')
    if not os.path.exists(dicomSummaryFile):
        print(dicomSummaryFile)
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
        print("Generating image + segmentation summary file.")
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

    print("Radiomic feature extraction parameters:")
    print("summaryFilePath: ", summaryFilePath)
    print("pyradiomicsParamFilePath: ", pyradiomicsParamFilePath)
    print("idColumnName: ", idColumnName)
    print("segmentationDirPath: ", segmentationDirPath)
    print("segmentationLabel: ", segmentationLabel)
    print("radOutputFilePath: ", radOutputFilePath)
    print("roiNames: ", roiNames)
    print("parallel: ", parallel)
    
    
    # Check if radiomic features are already present and user has not requested them to be updated
    if os.path.exists(radOutputFilePath) and not UPDATE:
        print("Radiomic features have already been extracted. Loading from existing spreadsheet.")
        radiomicFeatures = pd.read_csv(radOutputFilePath)
    else:
        print("Extracting radiomic features")
        radiomicFeatures = ctsegRadiomicFeatureExtractionParallel(summaryFilePath, dicomSummaryJSON, pyradiomicsParamFilePath, idColumnName,
                                                                  imageDirPath, segmentationDirPath, segmentationLabel, roiNames, radOutputFilePath, parallel)
    
    #### NEGATIVE CONTROL RADIOMIC FEATURE SELECTION ####
    if config['radiomic_extraction.negative_control']:
        ncRadOutputFilePath = os.path.join(outputDir, "features/" + config['meta.experiment'] + "_negative_control_radiomic_features.csv")
        if os.path.exists(ncRadOutputFilePath) and not UPDATE:
            print("Negative control features have already been extracted. Loading from existing spreadsheet.")
            negControlRadFeatures = pd.read_csv(ncRadOutputFilePath)
        else:
            print("Extracting negative control radiomic features.")
            negControlRadFeatures = ctsegNegativeControlRadiomicFeatureExtractionParallel(
                                                    summaryFilePath, dicomSummaryJSON, pyradiomicsParamFilePath, idColumnName,
                                                    imageDirPath, segmentationDirPath, segmentationLabel, roiNames, ncRadOutputFilePath, parallel)

# an argument parsing function that will use argparse to get all the arguments
# needed to run the pipeline
# returns a dictionary of the arguments and their values

def parse_args():
    parser = argparse.ArgumentParser(description='Process some cool radiomic features.')
    parser.add_argument("--configfile", type=str, help="Path to file with settings for radiogenomic pipeline")
    
    parser.add_argument("--experiment", help="Name of the experiment")
    parser.add_argument("--name", help="Name of the pipeline")
    
    parser.add_argument("--top_dir", help="Path to the top directory")
    parser.add_argument("--image_dir", help="Path to the image directory")
    parser.add_argument("--segmentation_dir", help="Path to the segmentation directory")
    parser.add_argument("--output_dir", help="Path to the output directory")
    
    parser.add_argument("--clinical_id_column_label", help="Label of the ID column", default="Case_ID")
    parser.add_argument("--clinical_outcome_status", help="Label of the outcome status column", default=None)
    parser.add_argument("--clinical_na_value", help="Value to use for missing data", default="Unknown")
    
    parser.add_argument("--radiomic_id_column_label", help="Label of the ID column", default="patient_ID")
    parser.add_argument("--segmentation_label", help="Label of the segmentation", default=1)
    parser.add_argument("--negative_control", help="Whether to use negative control or not", default=True)
    parser.add_argument("--parallel", help="[True/False] Whether to use parallel processing or not", default=True)
    parser.add_argument("--segmentation_modality", help="[SEG, RTSTRUCT, NIFTI]  Modality of the segmentation")
    parser.add_argument("--roi_names", help="List of ROI names to extract features from")
    parser.add_argument("--pyrad_param_file", help="Path to pyradiomics parameter file")
    parser.add_argument("--quality_checks", help="Whether to perform quality checks or not")

    parser.add_argument("--update", help="Flag to force rerun all steps of pipeline.",
                        action="store_true")
    args = parser.parse_args()

    if args.configfile is not None:
        CONFIGFILE = args.configfile
        with open(CONFIGFILE, 'r') as f:
            configDict = yaml.load(f, Loader=yaml.FullLoader)
        print("Loading configuration file: ", CONFIGFILE, "\n") 
    elif all([
        args.experiment is not None,
        args.name is not None,
        args.image_dir is not None,
        args.segmentation_dir is not None,
        args.segmentation_modality is not None,
        args.roi_names is not None,
        args.pyrad_param_file is not None,
        args.quality_checks is not None]):
        
        ### OUTPUT SETUP ### removing this to redo 
        # # Make output directory for the experiment
        # if config['file_paths.output_dir'] == None:
        #     outputDir = os.path.join(config['file_paths.top_dir'], "radiogenomic_output/", config["meta.experiment"])
        # else:  
        #     outputDir = os.path.join(os.getcwd(), config['file_paths.output_dir'], config["meta.experiment"])
        if args.output_dir is not None:
            output_dir = os.path.join(
                args.output_dir,
                args.experiment)
        else:
            output_dir = os.path.join(
                os.dirname(args.image_dir),
                "radiogenomic_output",
                args.experiment)
        
        # make a dictionary of the arguments and their values
        configDict = {
            "meta": {
                "name": args.name,
                "experiment": args.experiment
            },
            "file_paths": {
                "image_dir": args.image_dir,
                "segmentation_dir": args.segmentation_dir,
                "clinical_feature_file": None,
                "output_dir": output_dir
            },
            "clinical_data": {
                "outcome_status": args.clinical_outcome_status,
                "id_column_label": args.clinical_id_column_label,
                "na_value": args.clinical_na_value
            },
            "radiomic_extraction": {
                "segmentation_modality": args.segmentation_modality,
                "id_column_label": args.radiomic_id_column_label,
                "segmentation_label": args.segmentation_label,
                "roi_names": args.roi_names,
                "pyrad_param_file": args.pyrad_param_file,
                "negative_control": args.negative_control,
                "parallel": args.parallel
            },
            "quality_checks": args.quality_checks,
            "update": args.update
        }
    else:
        # get list of missing mandatory arguments 
        missing_args = []
        if args.experiment is None:
            missing_args.append("experiment")
        if args.name is None:
            missing_args.append("name")
        if args.image_dir is None:
            missing_args.append("image_dir")
        if args.segmentation_dir is None:
            missing_args.append("segmentation_dir")
        if args.segmentation_modality is None:
            missing_args.append("segmentation_modality")
        if args.roi_names is None:
            missing_args.append("roi_names")
        if args.pyrad_param_file is None:
            missing_args.append("pyrad_param_file")
        if args.quality_checks is None:
            missing_args.append("quality_checks")
        
        # construct error message
        error_msg = "ERROR: Missing required arguments. Please provide either a config file or all required parameters in the command line. See --help for more information. \n"
        error_msg += "Missing mandatory arguments: " + ", ".join(missing_args) + "."
        # throw error
        sys.exit(error_msg)
    return configDict   


if __name__ == "__main__":
    main()
