from argparse import ArgumentParser
import os

from readii.metadata import *
from readii.feature_extraction import *
from readii.utils import get_logger
import radiomics
from func_timeout import func_timeout, FunctionTimedOut
import yaml

logger = get_logger()


def parser():
    """Function to take command-line arguments and set them up for the pipeline run
    """
    parser = ArgumentParser("READII Feature Extraction Pipeline")

    # arguments
    parser.add_argument("data_directory", type=str,
                        help="Path to top-level directory of image dataset. Same as med-imagetools.")
    
    parser.add_argument("output_directory", type=str,
                       help="Path to output directory to save radiomic features and metadata.")
    
    parser.add_argument("--roi_names", type=str, default=None,
                        help="Name of region of interest in RTSTRUCT to perform extraction on.")
    
    parser.add_argument("--pyradiomics_setting", type=str, default=None,
                        help="Path to PyRadiomics configuration YAML file. If none provided, will use \
                              default in src/readii/data/.")
    
    parser.add_argument("--negative_control", type=str, default=None,
                        help="List of negative control types to run feature extraction on. Input as comma-separated list with no spaces.  \
                              Options: randomized_full,randomized_roi,randomized_non_roi,shuffled_full,shuffled_roi,shuffled_non_roi,randomized_sampled_full,randomized_sampled_roi,randomized_sampled_non_roi")

    parser.add_argument("--parallel", action="store_true",
                        help="Whether to run feature extraction in a parallel process. False by default.")

    parser.add_argument("--random_seed", type=int,
                        help="Value to set random seed to for reproducible negative controls")

    return parser.parse_known_args()[0]



def main():
    """Function to run READII radiomic feature extraction pipeline.
    """
    args = parser()
    pretty_args = '\n\t'.join([f"{k}: {v}" for k, v in vars(args).items()])
    logger.debug(
        f"Arguments:\n\t{pretty_args}"
    )
    
    # radiomics.setVerbosity(10)
    logger.info("Starting readii pipeline...")

     # Set up output directory
    outputDir = os.path.join(args.output_directory, "readii_outputs")
    if not os.path.exists(outputDir):
        logger.info(f"Directory {outputDir} does not exist. Creating...")
        os.makedirs(outputDir)
    else:
        logger.warning(f"Directory {outputDir} already exists. Will overwrite contents.")

    # Find med-imagetools output files
    logger.info("Finding med-imagetools outputs...")
    parentDirPath, datasetName = os.path.split(args.data_directory)
    imageFileListPath = os.path.join(parentDirPath + "/.imgtools/imgtools_" + datasetName + ".csv")
    if not os.path.exists(imageFileListPath):
        # Can we run med-imagetools in here?
        logger.error(
            f"Expected file {imageFileListPath} not found. Check the data_directory argument or run med-imagetools."
        )
        raise FileNotFoundError("Output for med-imagetools not found for this image set. Check the data_directory argument or run med-imagetools.")

    logger.info(f"Getting segmentation type...")
    try:
        # Get segType from imageFileList to generate the image metadata file and set up feature extraction
        segType = getSegmentationType(imageFileListPath)
    except RuntimeError as e:
        logger.error(str(e))
        logger.error("Feature extraction not complete.")
        exit()


    # Check if image metadata file has already been created
    imageMetadataPath = createImageMetadataFile(outputDir,
                                                parentDirPath,
                                                datasetName,
                                                segType,
                                                imageFileListPath)

    # Throwing TypeError: not all arguments converted during string formatting
    logger.info(f"Starting radiomic feature extraction for negative control: {args.negative_control}")
    try:
        # Timeout is 20 hours = 72,000 seconds
        ncRadiomicFeatures = func_timeout(timeout = 72000, func = radiomicFeatureExtraction, 
                                          args = (imageMetadataPath, parentDirPath, args.roi_names, args.pyradiomics_setting, 
                                                  outputDir, args.negative_control, args.random_seed, args.parallel))
    except FunctionTimedOut:
        logger.info(f"Passed radiomic features took too long to run. Removing Wavelet filter and running again.")
        
        # Log that this sample timed out for wavelet
        waveletLogFilePath = os.path.join(os.path.dirname(outputDir), "wavelet_logs/no_wavelet_features_list.log")
        os.makedirs(os.path.dirname(waveletLogFilePath), exist_ok=True)
        if not os.path.exists(waveletLogFilePath):
            waveletLogFile = open(waveletLogFilePath, "w")
        else:
            waveletLogFile = open(waveletLogFilePath, "a")
        
        waveletLogFile.write(datasetName)
        waveletLogFile.close()

        # Update the pyradiomics settings by removing the Wavelet filtering
        with open(args.pyradiomics_setting, 'r') as pyradFile:
            dictPyradiomicsSettings = yaml.safe_load(pyradFile)

        # Drop wavelet filter
        del dictPyradiomicsSettings['imageType']['Wavelet']


        ncRadiomicFeatures = radiomicFeatureExtraction(imageMetadataPath = imageMetadataPath,
                                                       imageDirPath = parentDirPath,
                                                       roiNames = args.roi_names,
                                                       pyradiomicsParamFilePath = dictPyradiomicsSettings,
                                                       outputDirPath = outputDir,
                                                       negativeControl = args.negative_control,
                                                       randomSeed=args.random_seed,
                                                       parallel = args.parallel)
    except Exception as e:
        logger.error(e)
        logger.error("Feature extraction not complete.")

    logger.info("Pipeline complete.")

if __name__ == "__main__":
    main()
