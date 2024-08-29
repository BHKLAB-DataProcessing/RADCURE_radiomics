from argparse import ArgumentParser
import os

from readii.metadata import *
from readii.feature_extraction import *
from readii.utils import get_logger

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
    
    parser.add_argument("--parallel", action="store_true",
                        help="Whether to run feature extraction in a parallel process. False by default.")

    return parser.parse_known_args()[0]



def main():
    """Function to run READII radiomic feature extraction pipeline.
    """
    args = parser()
    pretty_args = '\n\t'.join([f"{k}: {v}" for k, v in vars(args).items()])
    logger.debug(
        f"Arguments:\n\t{pretty_args}"
    )

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


    logger.info(f"Starting radiomic feature extraction...")
    radiomicFeatures = radiomicFeatureExtraction(imageMetadataPath = imageMetadataPath,
                                                 imageDirPath = parentDirPath,
                                                 roiNames = args.roi_names,
                                                 pyradiomicsParamFilePath = args.pyradiomics_setting,
                                                 outputDirPath = outputDir,
                                                 parallel = args.parallel)

    logger.info(f"Pipeline complete.")

if __name__ == "__main__":
    main()