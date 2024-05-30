from argparse import ArgumentParser
import os

from readii.metadata import *
from readii.feature_extraction import *

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

     # Set up output directory
    outputDir = os.path.join(args.output_directory, "readii_outputs")
    if not os.path.exists(outputDir):
        print("Creating output directory:", outputDir)
        os.makedirs(outputDir)

    # Find med-imagetools output files
    print("Finding med-imagetools outputs...")
    parentDirPath, datasetName = os.path.split(args.data_directory)
    imageFileListPath = os.path.join(parentDirPath + "/.imgtools/imgtools_" + datasetName + ".csv")
    if not os.path.exists(imageFileListPath):
        # Can we run med-imagetools in here?
        raise FileNotFoundError("Output for med-imagetools not found for this image set. Check the data_directory argument or run med-imagetools.")

    print("Getting segmentation type...")
    try:
        # Get segType from imageFileList to generate the image metadata file and set up feature extraction
        segType = getSegmentationType(imageFileListPath)
    except RuntimeError as e:
        print(str(e))
        print("Feature extraction not complete.")
        exit()

    # Check if image metadata file has already been created
    imageMetadataPath = os.path.join(outputDir, "ct_to_seg_match_list_" + datasetName + ".csv")
    print("Matching CT to segmentations...")
    # Generate image metadata file by matching CT and segmentations in imageFileList from med-imagetools
    matchCTtoSegmentation(imgFileListPath = imageFileListPath,
                            segType = segType,
                            outputDirPath = outputDir)


    print("Starting radiomic feature extraction...")
    ncRadFeatOutPath = os.path.join(outputDir, "features/", "radiomicfeatures_" + args.negative_control + "_" + datasetName + ".csv")

    print("Starting radiomic feature extraction for negative control: ", args.negative_control)
    ncRadiomicFeatures = radiomicFeatureExtraction(imageMetadataPath = imageMetadataPath,
                                                   imageDirPath = parentDirPath,
                                                   roiNames = args.roi_names,
                                                   pyradiomicsParamFilePath = args.pyradiomics_setting,
                                                   outputDirPath = outputDir,
                                                   negativeControl = args.negative_control,
                                                   randomSeed=args.random_seed,
                                                   parallel = args.parallel)

if __name__ == "__main__":
    main()