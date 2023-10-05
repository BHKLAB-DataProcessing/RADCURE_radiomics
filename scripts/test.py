import argparse
import sys

def main():
    parser = argparse.ArgumentParser(description='Process some cool radiomic features.')
        
    # Required arguments
    # either a --configfile [path to config file] or 
    # --config [PARAMETERS IN COMMAND LINE]
    # 
    # where [PARAMETERS IN COMMAND LINE] is:
    #       --experiment_name [experiment name] 
    #       --top_dir [path to top directory]
    #       --image_dir [path to image directory] 
    #       --segmentation_dir [path to segmentation directory] 
    #       --output_dir [path to output directory]
    #       --id_column_label [???]
    #       --segmentation_modality [SEG, RTSTRUCT, or NIFTI]
    #       --segmentation_label [???]
    #       --roi_names [list of ROI names to extract features from]
    #       --pyrad_param_file [path to pyradiomics parameter file]
    #       --negative_control [True or False]
    #       --parallel [True or False]
    #       --quality_checks [True or False]
    
    # Optional arguments
    #       --update [True or False]
   
    # if --configfile is not None:
    #     # load config file
    #     # check for all required parameters
    #     # if any are missing, throw an error
    #     # if all are present, run the pipeline
    #     pass
    # elif --config is not None:
    #     # check for all required parameters
    #     # if any are missing, throw an error
    #     # if all are present, run the pipeline
    #     pass
    # else:
    #     # throw an error
    #     pass

    # check if configfile is not None
    # if it is not None, load the config file
    # if it is None, check if the other required parameters are present
    # if they are, run the pipeline
    # if they are not, throw an error    
    parser.add_argument("--configfile", type=str, help="Path to file with settings for radiogenomic pipeline")
    parser.add_argument("--experiment_name", help="Name of the experiment")
    parser.add_argument("--top_dir", help="Path to the top directory")
    parser.add_argument("--image_dir", help="Path to the image directory")
    parser.add_argument("--segmentation_dir", help="Path to the segmentation directory")
    parser.add_argument("--output_dir", help="Path to the output directory")
    parser.add_argument("--id_column_label", help="Label of the ID column")
    parser.add_argument("--segmentation_modality", help="[SEG, RTSTRUCT, NIFTI]  Modality of the segmentation")
    parser.add_argument("--segmentation_label", help="Label of the segmentation")
    parser.add_argument("--roi_names", help="List of ROI names to extract features from")
    parser.add_argument("--pyrad_param_file", help="Path to pyradiomics parameter file")
    parser.add_argument("--negative_control", help="Whether to use negative control or not")
    parser.add_argument("--parallel", help="[True/False] Whether to use parallel processing or not")
    parser.add_argument("--quality_checks", help="Whether to perform quality checks or not")
    parser.add_argument("--update", help="Whether to update the pipeline or not")
    
    args = parser.parse_args()

    if args.configfile is not None:
        # load config file
        # check for all required parameters
        # if any are missing, throw an error
        # if all are present, run the pipeline
        pass
    elif all([args.experiment_name, args.top_dir, args.image_dir, args.segmentation_dir, args.output_dir, args.id_column_label, args.segmentation_modality, args.segmentation_label, args.roi_names, args.pyrad_param_file, args.negative_control, args.parallel, args.quality_checks]):
        # check for all required parameters
        # if all are present, run the pipeline
        pass
    else:
        # throw an error
        sys.exit("ERROR: Missing required arguments. Please provide either a config file or all required parameters in the command line. See --help for more information.")
        
if __name__ == '__main__':
    main()
