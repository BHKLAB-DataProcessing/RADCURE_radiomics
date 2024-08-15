import pandas as pd
from pathlib import Path
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
HTTP = HTTPRemoteProvider()

#######################################################################
# CONFIGURATION
#######################################################################
# Need to configure default resources memory and disk space for each group to be used in the cluster

# read in metadata/RADCURE_PatientIDs.csv 
# PATIENT_IDS = pd.read_csv("metadata/RADCURE_RTSTRUCT_PatientIDs.csv")

# # get 'Patient ID` column and convert to list
# patients, desc, series, file = glob_wildcards("rawdata/radiomics/RADCURE/{patient_id}/{desc}/{series}/{file}.dcm")
# PATIENT_IDS = list(set(patients))

with open("metadata/patients_rtstruct.txt", 'r') as file:
    PATIENT_IDS = file.read().splitlines()

# PATIENT_IDS = "PATIENT_IDS[:100]"

configfile: "config/snakeconfig.yaml"

RANDOM_SEED = config['RANDOM_SEED']
READII_ROI_REGEX = config['READII_ROI_REGEX']
NEG_CONTROLS = config['NEG_CONTROLS']
PYRAD_SETTING = config['PYRAD_SETTING']


envs = Path("envs")
medimagetools_docker = "docker://bhklab/med-imagetools:1.2.0.2"
readii_docker = "docker://bhklab/readii:1.4.2"

rule all:
    input:
        radFeatures = expand("results/{patient_id}/readii_outputs/features/radiomicfeatures_{patient_id}.csv", patient_id=PATIENT_IDS),
        #ncRadFeatures =  expand("results/{patient_id}/readii_outputs/features/radiomicfeatures_{negative_control}_{patient_id}.csv", 
                                # patient_id=PATIENT_IDS, negative_control=NEG_CONTROLS)
        maeObject = "results/RADCURE_readii_radiomic_MAE.rds"
        

rule runMedImageTools:
    input: 
        inputDir="rawdata/radiomics/RADCURE/{patient_id}"
    output: 
        csv_file="rawdata/radiomics/RADCURE/.imgtools/imgtools_{patient_id}.csv",
        json_file="rawdata/radiomics/RADCURE/.imgtools/imgtools_{patient_id}.json",
        edge_file="rawdata/radiomics/RADCURE/.imgtools/imgtools_{patient_id}_edges.csv"
        # outputDir=temp(directory("data/med-imageout/{patient_id}"))
    group:
        "readii"
    conda:
        "envs/medimage.yaml"
    # container:
    #     readii_docker 
    threads: 
        1
    resources:
        mem_mb=500,
        disk_mb=500
    shell:
        """
        autopipeline {input.inputDir} /tmp --update --dry_run
        """


rule runREADII:
    input:
        inputDir="rawdata/radiomics/RADCURE/{patient_id}",
        med_image_csv_file="rawdata/radiomics/RADCURE/.imgtools/imgtools_{patient_id}.csv",
        med_image_csv_edge_file="rawdata/radiomics/RADCURE/.imgtools/imgtools_{patient_id}_edges.csv",
        PYRAD_SETTING = local(PYRAD_SETTING)
    output:
        radFeatures="results/{patient_id}/readii_outputs/features/radiomicfeatures_{patient_id}.csv"
    group:
        "readii"
    params:
        roi_names=READII_ROI_REGEX,
    conda:
        "envs/readii.yaml"
    # container:
    #     readii_docker
    resources:
        mem_mb=500,
        disk_mb=500
    threads: 
        1
    log:
        "logs/{patient_id}/readii/{patient_id}.log"
    shell:
        """
        OUTPUT_DIR=$(dirname $(dirname $(dirname {output.radFeatures})))
        python scripts/readii_pipeline.py \
            --data_directory {input.inputDir} \
            --output_directory $OUTPUT_DIR \
            --roi_names {params.roi_names} \
            --pyradiomics_setting {input.PYRAD_SETTING} | tee {log}
        """


rule runREADIINegativeControl:
    input:
        rules.runREADII.output.radFeatures, # force the negative control to wait for the radiomic features to be generated
        inputDir="rawdata/radiomics/RADCURE/{patient_id}",
        med_image_csv_file="rawdata/radiomics/RADCURE/.imgtools/imgtools_{patient_id}.csv",
        med_image_csv_edges_file="rawdata/radiomics/RADCURE/.imgtools/imgtools_{patient_id}_edges.csv",
        PYRAD_SETTING = local(PYRAD_SETTING),
    output:
        radFeatures_negcontrols = "results/{patient_id}/readii_outputs/features/radiomicfeatures_{negative_control}_{patient_id}.csv"
    group:
        "readii"
    params:
        RANDOM_SEED = RANDOM_SEED,
        roi_names=READII_ROI_REGEX,
        negative_controls = "{negative_control}"
    conda:
        "envs/readii.yaml"
    # container:
    #     readii_docker
    threads: 
        1
    resources:
        mem_mb=500,
        disk_mb=500
    log:
        "logs/{patient_id}/readii/{patient_id}_{negative_control}.log"
    shell:
        """
        OUTPUT_DIR=$(dirname $(dirname $(dirname {output.radFeatures_negcontrols})))
        python scripts/readii_negative_control_pipeline.py \
            --data_directory {input.inputDir} \
            --output_directory $OUTPUT_DIR \
            --roi_names {params.roi_names}  \
            --pyradiomics_setting {input.PYRAD_SETTING} \
            --negative_control {params.negative_controls} \
            --random_seed {params.RANDOM_SEED}
        """

       
rule combineRadiomicFeatures:
    input:
        all_pat_radiomic_features = expand("results/{patient_id}/readii_outputs/features/radiomicfeatures_{patient_id}.csv", patient_id=PATIENT_IDS)
    output:
        combined_radiomic_features ="results/snakemake_RADCURE/features/radiomicfeatures_RADCURE.csv",
        radiomicDir = directory("results/snakemake_RADCURE/features")
    run:
    # combine all csvs from each input into respective output files, only include the header once 
        with open(output.combined_radiomic_features, "w") as radiomic_features:
            for i, radiomic_csv in enumerate(input.all_pat_radiomic_features):
                with open(radiomic_csv) as radiomic_csv:
                    if i == 0:
                        radiomic_features.write(radiomic_csv.read())
                    else:
                        radiomic_csv.readline()
                        radiomic_features.write(radiomic_csv.read())



rule combineNegativeControlFeatures:
    input:
        all_pat_negative_control_features = expand(
            "results/{patient_id}/readii_outputs/features/radiomicfeatures_{negative_control}_{patient_id}.csv",  
            patient_id=PATIENT_IDS, negative_control="{negative_control}")
        # radFeatures_negcontrols = "results/{patient_id}/readii_outputs/features/radiomicfeatures_{negative_control}_{patient_id}.csv"
    output:
        combined_negative_control_features = "results/snakemake_RADCURE/features/radiomicfeatures_{negative_control}_RADCURE.csv",
    run:
    # combine all csvs from each input into respective output files, only include the header once 
        with open(output.combined_negative_control_features, "w") as nc_radiomic_features:
            for i, nc_radiomic_csv in enumerate(input.all_pat_negative_control_features):
                with open(nc_radiomic_csv) as nc_radiomic_csv:
                    if i == 0:
                        nc_radiomic_features.write(nc_radiomic_csv.read())
                    else:
                        nc_radiomic_csv.readline()
                        nc_radiomic_features.write(nc_radiomic_csv.read())


rule makeMAE:
    input:
        combined_negative_control_features = expand("results/snakemake_RADCURE/features/radiomicfeatures_{negative_control}_RADCURE.csv", negative_control=NEG_CONTROLS),
        combined_radiomic_features ="results/snakemake_RADCURE/features/radiomicfeatures_RADCURE.csv",
        clinical="rawdata/clinical/clinical_RADCURE.xlsx",
        PYRAD_SETTING = local(PYRAD_SETTING),
        # radiomicDir="results/snakemake_RADCURE/features",
    output:
        outputFileName="results/RADCURE_readii_radiomic_MAE.rds"
    params:
        findFeature="firstorder_10Percentile",
        clinicalPatIDCol="patient_id",
        radiomicPatIDCol="patient_ID",
    conda:
        "envs/makeMAE.yaml"
    shell:
        """
        ./scripts/makeRadiogenomicMAE.R {input.clinical} {input.combined_radiomic_features} {output.outputFileName} {input.PYRAD_SETTING} \
            --radiomic_find_feature {params.findFeature} \
            --clinical_patient_id {params.clinicalPatIDCol} \
            --radiomic_patient_id {params.radiomicPatIDCol}
        """


rule getClinicalData:
    input:
        file = HTTP.remote("https://wiki.cancerimagingarchive.net/download/attachments/70226325/RADCURE_TCIA_Clinical%20June%2013%202023.xlsx?api=v2")
    output:
        clinical_file = "rawdata/clinical/clinical_RADCURE.xlsx"
    shell:
        """
        mv {input.file} {output.clinical_file}
        """

