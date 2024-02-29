from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
import pandas as pd
# GS = GSRemoteProvider()

# GS_PREFIX = "orcestradata/"

# ids, = glob_wildcards("rawdata/radiomics/RADCURE/{patient_id}")

# # read in metadata/RADCURE_PatientIDs.csv 
# PATIENT_IDS = pd.read_csv("metadata/RADCURE_PatientIDs.csv")
# # get 'Patient ID` column and convert to list
# PATIENT_IDS = PATIENT_IDS["PatientID"].unique().tolist()
PATIENT_IDS = "RADCURE-0314", "RADCURE-0317"
NEG_CONTROLS= "randomized_full", "randomized_roi", "randomized_non_roi", \
              "shuffled_full", "shuffled_roi", "shuffled_non_roi", \
              "randomized_sampled_full", "randomized_sampled_roi", "randomized_sampled_non_roi"
PYRAD_SETTING = "scripts/radiomic_extraction/pyradiomics/pyrad_settings/uhn-radcure-challenge_params.yaml"

rule all:
    input: 
        radFeatures = expand("results/{patient_id}/readii_outputs/features/radiomicfeatures_{patient_id}.csv", patient_id=PATIENT_IDS),
        radFeatures_negcontrols = expand("results/{patient_id}/readii_outputs/features/radiomicfeatures_{negative_control}_{patient_id}.csv", patient_id=PATIENT_IDS, negative_control=NEG_CONTROLS),
        combined_radiomic_features = "results/snakemake_RADCURE/features/radiomicfeatures_RADCURE.csv",
        combined_negative_control_features = expand("results/snakemake_RADCURE/features/radiomicfeatures_{negative_control}_RADCURE.csv", negative_control=NEG_CONTROLS),
        maeObject = "results/RADCURE_readii_radiomic_MAE.rds"

rule runMedImageTools:
    input: 
        inputDir="rawdata/radiomics/RADCURE/{patient_id}"
    output: 
        csv_file="rawdata/radiomics/RADCURE/.imgtools/imgtools_{patient_id}.csv",
        json_file="rawdata/radiomics/RADCURE/.imgtools/imgtools_{patient_id}.json",
        outputDir=directory("data/med-imageout/{patient_id}")
    conda:
        "envs/medimage.yaml"
    retries:
        1
    threads:
        1
    shell:
        """
        autopipeline {input.inputDir} {output.outputDir} --update --dry_run
        """


rule runREADII:
    input:
        inputDir="rawdata/radiomics/RADCURE/{patient_id}",
        med_image_csv_file="rawdata/radiomics/RADCURE/.imgtools/imgtools_{patient_id}.csv"
    output:
        # negative_control_files,
        outputDir=directory("results/{patient_id}"),
        radFeatures="results/{patient_id}/readii_outputs/features/radiomicfeatures_{patient_id}.csv"

    params:
        roi_names="GTVp*",
        pyrad_setting=PYRAD_SETTING,
    conda:
        "envs/readii.yaml"
    retries:
        3
    threads:
        1
    log:
        "logs/readii/{patient_id}.log"
    shell:
        """
        readii {input.inputDir} {output.outputDir} \
        --roi_names {params.roi_names} \
        --pyradiomics_setting {params.pyrad_setting} \
        --update > {log} 2>&1
        """


rule runREADIINegativeControl:
    input:
        inputDir="rawdata/radiomics/RADCURE/{patient_id}",
        outputDir="results/{patient_id}",
    output:
        radFeatures_negcontrols = "results/{patient_id}/readii_outputs/features/radiomicfeatures_{negative_control}_{patient_id}.csv"
    params:
        roi_names="GTVp*",
        pyrad_setting=PYRAD_SETTING,
        negative_controls = "{negative_control}"
    conda:
        "envs/readii.yaml"
    retries:
        3
    threads:
        1
    log:
        "logs/readii/{patient_id}_{negative_control}.log"
    shell:
        """
        readii {input.inputDir} {input.outputDir} \
        --roi_names {params.roi_names} \
        --pyradiomics_setting {params.pyrad_setting} \
        --negative_controls {params.negative_controls} \
        --update > {log} 2>&1
        """

       
rule combineRadiomicFeatures:
    input:
        all_pat_radiomic_features = expand("results/{patient_id}/readii_outputs/features/radiomicfeatures_{patient_id}.csv", patient_id=PATIENT_IDS),
        # all_negative_control_csv = expand("results/{patient_id}/readii_outputs/features/radiomicfeatures_{negative_controls}_{patient_id}_.csv",  patient_id=PATIENT_IDS, negative_control=NEG_CONTROLS)
    output:
        combined_radiomic_features ="results/snakemake_RADCURE/features/radiomicfeatures_RADCURE.csv",
        # negative_control_radiomic_features="results/radiomic_output/snakemake_RADCURE/features/snakemake_RADCURE_negative_control_radiomic_features.csv"
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
        all_pat_negative_control_features = expand("results/{patient_id}/readii_outputs/features/radiomicfeatures_{negative_control}_{patient_id}.csv",  patient_id=PATIENT_IDS, negative_control="{negative_control}")
    output:
        combined_negative_control_features = "results/snakemake_RADCURE/features/radiomicfeatures_{negative_control}_RADCURE.csv"
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
        clinical="rawdata/clinical/clinical_RADCURE.xlsx",
        radiomic="results/snakemake_RADCURE/features",
    output:
        outputFileName="results/RADCURE_readii_radiomic_MAE.rds"
    params:
        pyrad=PYRAD_SETTING,
        findFeature="firstorder_10Percentile",
        clinicalPatIDCol="patient_id",
        radiomicPatIDCol="patient_ID",
    conda:
        "envs/makeMAE.yaml"
    script:
        "scripts/makeRadiogenomicMAE.R"


rule getClinicalData:
    output:
        clinical_file = "rawdata/clinical/clinical_RADCURE.xlsx"
    shell:
        """
        wget -O {output} "https://wiki.cancerimagingarchive.net/download/attachments/70226325/RADCURE_TCIA_Clinical%20June%2013%202023.xlsx?api=v2"
        """

