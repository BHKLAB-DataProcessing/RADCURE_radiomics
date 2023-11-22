from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
import pandas as pd
GS = GSRemoteProvider()

GS_PREFIX = "orcestradata/"

# ids, = glob_wildcards("rawdata/radiomics/RADCURE/{patient_id}")

# # read in metadata/RADCURE_PatientIDs.csv 
PATIENT_IDS = pd.read_csv("metadata/RADCURE_PatientIDs.csv")
# # get 'Patient ID` column and convert to list
PATIENT_IDS = PATIENT_IDS["PatientID"].unique().tolist()
# PATIENT_IDS = "RADCURE-0020"
rule all:
    input: 
        # "results/radiomic_output/snakemake_RADCURE/features/snakemake_RADCURE_radiomic_features.csv",
        "results/RADCURE_radiomic_MAE.rds"
        # imagesjson = expand("rawdata/radiomics/RADCURE/.imgtools/imgtools_{patient_id}.json", patient_id=PATIENT_IDS[1:10]),  
        # imagescsv = expand(".imgtools/imgtools_{patient_id}.csv", patient_id=PATIENT_IDS),


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

rule extractRadiomicFeatures:
    input: 
        imagescsv = "rawdata/radiomics/RADCURE/.imgtools/imgtools_{patient_id}.csv",
        imagesjson = "rawdata/radiomics/RADCURE/.imgtools/imgtools_{patient_id}.json",
        image_dir = "rawdata/radiomics/RADCURE/{patient_id}",
        segmentation_dir = "rawdata/radiomics/RADCURE/{patient_id}"      
    output:
        features="results/radiomic_output/{patient_id}/features/{patient_id}_radiomic_features.csv",
        negative_control_features="results/radiomic_output/{patient_id}/features/{patient_id}_negative_control_radiomic_features.csv",
        output_dir=directory("results/radiomic_output/{patient_id}")
    conda:
        "envs/radiomicExtraction.yaml"
    retries: 3
    params:
        # config = "scripts/radiomic_extraction/RADCURE_config.yaml",
        #### OR 
        name = "snakemake_RADCURE",
        segmentation_modality = "RTSTRUCT",
        roi_names = "GTVp*",
        pyrad_param_file = "scripts/radiomic_extraction/pyradiomics/pyrad_settings/settings_original_allFeatures.yaml",
        quality_checks = "False"
    log:
        "logs/radiomic_extraction/{patient_id}.log"
    shell:
        # "python3 scripts/radiomic_extraction/radiogenomic_pipeline.py {params.config}"
        """
        python scripts/radiomic_extraction/radiogenomic_pipeline.py \
            --image_dir {input.image_dir} \
            --segmentation_dir $(dirname {input.segmentation_dir}) \
            --output_dir $(dirname {output.output_dir}) \
            --name {params.name} \
            --experiment {wildcards.patient_id} \
            --segmentation_modality {params.segmentation_modality} \
            --roi_names {params.roi_names} \
            --pyrad_param_file {params.pyrad_param_file} \
            --quality_checks {params.quality_checks} \
            --update > {log} 2>&1
        """
       
rule combineRadiomicFeatures:
    input:
        all_radiomic_csv = expand("results/radiomic_output/{patient_id}/features/{patient_id}_radiomic_features.csv", patient_id=PATIENT_IDS),
        all_negative_control_csv = expand("results/radiomic_output/{patient_id}/features/{patient_id}_negative_control_radiomic_features.csv",  patient_id=PATIENT_IDS)
    output:
        radiomic_features="results/radiomic_output/snakemake_RADCURE/features/snakemake_RADCURE_radiomic_features.csv",
        negative_control_radiomic_features="results/radiomic_output/snakemake_RADCURE/features/snakemake_RADCURE_negative_control_radiomic_features.csv"
    run:
    # combine all csvs from each input into respective output files, only include the header once 
        with open(output.radiomic_features, "w") as radiomic_features, open(output.negative_control_radiomic_features, "w") as negative_control_radiomic_features:
            for i, (radiomic_csv, negative_control_csv) in enumerate(zip(input.all_radiomic_csv, input.all_negative_control_csv)):
                with open(radiomic_csv) as radiomic_csv, open(negative_control_csv) as negative_control_csv:
                    if i == 0:
                        radiomic_features.write(radiomic_csv.read())
                        negative_control_radiomic_features.write(negative_control_csv.read())
                    else:
                        radiomic_csv.readline()
                        negative_control_csv.readline()
                        radiomic_features.write(radiomic_csv.read())
                        negative_control_radiomic_features.write(negative_control_csv.read())

rule makeMAE:
    input:
        clinical="clinical/clinical_RADCURE.xlsx",
        radiomic="results/radiomic_output/snakemake_RADCURE/features/snakemake_RADCURE_radiomic_features.csv",
        negativecontrol="results/radiomic_output/snakemake_RADCURE/features/snakemake_RADCURE_negative_control_radiomic_features.csv"
    output:
        outputFileName="results/RADCURE_radiomic_MAE.rds"
    params:
        pyrad="scripts/radiomic_extraction/pyradiomics/pyrad_settings/settings_original_allFeatures.yaml",
        findFeature="firstorder_10Percentile",
        clinicalPatIDCol="patient_id",
        radiomicPatIDCol="patient_ID",
    conda:
        "envs/makeMAE.yaml"
    script:
        "scripts/makeRadiogenomicMAE.R"

rule getClinicalData:
    output:
        clinical_file = "clinical/clinical_RADCURE.xlsx"
    shell:
        """
        wget -O {output} "https://wiki.cancerimagingarchive.net/download/attachments/70226325/RADCURE_TCIA_Clinical%20June%2013%202023.xlsx?api=v2"
        """

