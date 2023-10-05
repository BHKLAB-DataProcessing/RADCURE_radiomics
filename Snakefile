from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
GS = GSRemoteProvider()

GS_PREFIX = "orcestradata/radiomics/radcure_test_sample/images/"

# create patient ids for set range and make sure ID # is 4 digits long
PATIENT_IDS = [f"RADCURE-{str(i).zfill(4)}" for i in [20,65,99,112]]

rule all:
    input: 
        GS.remote(GS_PREFIX + "object_output/RADCURE_radiomic_MAE.rds")
        # imagesjson = expand(".imgtools/imgtools_{patient_id}.json", patient_id=PATIENT_IDS),  
        # imagescsv = expand(".imgtools/imgtools_{patient_id}.csv", patient_id=PATIENT_IDS),

rule run_medimagetools:
    input: 
        inputDir="{patient_id}"
    output: 
        csv_file=".imgtools/imgtools_{patient_id}.csv",
        json_file=".imgtools/imgtools_{patient_id}.json",
        outputDir=directory("data/med-imageout/{patient_id}")
    conda:
        "envs/medimage.yaml"
    shell:
        """
        autopipeline {input.inputDir} {output.outputDir} --update --dry_run
        """

rule extractRadiomicFeatures:
    input: 
        imagescsv = ".imgtools/imgtools_{patient_id}.csv",
        imagesjson = ".imgtools/imgtools_{patient_id}.json"        
    output:
        features="radiomic_output/{patient_id}/features/snakemake_RADCURE_radiomic_features.csv",
        negative_control_features="radiomic_output/{patient_id}/features/snakemake_RADCURE_negative_control_radiomic_features.csv"
    conda:
        "envs/radiomicExtraction.yaml"
    params:
        config = "scripts/radiomic_extraction/RADCURE_config.yaml", 
    shell:
        "python3 scripts/radiomic_extraction/radiogenomic_pipeline.py {params.config}"
       
rule combineRadiomicFeatures:
    input:
        all_radiomic_csv = expand("radiomic_output/{patient_id}/features/snakemake_RADCURE_radiomic_features.csv",  patient_id=PATIENT_IDS),
        all_negative_control_csv = expand("radiomic_output/{patient_id}/features/snakemake_RADCURE_negative_control_radiomic_features.csv",  patient_id=PATIENT_IDS)
    output:
        radiomic_features="radiomic_output/snakemake_RADCURE/features/snakemake_RADCURE_radiomic_features.csv",
        negative_control_radiomic_features="radiomic_output/snakemake_RADCURE/features/snakemake_RADCURE_negative_control_radiomic_features.csv"
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
        clinical=GS.remote(GS_PREFIX + "clinical/clinical_RADCURE.xlsx"),
        radiomic=GS.remote(GS_PREFIX + "radiomic_output/snakemake_RADCURE/features/snakemake_RADCURE_radiomic_features.csv"),
        negativecontrol=GS.remote(GS_PREFIX + "radiomic_output/snakemake_RADCURE/features/snakemake_RADCURE_negative_control_radiomic_features.csv")
    output:
        GS.remote(GS_PREFIX + "object_output/RADCURE_radiomic_MAE.rds")
    params:
        pyrad="scripts/radiomic_extraction/pyradiomics/pyrad_settings/settings_original_allFeatures.yaml",
        findFeature="firstorder_10Percentile",
        clinicalPatIDCol="patient_id",
        radiomicPatIDCol="patient_ID",
        outputFileName="data/RADCURE_radiomic_MAE.rds"
    container:
        "docker://jjjermiah/radcure_radiomics:0.1"
    conda:
        "envs/makeMAE.yaml"
    script:
        "scripts/makeRadiogenomicMAE.R"

rule get_clinical:
    output:
        GS.remote(GS_PREFIX + "clinical/clinical_RADCURE.xlsx")
    shell:
        """
        wget -O {output} "https://wiki.cancerimagingarchive.net/download/attachments/70226325/RADCURE_TCIA_Clinical%20June%2013%202023.xlsx?api=v2"
        """