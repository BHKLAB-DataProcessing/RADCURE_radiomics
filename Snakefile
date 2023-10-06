rule all:
    input: 
        # "/Users/katyscott/Documents/RADCURE/.imgtools/imgtools_RADCURE.csv",
        # "/Users/katyscott/Documents/RADCURE/.imgtools/imgtools_RADCURE.json",
        "data/radiomic_output/snakemake_RADCURE/features/snakemake_RADCURE_radiomic_features.csv",
        "data/radiomic_output/snakemake_RADCURE/features/snakemake_RADCURE_negative_control_radiomic_features.csv",
        # "data/RADCURE_radiomic_MAE.rds"

rule run_medimagetools:
    input: 
        inputDir="/home/bioinf/bhklab/jermiah/RadioGenomics/pipelines/RADCURE_radiomics/orcestradata/radiomics/radcure_test_sample/images",
        outputDir="/home/bioinf/bhklab/jermiah/RadioGenomics/pipelines/RADCURE_radiomics/orcestradata/radiomics/radcure_test_sample/images",
    output: 
        "/home/bioinf/bhklab/jermiah/RadioGenomics/pipelines/RADCURE_radiomics/orcestradata/radiomics/radcure_test_sample/.imgtools/imgtools_images.csv",
        "/home/bioinf/bhklab/jermiah/RadioGenomics/pipelines/RADCURE_radiomics/orcestradata/radiomics/radcure_test_sample/.imgtools/imgtools_images.json"
    conda:
        "envs/medimage.yaml"
    shell:
        """
        autopipeline {input.inputDir} {input.outputDir} --dry_run
        """

rule extractRadiomicFeatures:
    input: 
        "/home/bioinf/bhklab/jermiah/RadioGenomics/pipelines/RADCURE_radiomics/orcestradata/radiomics/radcure_test_sample/.imgtools/imgtools_images.csv",
        "/home/bioinf/bhklab/jermiah/RadioGenomics/pipelines/RADCURE_radiomics/orcestradata/radiomics/radcure_test_sample/.imgtools/imgtools_images.json",
        config= "scripts/radiomic_extraction/RADCURE_config.yaml"
    output:
        "data/radiomic_output/snakemake_RADCURE/features/snakemake_RADCURE_radiomic_features.csv",
        "data/radiomic_output/snakemake_RADCURE/features/snakemake_RADCURE_negative_control_radiomic_features.csv"
    conda:
        "envs/radiomicExtraction.yaml"
    shell:
        "python3 scripts/radiomic_extraction/radiogenomic_pipeline.py {input.config} --update"
        

rule makeMAE:
    input:
        clinical="data/clinical_RADCURE.xlsx",
        radiomic="data/radiomic_output/snakemake_RADCURE/features/snakemake_RADCURE_radiomic_features.csv",
        pyrad="scripts/radiomic_extraction/pyradiomics/pyrad_settings/settings_original_allFeatures.yaml",
        negativecontrol="data/radiomic_output/snakemake_RADCURE/features/snakemake_RADCURE_negative_control_radiomic_features.csv"
    params:
        findFeature="firstorder_10Percentile",
        clinicalPatIDCol="patient_id",
        radiomicPatIDCol="patient_ID",
        outputFileName="data/RADCURE_radiomic_MAE.rds"
    container:
        "docker://jjjermiah/radcure_radiomics:0.1"
    conda:
        "envs/makeMAE.yaml"
    output:
        "data/RADCURE_radiomic_MAE.rds"
    script:
        "scripts/makeRadiogenomicMAE.R"

rule get_clinical:
    output:
        "data/clinical_RADCURE.xlsx"
    shell:
        """
        wget -O {output} "https://wiki.cancerimagingarchive.net/download/attachments/70226325/RADCURE_TCIA_Clinical%20June%2013%202023.xlsx?api=v2"
        """