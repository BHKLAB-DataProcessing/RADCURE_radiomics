rule all:
    input: 
        "/Users/katyscott/Documents/RADCURE/.imgtools/imgtools_RADCURE.csv",
        "/Users/katyscott/Documents/RADCURE/.imgtools/imgtools_RADCURE.json",
        "data/RADCURE_radiomic_MAE.rds"


rule run_medimagetools:
    input: 
        inputDir="/Users/katyscott/Documents/RADCURE/RADCURE",
        outputDir="data/med-imageout"
    output: 
        "/Users/katyscott/Documents/RADCURE/.imgtools/imgtools_RADCURE.csv",
        "/Users/katyscott/Documents/RADCURE/.imgtools/imgtools_RADCURE.json"
    conda:
        "envs/medimage.yaml"
    shell:
        """
        autopipeline {input.inputDir} {input.outputDir} --dry_run
        """

rule makeMAE:
    input:
        clinical="data/clinical_RADCURE.xlsx",
        radiomic="data/TCIA_RADCUREv1_radiomic_features.csv",
        pyrad="data/settings_allfilters_allfeatures.yaml"
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