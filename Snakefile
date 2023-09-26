from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
GS = GSRemoteProvider()

GS_PREFIX = "orcestradata/radiomics/radcure_test_sample/images/"

# create patient ids for set range and make sure ID # is 4 digits long
PATIENT_IDS = [f"RADCURE-{str(i).zfill(4)}" for i in [20,65,99,112]]

rule all:
    input: 
        GS.remote(GS_PREFIX + "/object_output/RADCURE_radiomic_MAE.rds")
        # "/Users/katyscott/Documents/RADCURE/.imgtools/imgtools_images.csv",
        # "/Users/katyscott/Documents/RADCURE/.imgtools/imgtools_images.json",
        # "data/radiomic_output/snakemake_RADCURE/features/snakemake_RADCURE_radiomic_features.csv",
        # "data/radiomic_output/snakemake_RADCURE/features/snakemake_RADCURE_negative_control_radiomic_features.csv",
        # "data/RADCURE_radiomic_MAE.rds"
        
rule run_medimagetools:
    input: 
        # inputDir="/Users/katyscott/Documents/RADCURE/RADCURE",
        inputDir=GS.remote(GS_PREFIX + "{patient_id}")
    output: 
        # "/Users/katyscott/Documents/RADCURE/.imgtools/imgtools_RADCURE.csv",
        # "/Users/katyscott/Documents/RADCURE/.imgtools/imgtools_RADCURE.json",
        GS.remote(GS_PREFIX + "/.imgtools/imgtools_{patient_id}.csv"),
        GS.remote(GS_PREFIX + "/.imgtools/imgtools_{patient_id}.json"),
        outputDir=directory("data/med-imageout/{patient_id}")
    conda:
        "envs/medimage.yaml"
    shell:
        """
        autopipeline {input} {output} --update --dry_run
        """

rule extractRadiomicFeatures:
    input: 
        config = "scripts/radiomic_extraction/RADCURE_config.yaml", 
        imagescsv = GS.remote(GS_PREFIX + "/.imgtools/imgtools_{patient_id}.csv"),
        imagesjson = GS.remote(GS_PREFIX + "/.imgtools/imgtools_{patient_id}.json")        
    output:
        GS.remote(GS_PREFIX + "/radiomic_output/{patient_id}/features/snakemake_RADCURE_radiomic_features.csv"),
        GS.remote(GS_PREFIX + "/radiomic_output/{patient_id}/features/snakemake_RADCURE_negative_control_radiomic_features.csv")
    conda:
        "envs/radiomicExtraction.yaml"
    shell:
        "python3 scripts/radiomic_extraction/radiogenomic_pipeline.py {input.config}"
       
rule combineRFE:
    input:
        all_radiomic_csv = GS.remote(expand(GS_PREFIX + "/radiomic_output/{patient_id}/features/snakemake_RADCURE_radiomic_features.csv",  patient_id=PATIENT_IDS)),
        all_negative_control_csv = GS.remote(expand(GS_PREFIX + "/radiomic_output/{patient_id}/features/snakemake_RADCURE_negative_control_radiomic_features.csv",  patient_id=PATIENT_IDS))
    output:
        GS.remote(GS_PREFIX + "/radiomic_output/snakemake_RADCURE/features/snakemake_RADCURE_radiomic_features.csv"),
        GS.remote(GS_PREFIX + "/radiomic_output/snakemake_RADCURE/features/snakemake_RADCURE_negative_control_radiomic_features.csv")
# rule all:
#     input: 
#         GS.remote(GS_PREFIX + "/object_output/RADCURE_radiomic_MAE.rds")
#         # "/Users/katyscott/Documents/RADCURE/.imgtools/imgtools_images.csv",
#         # "/Users/katyscott/Documents/RADCURE/.imgtools/imgtools_images.json",
#         # "data/radiomic_output/snakemake_RADCURE/features/snakemake_RADCURE_radiomic_features.csv",
#         # "data/radiomic_output/snakemake_RADCURE/features/snakemake_RADCURE_negative_control_radiomic_features.csv",
#         # "data/RADCURE_radiomic_MAE.rds"
        

# rule run_medimagetools:
#     input: 
#         # inputDir="/Users/katyscott/Documents/RADCURE/RADCURE",
#         inputDir=GS.remote(GS_PREFIX + "/images")
#     output: 
#         # "/Users/katyscott/Documents/RADCURE/.imgtools/imgtools_RADCURE.csv",
#         # "/Users/katyscott/Documents/RADCURE/.imgtools/imgtools_RADCURE.json",
#         GS.remote(GS_PREFIX + "/.imgtools/imgtools_images.csv"),
#         GS.remote(GS_PREFIX + "/.imgtools/imgtools_images.json"),
#         outputDir=directory("data/med-imageout")
#     conda:
#         "envs/medimage.yaml"
#     shell:
#         """
#         autopipeline {input.inputDir} {output.outputDir} --dry_run
#         """

# rule extractRadiomicFeatures:
#     input: 
#         config = "scripts/radiomic_extraction/RADCURE_config.yaml", 
#         imagescsv = GS.remote(GS_PREFIX + "/.imgtools/imgtools_images.csv"),
#         imagesjson = GS.remote(GS_PREFIX + "/.imgtools/imgtools_images.json"),        
#     output:
#         GS.remote(GS_PREFIX + "/radiomic_output/snakemake_RADCURE/features/snakemake_RADCURE_radiomic_features.csv"),
#         GS.remote(GS_PREFIX + "/radiomic_output/snakemake_RADCURE/features/snakemake_RADCURE_negative_control_radiomic_features.csv")
#     conda:
#         "envs/radiomicExtraction.yaml"
#     shell:
#         "python3 scripts/radiomic_extraction/radiogenomic_pipeline.py {input.config}"
       

rule makeMAE:
    input:
        pyrad="scripts/radiomic_extraction/pyradiomics/pyrad_settings/settings_original_allFeatures.yaml",
        clinical=GS.remote(GS_PREFIX + "/clinical/clinical_RADCURE.xlsx"),
        radiomic=GS.remote(GS_PREFIX + "/radiomic_output/snakemake_RADCURE/features/snakemake_RADCURE_radiomic_features.csv"),
        negativecontrol=GS.remote(GS_PREFIX + "/radiomic_output/snakemake_RADCURE/features/snakemake_RADCURE_negative_control_radiomic_features.csv")
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
        GS.remote(GS_PREFIX + "/object_output/RADCURE_radiomic_MAE.rds")
    script:
        "scripts/makeRadiogenomicMAE.R"

rule get_clinical:
    output:
        GS.remote(GS_PREFIX + "/clinical/clinical_RADCURE.xlsx")
    shell:
        """
        wget -O {output} "https://wiki.cancerimagingarchive.net/download/attachments/70226325/RADCURE_TCIA_Clinical%20June%2013%202023.xlsx?api=v2"
        """