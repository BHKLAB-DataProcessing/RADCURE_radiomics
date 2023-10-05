# title

# Running pipeline


### Dependencies

Creating conda environment with snakemake


```bash
conda create -n snakemake -c conda-forge -c bioconda snakemake # (optional: singularity)
```


## Running using buckets


The bucket used is (GCP Project ID = `orcestra-388613`):

```bash
orcestradata/radiomics/radcure_test_sample/
```


```bash
snakemake -c8 \
--use-conda \
--default-remote-provider=GS \
--default-remote-prefix=orcestradata/radiomics/radcure_test_sample/images \
--keep-remote 
```


where the `RADCURE-####` folders are all in `orcestradata/radiomics/radcure_test_sample/images`


### Using Singularity

* containerization technology that allows users to run lightweight, isolated environments for specific use-cases.
* initially developed for Linux
* using on MacOS requires a linux VM (https://docs.sylabs.io/guides/3.0/user-guide/installation.html#mac)

``` bash
snakemake -c1 --use-singularity
```

### Using conda environments

```bash
snakemake -c1 --use-conda
```


### Creating directed acyclic graph (DAG)

```bash
snakemake -c1 --dag | dot -Tsvg > dag.svg
```



# Testing individual scripts

```bash
image_dir="orcestradata/radiomics/radcure_test_sample/images/RADCURE-0020"
segmentation_dir="orcestradata/radiomics/radcure_test_sample/images"
output_dir="results/radiomic_output"
python scripts/radiomic_extraction/radiogenomic_pipeline.py \
    --name "snakemake_RADCURE_pipeline" \
    --experiment "snakemake_RADCURE" \
    --image_dir ${image_dir} \
    --segmentation_dir ${segmentation_dir} \
    --output_dir ${output_dir} \
    --segmentation_modality RTSTRUCT \
    --roi_names GTVp.* \
    --pyrad_param_file scripts/radiomic_extraction/pyradiomics/pyrad_settings/settings_original_allFeatures.yaml \
    --quality_checks True 
```


