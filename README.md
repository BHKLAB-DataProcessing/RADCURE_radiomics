# title

# Running pipeline


### Dependencies

Creating conda environment with snakemake


```bash
conda create -n snakemake -c conda-forge -c bioconda snakemake # (optional: singularity)
```


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



