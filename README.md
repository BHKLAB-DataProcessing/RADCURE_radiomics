# title

# Running pipeline

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


