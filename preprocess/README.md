# Crossbuild analysis preprocessing Docker image

This Docker image bundles all tools required for crossbuild genomic analysis, including:

- bcftools with the liftover plugin to remap VCF coordinates
- Ensembl Variant Effect Predictor (VEP) pre-installed with required Perl modules
- CrossMap Python tool for coordinate liftover
- Supporting system libraries and utilities

# Quick start

Once you clone this repository, simply:

Build the Docker image:

```bash
docker build -t crossbuild .
```

Download VEP caches (required for offline annotation):

Inside the running container or on the host, run:

```bash
vep_install -a cf -s homo_sapiens -y GRCh38 -c /path/to/cache_dir
```

```bash
vep_install -a cf -s homo_sapiens -y GRCh37 -c /path/to/cache_dir
```
This downloads the cache files needed by VEP for offline use.

Use the tools inside the container:

- Run bcftools +liftover for VCF liftover operations.
- Run vep for variant effect prediction.
- Run CrossMap for alternative coordinate conversions.
