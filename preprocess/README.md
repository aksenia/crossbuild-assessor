# Crossbuild analysis preprocessing Docker image

This Docker image bundles all tools required for crossbuild genomic analysis, including:

- bcftools with the liftover plugin to remap VCF coordinates
- Ensembl Variant Effect Predictor (VEP) pre-installed with required Perl modules
- CrossMap Python tool for coordinate liftover
- Supporting system libraries and utilities

# Quick start

Once you clone this repository, simply:

Build the image:

```bash
docker build -t crossbuild .
singularity build crossbuild.sif docker-daemon://crossbuild:latest
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

## Use the tools inside the container:

- Run bcftools +liftover for VCF liftover operations.
- Run vep for variant effect prediction.
- Run CrossMap for alternative coordinate conversions.

## Run the full data preprocessing pipeline: 

Set the `--cores` according to resources on your system, mind `-n` flag for dry-run.

in Singularity 

```bash
# to run safely without writing to home directory
mkdir tmp_home
# set the tmp_home as HOME var before running your jobs
singularity exec --no-home   -B $(pwd)/tmp_home:/tmp_home -B $(pwd)/data:/data -B $(pwd)/snake/config.yaml:/app/snake/config.yaml crossbuild.sif bash -c 'HOME=/tmp_home snakemake --snakefile /app/snake/Snakefile --configfile /app/snake/config.yaml -np --cores 1'
```
in Docker

```bash
docker run --rm -it  -v $(pwd)/data:/data -v $(pwd)/snake/config.yaml:/app/snake/config.yaml crossbuild snakemake --snakefile /app/snake/Snakefile --configfile /app/snake/config.yaml -np --cores 1
```

## Large datasets

For large datasets you might want to submit a job to the cluster, or if you want to only submit computationally intensive jobs such as vep rules, you can do targeted rule execution in snakemake: 

```bash
# your Docker/singularity command with correct mounts (see below) \
	'snakemake --snakefile /app/snake/Snakefile --configfile /app/snake/config.yaml --target-jobs vep_hg19:sample=sampleid  -p --cores 1'
```


## Directory mounting 

To keep the container small, we need to mount all the reference and data directories needed. For convenience one can use an utility script to produce the final command with correct mounts. All the paths are on the *host* system:  

in Singularity

```bash
singularity exec crossbuild.sif python build_command.py --tmp_home $(pwd)/tmp_home --data $(pwd)/data/sampleid --results $(pwd)/output/sampleid  --vep_cache_hg19 /path/to/VEP/cache_hg19   --vep_cache_hg38 /path/to/VEP/cache_hg38   --hg19_fa /path/to/hg19.fa   --hg38_fa /path/to/hg38.fa   --chain_file /path/to/hg19ToHg38.over.chain  --config_yaml $(pwd)/data/sampleid/config.yaml  --cores 4
```


in Docker

```bash
docker run --rm -it crossbuild python build_command.py --data $(pwd)/data/sampleid --results $(pwd)/output/sampleid  --vep_cache_hg19 /path/to/VEP/cache_hg19 --vep_cache_hg38 /path/to/VEP/cache_hg38 --hg19_fa /path/to/hg19.fa --hg38_fa /path/to/hg38.fa --chain_file /path/to/hg19ToHg38.over.chain --config_yaml $(pwd)/data/sampleid/config.yaml --cores 4   --engine docker
``` 
