#!/bin/bash
set -e

# run
# docker run -v $HOME/.vep:/home/crossbuild/.vep:ro crossbuild

VEP_CACHE_DIR="$HOME/.vep"

mkdir -p "$VEP_CACHE_DIR"

echo "Downloading VEP caches (GRCh37 and GRCh38)..."

# Example download for homo_sapiens_merged GRCh37
curl -L https://ftp.ensembl.org/pub/release-113/variation/indexed_vep_cache/homo_sapiens_merged/113_GRCh37.tar.gz | tar -xz -C "$VEP_CACHE_DIR"

# Example download for homo_sapiens_merged GRCh38
curl -L https://ftp.ensembl.org/pub/release-113/variation/indexed_vep_cache/homo_sapiens_merged/113_GRCh38.tar.gz | tar -xz -C "$VEP_CACHE_DIR"

# Add more cache downloads here as needed

echo "VEP caches downloaded to $VEP_CACHE_DIR"
