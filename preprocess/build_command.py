import os
import argparse

def build_command(user_paths, cores=1, engine="singularity"):
    # Required mounts for both engines, absolute paths on host
    required_mounts = {
        "tmp_home": os.path.abspath("./tmp_home"),
        "data": os.path.abspath("./data/test-small"),
        "results": os.path.abspath("./data/test-small/preprocess"),
        "vep_cache_hg19": user_paths["vep_cache_hg19"],
        "vep_cache_hg38": user_paths["vep_cache_hg38"],
        "hg19_fa": user_paths["hg19_fa"],
        "hg19_fa_fai": user_paths["hg19_fa"] + ".fai",
        "hg38_fa": user_paths["hg38_fa"],
        "hg38_fa_fai": user_paths["hg38_fa"] + ".fai",
        "chain_file": user_paths["chain_file"],
        "config_yaml": os.path.abspath("./data/test-small/config.yaml")
    }

    if engine == "singularity":
        binds = [
            f"{required_mounts['tmp_home']}:/tmp_home",
            f"{required_mounts['data']}:/data",
            f"{required_mounts['results']}:/results",
            f"{required_mounts['vep_cache_hg19']}:/ref/VEP/cache_hg19",
            f"{required_mounts['vep_cache_hg38']}:/ref/VEP/cache_hg38",
            f"{required_mounts['hg19_fa']}:/ref/hg19.fa",
            f"{required_mounts['hg19_fa_fai']}:/ref/hg19.fa.fai",
            f"{required_mounts['hg38_fa']}:/ref/hg38.fa",
            f"{required_mounts['hg38_fa_fai']}:/ref/hg38.fa.fai",
            f"{required_mounts['chain_file']}:/ref/hg19ToHg38.over.chain",
            f"{required_mounts['config_yaml']}:/app/snake/config.yaml"
        ]

        binds_str = " -B ".join(binds)
        cmd = (
            f"singularity exec --no-home -B {binds_str} "
            f"crossbuild.sif bash -c 'HOME=/tmp_home snakemake "
            f"--snakefile /app/snake/Snakefile --configfile /app/snake/config.yaml "
            f"-p --cores {cores}'"
        )
    elif engine == "docker":
        # Docker uses -v for mounts
        binds = [
            f"{required_mounts['tmp_home']}:/tmp_home:rw",
            f"{required_mounts['data']}:/data:rw",
            f"{required_mounts['results']}:/results:rw",
            f"{required_mounts['vep_cache_hg19']}:/ref/VEP/cache_hg19:ro",
            f"{required_mounts['vep_cache_hg38']}:/ref/VEP/cache_hg38:ro",
            f"{required_mounts['hg19_fa']}:/ref/hg19.fa:ro",
            f"{required_mounts['hg19_fa_fai']}:/ref/hg19.fa.fai:ro",
            f"{required_mounts['hg38_fa']}:/ref/hg38.fa:ro",
            f"{required_mounts['hg38_fa_fai']}:/ref/hg38.fa.fai:ro",
            f"{required_mounts['chain_file']}:/ref/hg19ToHg38.over.chain:ro",
            f"{required_mounts['config_yaml']}:/app/snake/config.yaml:ro"
        ]

        binds_str = " ".join([f"-v {b}" for b in binds])
        container_image = "crossbuild:latest"  # Adjust docker image name/tag accordingly
        cmd = (
            f"docker run --rm -it {binds_str} "
            f"{container_image} bash -c 'HOME=/tmp_home snakemake "
            f"--snakefile /app/snake/Snakefile --configfile /app/snake/config.yaml "
            f"-p --cores {cores}'"
        )
    else:
        raise ValueError("Engine must be 'singularity' or 'docker'")

    return cmd


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Build container run command for Snakemake")
    parser.add_argument("--vep_cache_hg19", required=True, help="Absolute host path for VEP hg19 cache (mounted as /ref/VEP/cache_hg19)")
    parser.add_argument("--vep_cache_hg38", required=True, help="Absolute host path for VEP hg38 cache (mounted as /ref/VEP/cache_hg38)")
    parser.add_argument("--hg19_fa", required=True, help="Absolute host path for hg19 fasta (mounted as /ref/hg19.fa)")
    parser.add_argument("--hg38_fa", required=True, help="Absolute host path for hg38 fasta (mounted as /ref/hg38.fa)")
    parser.add_argument("--chain_file", required=True, help="Absolute host path for liftover chain file (mounted as /ref/hg19ToHg38.over.chain)")
    parser.add_argument("--cores", type=int, default=1, help="Number of CPU cores for Snakemake")
    parser.add_argument("--engine", choices=["singularity", "docker"], default="singularity", help="Container engine to build command for")
    args = parser.parse_args()

    user_paths = {
        "vep_cache_hg19": os.path.abspath(args.vep_cache_hg19),
        "vep_cache_hg38": os.path.abspath(args.vep_cache_hg38),
        "hg19_fa": os.path.abspath(args.hg19_fa),
        "hg38_fa": os.path.abspath(args.hg38_fa),
        "chain_file": os.path.abspath(args.chain_file),
    }

    command = build_command(user_paths, args.cores, args.engine)
    print(command)
