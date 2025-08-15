import os
import argparse

def build_command(user_paths, cores=1, engine="singularity"):
    """
    Build a container run command (Singularity or Docker) for Snakemake with exact mounts.

    user_paths: dict with keys -
        data, results, vep_cache_hg19, vep_cache_hg38,
        hg19_fa, hg38_fa, chain_file, config_yaml (local path to config.yaml)
        tmp_home (only required for singularity)

    cores: int, number of CPU cores to use
    engine: 'singularity' or 'docker'
    """
    # Absolute host paths for required mounts
    mounts = {
        "data": os.path.abspath(user_paths["data"]),
        "results": os.path.abspath(user_paths["results"]),
        "vep_cache_hg19": os.path.abspath(user_paths["vep_cache_hg19"]),
        "vep_cache_hg38": os.path.abspath(user_paths["vep_cache_hg38"]),
        "hg19_fa": os.path.abspath(user_paths["hg19_fa"]),
        "hg19_fa_fai": os.path.abspath(user_paths["hg19_fa"]) + ".fai",
        "hg38_fa": os.path.abspath(user_paths["hg38_fa"]),
        "hg38_fa_fai": os.path.abspath(user_paths["hg38_fa"]) + ".fai",
        "chain_file": os.path.abspath(user_paths["chain_file"]),
        # Mount the local config.yaml into /app/snake/config.yaml inside container
        "config_yaml": os.path.abspath(user_paths["config_yaml"])
    }

    if engine == "singularity":
        tmp_home = os.path.abspath(user_paths["tmp_home"])
        binds = [
            f"{tmp_home}:/tmp_home",
            f"{mounts['data']}:/data",
            f"{mounts['results']}:/results",
            f"{mounts['vep_cache_hg19']}:/ref/VEP/cache_hg19",
            f"{mounts['vep_cache_hg38']}:/ref/VEP/cache_hg38",
            f"{mounts['hg19_fa']}:/ref/hg19.fa",
            f"{mounts['hg19_fa_fai']}:/ref/hg19.fa.fai",
            f"{mounts['hg38_fa']}:/ref/hg38.fa",
            f"{mounts['hg38_fa_fai']}:/ref/hg38.fa.fai",
            f"{mounts['chain_file']}:/ref/hg19ToHg38.over.chain",
            f"{mounts['config_yaml']}:/app/snake/config.yaml"
        ]
        binds_str = " -B ".join(binds)
        comment = ("# Singularity command to run Snakemake inside container with all required mounts "
                   "(including tmp_home). Please create the temporary home directory on host "
                   "to avoid permission issues during the run.")
        cmd = (
            f"{comment}\n"
            f"singularity exec --no-home -B {binds_str} "
            f"crossbuild.sif bash -c 'HOME=/tmp_home snakemake "
            f"--snakefile /app/snake/Snakefile --configfile /app/snake/config.yaml "
            f"-p --cores {cores}'"
        )
    elif engine == "docker":
        binds = [
            f"{mounts['data']}:/data:rw",
            f"{mounts['results']}:/results:rw",
            f"{mounts['vep_cache_hg19']}:/ref/VEP/cache_hg19:ro",
            f"{mounts['vep_cache_hg38']}:/ref/VEP/cache_hg38:ro",
            f"{mounts['hg19_fa']}:/ref/hg19.fa:ro",
            f"{mounts['hg19_fa_fai']}:/ref/hg19.fa.fai:ro",
            f"{mounts['hg38_fa']}:/ref/hg38.fa:ro",
            f"{mounts['hg38_fa_fai']}:/ref/hg38.fa.fai:ro",
            f"{mounts['chain_file']}:/ref/hg19ToHg38.over.chain:ro",
            f"{mounts['config_yaml']}:/app/snake/config.yaml:ro"
        ]
        binds_str = " ".join([f"-v {b}" for b in binds])
        container_image = "crossbuild:latest"  # Adjust docker image tag as needed
        comment = "# Docker command to run Snakemake inside container with all required mounts (tmp_home not mounted)"
        cmd = (
            f"{comment}\n"
            f"docker run --rm -it {binds_str} "
            f"{container_image} bash -c 'snakemake "
            f"--snakefile /app/snake/Snakefile --configfile /app/snake/config.yaml "
            f"-p --cores {cores}'"
        )
    else:
        raise ValueError("Engine must be 'singularity' or 'docker'")

    return cmd


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Build container run command for Snakemake with appropriate mounts")
    parser.add_argument("--tmp_home", help="Host path for tmp_home directory (required only for singularity)")
    parser.add_argument("--data", required=True, help="Host path for input data (mounted to /data)")
    parser.add_argument("--results", required=True, help="Host path for results directory (mounted to /results)")
    parser.add_argument("--vep_cache_hg19", required=True, help="Host path for VEP hg19 cache (mounted to /ref/VEP/cache_hg19)")
    parser.add_argument("--vep_cache_hg38", required=True, help="Host path for VEP hg38 cache (mounted to /ref/VEP/cache_hg38)")
    parser.add_argument("--hg19_fa", required=True, help="Host path for hg19 fasta (mounted to /ref/hg19.fa)")
    parser.add_argument("--hg38_fa", required=True, help="Host path for hg38 fasta (mounted to /ref/hg38.fa)")
    parser.add_argument("--chain_file", required=True, help="Host path for liftover chain file (mounted to /ref/hg19ToHg38.over.chain)")
    parser.add_argument("--config_yaml", required=True, help="Local host path to config.yaml (mounted to /app/snake/config.yaml)")
    parser.add_argument("--cores", type=int, default=1, help="Number of CPU cores for Snakemake")
    parser.add_argument("--engine", choices=["singularity", "docker"], default="singularity", help="Container engine to build command for")
    args = parser.parse_args()

    if args.engine == "singularity" and not args.tmp_home:
        parser.error("--tmp_home is required for singularity engine")

    user_paths = {
        "tmp_home": args.tmp_home,
        "data": args.data,
        "results": args.results,
        "vep_cache_hg19": args.vep_cache_hg19,
        "vep_cache_hg38": args.vep_cache_hg38,
        "hg19_fa": args.hg19_fa,
        "hg38_fa": args.hg38_fa,
        "chain_file": args.chain_file,
        "config_yaml": args.config_yaml,
    }

    command = build_command(user_paths, args.cores, args.engine)
    print(command)
