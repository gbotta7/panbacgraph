source /homes2/gianfranco/miniforge3/etc/profile.d/conda.sh
conda activate snakemake

# snakemake --use-conda --cores 1 --conda-create-envs-only # build conda envs outside singularity
snakemake --use-conda --cores 24 --rerun-incomplete 2>&1 | tee snakemake_$(date +%Y%m%d_%H%M%S).log