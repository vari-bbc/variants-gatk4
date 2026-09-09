#!/bin/bash
#SBATCH --export=NONE
#SBATCH -J variants_workflow
#SBATCH -o logs/variants_workflow.o
#SBATCH -e logs/variants_workflow.e
#SBATCH --time 95:00:00
#SBATCH --mem=8G

#set -euo pipefail

cd ${SLURM_SUBMIT_DIR}

snakemake_module="bbc2/snakemake/snakemake-7.32.3"

module load $snakemake_module

# make logs dir if it does not exist already. Without this, logs/ is automatically generate only after the first run of the pipeline
logs_dir="logs/"
[[ -d $logs_dir ]] || mkdir -p $logs_dir

#snakemake --snakefile 'Snakefile' --dag | dot -Tpng > $logs_dir/dag.png
#snakemake --snakefile 'Snakefile' --filegraph | dot -Tpng > $logs_dir/filegraph.png
snakemake --snakefile 'Snakefile' --rulegraph | dot -Tpng > $logs_dir/rulegraph.png

echo "Start snakemake workflow." >&1                   
echo "Start snakemake workflow." >&2                   

snakemake \
-p \
--latency-wait 20 \
--snakefile 'Snakefile' \
--use-envmodules \
--jobs 100 \
--cluster "mkdir -p $logs_dir/{rule}; sbatch \
-p ${SLURM_JOB_PARTITION} \
--export=ALL \
--ntasks {threads} \
--mem={resources.mem_gb}G \
-t 48:00:00 \
-o $logs_dir/{rule}/{resources.log_prefix}.o \
-e $logs_dir/{rule}/{resources.log_prefix}.e" # SLURM hangs if output dir does not exist, so we create it before running sbatch on the snakemake jobs.

echo "snakemake workflow done." >&1                   
echo "snakemake workflow done." >&2       
