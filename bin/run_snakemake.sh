#PBS -l walltime=120:00:00
#PBS -l mem=8gb
#PBS -N variants_workflow
#PBS -o logs/variants_workflow.o
#PBS -e logs/variants_workflow.e

#set -euo pipefail

cd ${PBS_O_WORKDIR}

snakemake_module="bbc/snakemake/snakemake-6.15.0"

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
--cluster "ssh ${PBS_O_LOGNAME}@submit 'module load $snakemake_module; cd ${PBS_O_WORKDIR}; qsub \
-q ${PBS_O_QUEUE} \
-V \
-l nodes=1:ppn={threads} \
-l mem={resources.mem_gb}gb \
-l walltime=48:00:00 \
-o {log.stdout} \
-e {log.stderr}'"

echo "snakemake workflow done." >&1                   
echo "snakemake workflow done." >&2       
