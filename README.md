# Variant calling workflow using GATK4

For now, for simplicity, we assume that each sample is a 'platform unit'.

# Usage

1. Fill out/Edit the `bin/config.yaml` file.

2. Make samplesheet. You can go to the `bin/` folder and run `bash make_samples_template.sh`.
Rename the output file `samples_template.tsv` as `samples.tsv`.

3. A file named `grouped_contigs.tsv` with a column named 'name' and a second column named 'contigs' containing comma-delimited sequence names (chromosome/contig names) is needed. You may first load a R module by running `module load bbc2/R/alt/R-4.2.1-setR_LIBS_USER`, then run `Rscript group_contigs.R` to generate this file based on the reference sequence listed in `config.yaml`.

4. Submit job with `sbatch -p <parition name: e.g., long or bbc> bin/run_snakemake.sh`. 

# Workflow

![Workflow](./logs/rulegraph.png) 
