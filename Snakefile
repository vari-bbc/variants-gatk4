import pandas as pd
import numpy as np
import os
import re
import itertools
from snakemake.utils import validate, min_version
##### set minimum snakemake version #####
min_version("6.15.0")

##### Editable variables #####

configfile: "bin/config.yaml"

ref_fa = config['ref']['sequence']
ref_fai = config['ref']['fai']
gtf = config['ref']['annotation']
bwa_index = config['ref']['index']

mito_chrom=config['ref']['mito_chr']
fai_parsed = pd.read_table(ref_fai, names=['chr','len','offset','bases_per_line','bytes_per_line'])
if mito_chrom not in fai_parsed['chr'].values:
    raise Exception('{mito} not found in reference fai file.'.format(mito=mito_chrom))

chroms_no_mito = ' '.join(fai_parsed[fai_parsed['chr'] != mito_chrom]['chr'].values)

##### load config and sample sheets #####
samplesheet="bin/samples.tsv"
units = pd.read_table(samplesheet, dtype={"sample" : str})
units['se_or_pe'] = ["SE" if x else "PE" for x in units['fq2'].isnull()]

samples = units[["sample","se_or_pe"]].drop_duplicates()
if not samples['sample'].is_unique:
    raise Exception('A sample has more than one combination of and se_or_pe.')

contigs_file = "bin/grouped_contigs.tsv"

# read in file with col1 as contig group name and col2 as a commma separated list of contigs. This was written for use with GATK for splitting the variant calling steps by chromosome/contig. We group the unplaced contigs together since those tend to be small.
contig_groups = pd.read_table(contigs_file)

# if SE data, error out
if (units['fq2'].isnull().any()):
    raise Exception('SE not supported.')

snakemake_dir = os.getcwd() + "/"

# make a tmp directory for analyses
tmp_dir = os.path.join(snakemake_dir, "tmp")
if not os.path.exists(tmp_dir):
    os.mkdir(tmp_dir)

rule all:
    input:
        "analysis/multiqc/multiqc_report.html",
        #expand("analysis/final/06_filter_vcf/all.merged.filt.PASS.{var_type}.vcf.gz.vt_peek.txt", var_type=["SNP", "INDEL"]),
        "analysis/final/07_annot/snpEff/all.merged.filt.PASS.snpEff.vcf.gz",
        #"analysis/final/07_annot/Funcotator/all.merged.filt.PASS.Funcotator.vcf.gz",
        #"analysis/final/07_annot/Funcotator/target_genes_only/all.merged.filt.PASS.Funcotator.targetsOnly.filt.xlsx",
        expand("analysis/final/07_annot/Funcotator/target_genes_only/table/{sample}.merged.filt.PASS.Funcotator.targetsOnly.filt.xlsx", sample=samples['sample'].unique()),
        "analysis/final/07_annot/Funcotator/target_genes_only/table/all.merged.filt.PASS.Funcotator.targetsOnly.filt.xlsx",
        expand("analysis/final/07_annot/Funcotator/target_genes_only/custom_filters/{sample}.merged.filt.PASS.Funcotator.targetsOnly.filt.v3.xlsx", sample=samples['sample'].unique())


def get_orig_fastq(wildcards):
    if wildcards.read == "R1":
            fastq = expand("raw_data/{fq}", fq = units[units["sample"] == wildcards.sample]["fq1"].values)
    elif wildcards.read == "R2":
            fastq = expand("raw_data/{fq}", fq = units[units["sample"] == wildcards.sample]["fq2"].values)
    return fastq 

rule rename_fastqs:
    """
    Rename fastqs by biologically meaningful name. Concatenate different runs of same library.
    """
    input:
        get_orig_fastq
    output:
        "analysis/renamed_data/{sample}_{read}.fastq.gz"
    benchmark:
        "benchmarks/rename_fastqs/{sample}_{read}.txt"
    params:
    threads: 1
    resources:
        mem_gb=8,
        log_prefix=lambda wildcards: "_".join(wildcards)
    envmodules:
    shell:
        """
        if [ `printf '{input}' | wc -w` -gt 1 ]
        then
            cat {input} > {output}
        else
            ln -sr {input} {output}
        fi

        """

rule fastqc:
    """
    Run fastqc on raw_data/ files.
    """
    input:
        "analysis/renamed_data/{fq_pref}.fastq.gz"
    output:
        html="analysis/fastqc/{fq_pref}_fastqc.html",
        zip="analysis/fastqc/{fq_pref}_fastqc.zip"
    params:
        outdir="analysis/fastqc/"
    benchmark:
        "benchmarks/fastqc/{fq_pref}.txt"
    envmodules:
        config['modules']['fastqc']
    threads: 1
    resources:
        mem_gb = 32,
        log_prefix=lambda wildcards: "_".join(wildcards)
    shell:
        """
        fastqc --outdir {params.outdir} {input}
        """

rule fastq_screen:
    """
    Run fastq_screen to detect any contamination from other species or excessive rRNA.
    """
    input:
        "analysis/renamed_data/{fq_pref}.fastq.gz"
    output:
        html = "analysis/fastq_screen/{fq_pref}_screen.html",
        txt = "analysis/fastq_screen/{fq_pref}_screen.txt",
    params:
    benchmark:
        "benchmarks/fastq_screen/{fq_pref}.txt"
    envmodules:
        config['modules']['fastq_screen']
    threads: 8
    resources:
        mem_gb = 32,
        log_prefix=lambda wildcards: "_".join(wildcards)
    shell:
        """
        fastq_screen --threads {threads} --outdir analysis/fastq_screen/ {input}
        """

rule trim_galore_PE:
    """
    Run trim_galore on paired-end reads.
    """
    input:
        expand("analysis/renamed_data/{{sample}}_R{read}.fastq.gz", read=["1","2"])
    output:
        expand("analysis/trim_galore/{{sample}}_R{ext}", ext=["1_val_1.fq.gz","2_val_2.fq.gz"]),
        expand("analysis/trim_galore/{{sample}}_R1{ext}", ext=[".fastq.gz_trimming_report.txt","_val_1_fastqc.html"]),
        expand("analysis/trim_galore/{{sample}}_R2{ext}", ext=[".fastq.gz_trimming_report.txt","_val_2_fastqc.html"]),
    params:
    benchmark:
        "benchmarks/trim_galore/{sample}.txt"
    envmodules:
        config['modules']['trim_galore']
    threads: 4
    resources:
        mem_gb = 64,
        log_prefix=lambda wildcards: "_".join(wildcards)
    shell:
        """
        trim_galore --paired {input} --output_dir analysis/trim_galore/ --cores {threads} --fastqc
        """

rule trim_galore_SE:
    """
    Run trim_galore on single-end reads.
    """
    input:
        "analysis/renamed_data/{sample}_R1.fastq.gz"
    output:
        temp("analysis/trim_galore/{sample}_R1_trimmed.fq.gz"),
        "analysis/trim_galore/{sample}_R1.fastq.gz_trimming_report.txt",
        expand("analysis/trim_galore/{{sample}}_R1_trimmed_fastqc{ext}", ext=['.html','.zip']),
    params:
    benchmark:
        "benchmarks/trim_galore/{sample}.txt"
    envmodules:
        config['modules']['trim_galore']
    threads: 4
    resources:
        mem_gb = 64,
        log_prefix=lambda wildcards: "_".join(wildcards)
    shell:
        """
        trim_galore {input} --output_dir analysis/trim_galore/ --cores {threads} --fastqc
        """

rule bwamem:
    input:
        lambda wildcards: expand("analysis/trim_galore/{sample}_R{read}_val_{read}.fq.gz", read=["1","2"], sample=wildcards.sample) if samples[samples['sample']==wildcards.sample]['se_or_pe'].values=="PE" else "analysis/trim_galore/{sample}_R1_trimmed.fq.gz"
    output:
        outbam="analysis/bwamem/{sample}.bam",
        outbai="analysis/bwamem/{sample}.bam.bai",
        idxstat="analysis/bwamem/{sample}.bam.idxstat",
        samblaster_err="analysis/bwamem/{sample}.samblaster.e",
    benchmark:
        "benchmarks/bwamem/{sample}.txt"
    params:
        bwa_RG='"@RG\\tID:{sample}\\tLB:{sample}\\tPL:ILLUMINA\\tSM:{sample}"',
        bwa_idx=bwa_index,
        samblaster_params=lambda wildcards: "" if samples[samples['sample']==wildcards.sample]['se_or_pe'].values=="PE" else "--ignoreUnmated"
    threads: 16
    envmodules:
        config['modules']['bwa'],
        config['modules']['samblaster'],
        config['modules']['samtools'],
    resources:
        mem_gb=180,
        log_prefix=lambda wildcards: "_".join(wildcards)
    shell:
        """
        bwa mem \
        -t {threads} \
        -R {params.bwa_RG} \
        {params.bwa_idx} \
        {input} | \
        samblaster {params.samblaster_params} 2>{output.samblaster_err} | \
        samtools sort \
        -m 6G \
        -@ {threads} \
        -O "BAM" \
        -o {output.outbam} \
        -

        echo "END bwamem"
        echo "END bwamem" 1>&2

        samtools index -@ {threads} {output.outbam}

        echo "END indexing"
        echo "END indexing" 1>&2
        
        samtools idxstats {output.outbam} > {output.idxstat}

        echo "END idxstats"
        echo "END idxstats" 1>&2
 
        """

rule haplotypecaller:
    input:
        bam=lambda wildcards: "analysis/bwamem/{sample}.bam" if wildcards.round == "bqsr" else "analysis/final/00_BQSR/{sample}.bqsr.bam"
    output:
        "analysis/{round}/01_haplotypecaller/{sample}.{contig_group}.g.vcf.gz"
    benchmark:
        "benchmarks/{round}/01_haplotypecaller/{sample}.{contig_group}.txt"
    params:
        dbsnp=lambda wildcards: f'--dbsnp {config["ref"]["known_snps"]}' if config["ref"]["known_snps"] != "" else "",
        ref_fasta=config["ref"]["sequence"],
        contigs = lambda wildcards: "-L " + contig_groups[contig_groups.name == wildcards.contig_group]['contigs'].values[0].replace(",", " -L "),
        
    envmodules:
        config["modules"]["gatk"]
    threads: 4
    resources: 
        mem_gb = 80,
        log_prefix=lambda wildcards: "_".join(wildcards)
    shell:
        """
        gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" \
        HaplotypeCaller \
        -R {params.ref_fasta} \
        -I {input.bam} \
        -O {output} \
        -ERC GVCF \
        --native-pair-hmm-threads {threads} \
        {params.dbsnp} \
        {params.contigs}
        """

rule flagstat:
    """
    Run samtools flagstat.
    """
    input:
        "analysis/{align_dirname}/{bam_name}.bam"
    output:
        "analysis/{align_dirname}/flagstat/{bam_name}.flagstat"
    params:
    benchmark:
        "benchmarks/{align_dirname}/flagstat/{bam_name}.txt"
    envmodules:
        config['modules']['samtools']
    threads: 8
    resources: 
        mem_gb = 32,
        log_prefix=lambda wildcards: "_".join(wildcards)
    shell:
        """
        samtools flagstat -@ {threads} {input} > {output}  
        """

rule CollectInsertSizeMetrics:
    """
    Run Picard CollectInsertSizeMetrics.
    """
    input:
        "analysis/{align_dirname}/{bam_name}.bam"
    output:
        out="analysis/{align_dirname}/CollectInsertSizeMetrics/{bam_name}.insert_size_metrics.txt",
        hist="analysis/{align_dirname}/CollectInsertSizeMetrics/{bam_name}.insert_size_histogram.pdf"
    params:
        temp="./analysis/{align_dirname}/CollectInsertSizeMetrics/"
    benchmark:
        "benchmarks/{align_dirname}/CollectInsertSizeMetrics/{bam_name}.txt"
    envmodules:
        config['modules']['picard']
    threads: 4
    resources: 
        mem_gb = 64,
        log_prefix=lambda wildcards: "_".join(wildcards)
    shell:
        """
        java -Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir={params.temp} -jar $PICARD CollectInsertSizeMetrics I={input} O={output.out} H={output.hist}  
        """

rule CollectAlignmentSummaryMetrics:
    """
    Run Picard CollectAlignmentSummaryMetrics.
    """
    input:
        "analysis/{align_dirname}/{bam_name}.bam"
    output:
        out="analysis/{align_dirname}/CollectAlignmentSummaryMetrics/{bam_name}.aln_metrics.txt",
    params:
        temp="./analysis/{align_dirname}/CollectAlignmentSummaryMetrics/",
        reffasta=config["ref"]["sequence"]
    benchmark:
        "benchmarks/{align_dirname}/CollectAlignmentSummaryMetrics/{bam_name}.txt"
    envmodules:
        config['modules']['picard']
    threads: 4
    resources: 
        mem_gb = 64,
        log_prefix=lambda wildcards: "_".join(wildcards)
    shell:
        """
        java -Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir={params.temp} -jar $PICARD CollectAlignmentSummaryMetrics I={input} O={output.out} R={params.reffasta}  
        """

rule combinevar:
    input:
        lambda wildcards: expand("analysis/{{round}}/01_haplotypecaller/{sample}.{contig_group}.g.vcf.gz", sample=samples['sample'].unique(), contig_group=wildcards.contig_group)

    output:
        touch=touch("analysis/{round}/02_combinevar/{contig_group}.done"),
        genomicsdb=directory("analysis/{round}/02_combinevar/{contig_group}.genomicsdb"),
    benchmark:
        "benchmarks/{round}/02_combinevar/{contig_group}.txt"
    params:
        sample_gvcfs = lambda wildcards, input: list(map("-V {}".format, input)),
        contigs = lambda wildcards: "-L " + contig_groups[contig_groups.name == wildcards.contig_group]['contigs'].values[0].replace(",", " -L "),
    envmodules:
        config["modules"]["gatk"]
    threads: 4
    resources: 
        mem_gb = 80,
        log_prefix=lambda wildcards: "_".join(wildcards)
    shell:
        """
        gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" \
        GenomicsDBImport \
        {params.sample_gvcfs} \
        --genomicsdb-workspace-path {output.genomicsdb} \
        {params.contigs}
        """

rule jointgeno:
    input:
        "analysis/{round}/02_combinevar/{contig_group}.done"
    output:
        vcf="analysis/{round}/03_jointgeno/all.{contig_group}.vcf.gz",
    benchmark:
        "benchmarks/{round}/03_jointgeno/all.{contig_group}.txt"
    params:
        ref_fasta=config["ref"]["sequence"],
        genomicsdb="analysis/{round}/02_combinevar/{contig_group}.genomicsdb"
    envmodules:
        config["modules"]["gatk"]
    threads: 4
    resources: 
        mem_gb = 80,
        log_prefix=lambda wildcards: "_".join(wildcards)
    shell:
        """
        gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" \
        GenotypeGVCFs \
        -R {params.ref_fasta} \
        -V gendb://{params.genomicsdb} \
        -O {output.vcf}
        """

rule sortVCF:
    """
    Sort the output VCFs from joint genotyping. Merging errors out sometimes if we do not do this step.
    """
    input:
        vcf="analysis/{round}/03_jointgeno/all.{contig_group}.vcf.gz",
    output:
        sorted_vcf="analysis/{round}/04_sortvcf/all.{contig_group}.sort.vcf.gz"
    benchmark:
        "benchmarks/{round}/04_sortvcf/all.{contig_group}.txt"
    params:
        dictionary=config['ref']['dict'],
    envmodules:
        config["modules"]["gatk"]
    threads: 4
    resources: 
        mem_gb = 80,
        log_prefix=lambda wildcards: "_".join(wildcards)
    shell:
        """
        gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" \
        SortVcf \
        -I {input.vcf} \
        -O {output.sorted_vcf} \
        -SD {params.dictionary} 
        """

rule merge_vcf:
    """
    Merge the contig group VCFs into one unified VCF.
    """
    input:
        expand("analysis/{{round}}/04_sortvcf/all.{contig_grp}.sort.vcf.gz", contig_grp=contig_groups.name)
    output:
        raw="analysis/{round}/05_merge_vcf/all.merged.vcf.gz",
        raw_tbi="analysis/{round}/05_merge_vcf/all.merged.vcf.gz.tbi",
        vt_peek_raw="analysis/{round}/05_merge_vcf/all.merged.vcf.gz.vt_peek.txt",
    benchmark:
        "benchmarks/{round}/05_merge_vcf/benchmark.txt"
    params:
        ref_fasta=config["ref"]["sequence"],
        dictionary=config['ref']['dict'],
        in_vcfs = lambda wildcards, input: ' '.join(['--INPUT ' + vcf for vcf in input]) 
    envmodules:
        config["modules"]["gatk"],
        config["modules"]["vt"],
    threads: 4
    resources: 
        mem_gb = 80,
        log_prefix=lambda wildcards: "_".join(wildcards)
    shell:
        """
        gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" \
        MergeVcfs \
        {params.in_vcfs} \
        --SEQUENCE_DICTIONARY {params.dictionary} \
        --OUTPUT {output.raw} 
        
        vt peek -r {params.ref_fasta} {output.raw} 2> {output.vt_peek_raw}
        """

def get_filt_params (wildcards):
    # parameters adapted from https://gencore.bio.nyu.edu/variant-calling-pipeline-gatk4/
    if (wildcards.var_type == "SNP"):
        return """--cluster-window-size 10 --cluster-size 3 \
        --filter-name 'QD' --filter-expression 'QD < 2.0' \
        --filter-name 'FS' --filter-expression 'FS > 60.0' \
        --filter-name 'MQ_filter' --filter-expression 'MQ < 40.0' \
        --filter-name 'SOR_filter' --filter-expression 'SOR > 4.0' \
        --filter-name 'MQRankSum_filter' --filter-expression 'MQRankSum < -12.5' \
        --filter-name 'ReadPosRankSum_filter' --filter-expression 'ReadPosRankSum < -8.0' \
        --genotype-filter-name 'GQ' --genotype-filter-expression 'GQ < 20'"""
    elif (wildcards.var_type == "INDEL"):
        return """--filter-name 'QD' --filter-expression 'QD < 2.0' \
        --filter-name 'FS' --filter-expression 'FS > 200.0' \
        --filter-name 'SOR_filter' --filter-expression 'SOR > 10.0' \
        --genotype-filter-name 'GQ' --genotype-filter-expression 'GQ < 20'"""
    else:
        raise Exception("var_type wildcard must be either SNP or INDEL.")

rule filter_vcf:
    """
    Do quality filters. Use different paramters depending on SNPs verus indels ('SNP' or 'INDEL').
    """
    input:
        "analysis/{round}/05_merge_vcf/all.merged.vcf.gz"
    output:
        raw="analysis/{round}/06_filter_vcf/all.merged.{var_type}.vcf.gz",
        filt="analysis/{round}/06_filter_vcf/all.merged.filt.{var_type}.vcf.gz",
        pass_only="analysis/{round}/06_filter_vcf/all.merged.filt.PASS.{var_type}.vcf.gz",
        vt_peek_pass="analysis/{round}/06_filter_vcf/all.merged.filt.PASS.{var_type}.vcf.gz.vt_peek.txt"
    benchmark:
        "benchmarks/{round}/06_filter_vcf/{var_type}.txt"
    params:
        ref_fasta=config["ref"]["sequence"],
        filt_params=get_filt_params
    envmodules:
        config["modules"]["gatk"],
        config["modules"]["vt"]
    threads: 4
    resources: 
        mem_gb = 80,
        log_prefix=lambda wildcards: "_".join(wildcards)
    shell:
        """
        gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" \
        SelectVariants \
        -R {params.ref_fasta} \
        -V {input} \
        --select-type-to-include {wildcards.var_type} \
        -O {output.raw}

        echo "SelectVariants 1 done."
        echo "SelectVariants 1 done." >&2
        
        gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" \
        VariantFiltration \
        --R {params.ref_fasta} \
        --V {output.raw} \
        {params.filt_params} \
        -O {output.filt} 
        
        echo "VariantFiltration done."
        echo "VariantFiltration done." >&2

        gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" \
        SelectVariants \
        -R {params.ref_fasta} \
        -V {output.filt} \
        --exclude-filtered \
        --set-filtered-gt-to-nocall \
        -O {output.pass_only} 
        
        echo "SelectVariants 2 done."
        echo "SelectVariants 2 done." >&2

        vt peek -r {params.ref_fasta} {output.pass_only} 2> {output.vt_peek_pass}
        """

rule BQSR:
    """
    Base quality score recalibrator
    """
    input:
        bam="analysis/bwamem/{sample}.bam",
        known_snps_vcf=config["ref"]["known_snps"] if config["ref"]["known_snps"] else "analysis/bqsr/06_filter_vcf/all.merged.filt.PASS.SNP.vcf.gz",
        known_indels_vcf=config["ref"]["known_indels"] if config["ref"]["known_indels"] else "analysis/bqsr/06_filter_vcf/all.merged.filt.PASS.INDEL.vcf.gz",
    output:
        recal_table="analysis/final/00_BQSR/{sample}.bqsr.table",
        recal_bam="analysis/final/00_BQSR/{sample}.bqsr.bam"
    benchmark:
        "benchmarks/final/00_BQSR/{sample}.txt"
    params:
        ref_fasta=config["ref"]["sequence"],
    envmodules:
        config["modules"]["gatk"],
    threads: 4
    resources: 
        mem_gb = 80,
        log_prefix=lambda wildcards: "_".join(wildcards)
    shell:
        """
        gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" \
        BaseRecalibrator \
        -R {params.ref_fasta} \
        -I {input.bam} \
        --known-sites {input.known_snps_vcf} \
        --known-sites {input.known_indels_vcf} \
        -O {output.recal_table} 
 
        echo "BaseRecalibrator done."
        echo "BaseRecalibrator done." >&2
        
        gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" \
        ApplyBQSR \
        -R {params.ref_fasta} \
        -I {input.bam} \
        -bqsr {output.recal_table} \
        -O {output.recal_bam}

        echo "ApplyBQSR done."
        echo "ApplyBQSR done." >&2

        """

rule qualimap:
    """
    Run qualimap on unfiltered alignments.
    """
    input:
        "analysis/bwamem/{sample}.bam",
    output:
        touch("analysis/qualimap/{sample}/done")
    benchmark:
        "benchmarks/qualimap/{sample}.txt"
    envmodules:
        config['modules']['qualimap']
    params:
    resources:
        mem_gb=100,
        log_prefix=lambda wildcards: "_".join(wildcards)
    threads: 8
    shell:
        """
        qualimap bamqc -bam {input} --java-mem-size={resources.mem_gb}G --paint-chromosome-limits -outdir analysis/qualimap/{wildcards.sample} -nt {threads}

        """

rule concat_variants:
    """
    Concatenate SNPs and indels.
    """
    input:
        vcfs=expand("analysis/final/06_filter_vcf/all.merged.filt.PASS.{var_type}.vcf.gz", var_type=["SNP", "INDEL"]),
    output:
        merged_vcf="analysis/final/07_annot/all.merged.filt.PASS.vcf.gz",
    benchmark:
        "benchmarks/final/07_annot/concat_variants/bench.txt"
    envmodules:
        config['modules']['bcftools']
    params:
    resources:
        mem_gb=100,
        log_prefix="concat_variants"
    threads: 8
    shell:
        """
        # merge
        bcftools concat --threads {threads} --allow-overlaps -O z -o {output.merged_vcf} {input.vcfs}
        bcftools index --tbi --threads {threads} {output.merged_vcf}
        
        """

rule normalize_variants:
    """
    Decompose and normalize variants.
    """
    input:
        genome_fa=config['ref']['sequence'],
        merged_vcf="analysis/final/07_annot/all.merged.filt.PASS.vcf.gz"
    output:
        merged_vcf="analysis/final/07_annot/all.merged.filt.PASS.decomp_and_norm.vcf.gz",
    benchmark:
        "benchmarks/final/07_annot/normalize_variants/bench.txt"
    envmodules:
        config['modules']['bcftools']
    params:
    resources:
        mem_gb=100,
        log_prefix="normalize_variants"
    threads: 8
    shell:
        """
        # see https://github.com/brentp/slivar/wiki/decomposing-and-subsetting-vcfs
        bcftools norm --threads {threads} --multiallelics - -w 10000 --fasta-ref {input.genome_fa} -O z -o {output.merged_vcf} {input.merged_vcf}
        bcftools index --tbi --threads {threads} {output.merged_vcf}
        
        """

rule snpEff:
    """
    Run snpEff.
    """
    input:
        merged_vcf="analysis/final/07_annot/all.merged.filt.PASS.vcf.gz"
    output:
        html="analysis/final/07_annot/snpEff/snpEff.html",
        genes="analysis/final/07_annot/snpEff/snpEff.genes.txt",
        annot_vcf="analysis/final/07_annot/snpEff/all.merged.filt.PASS.snpEff.vcf.gz"
    benchmark:
        "benchmarks/final/07_annot/snpEff/bench.txt"
    envmodules:
        config['modules']['snpEff'],
        config['modules']['bcftools']
    params:
        snpEff_dataDir=f'-dataDir {config["ref"]["snpEff_dataDir"]}' if config['ref']['snpEff_dataDir'] else "",
        snpEff_genome_id=config['ref']['snpEff_genome_ID']
    resources:
        mem_gb=100,
        log_prefix="snpEff"
    threads: 8
    shell:
        """
        
        # annotate
        java -Xms8g -Xmx80g -Djava.io.tmpdir=./tmp -jar $SNPEFF/snpEff.jar {params.snpEff_dataDir} -v -stats {output.html} {params.snpEff_genome_id} {input.merged_vcf} | bcftools view --threads {threads} -O z -o {output.annot_vcf} -
        bcftools index --tbi --threads {threads} {output.annot_vcf}

        """


rule Funcotator:
    """
    Run Funcotator.
    """
    input:
        merged_vcf="analysis/final/07_annot/all.merged.filt.PASS.decomp_and_norm.vcf.gz",
        ref_fa=config['ref']['sequence']
    output:
        annot_vcf="analysis/final/07_annot/Funcotator/all.merged.filt.PASS.Funcotator.vcf.gz"
    benchmark:
        "benchmarks/final/07_annot/Funcotator/bench.txt"
    envmodules:
        config["modules"]["gatk"],
        config['modules']['bcftools']
    params:
        data_sources=config['Funcotator']['data_source'],
        ref_ver=config['Funcotator']['ref_ver']
    resources:
        mem_gb=100,
        log_prefix="Funcotator"
    threads: 8
    shell:
        """
        gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" \
                Funcotator \
                --data-sources-path {params.data_sources} \
                --output {output.annot_vcf} \
                --output-file-format "VCF" \
                --ref-version {params.ref_ver} \
                --reference {input.ref_fa} \
                --variant {input.merged_vcf}

        """



#rule slivar:
#    """
#    Run slivar.
#    """
#    input:
#        merged_vcf="analysis/final/07_annot/all.merged.filt.PASS.decomp_and_norm.vcf.gz",
#        genome_fa=config['ref']['sequence'],
#        ensembl_gff3=config['slivar']['ensembl_gff3'],
#    output:
#        bcftools_csq="analysis/final/07_annot/slivar/all.merged.filt.PASS.csq.vcf.gz",
#        slivar="analysis/final/07_annot/slivar/all.merged.filt.PASS.csq.slivar.vcf.gz"
#    benchmark:
#        "benchmarks/final/07_annot/slivar/bench.txt"
#    envmodules:
#        config['modules']['slivar'],
#        config['modules']['bcftools']
#    params:
#    resources:
#        mem_gb=100,
#        log_prefix="slivar"
#    threads: 8
#    shell:
#        """
#        # see https://github.com/brentp/slivar/wiki/rare-disease
#        bcftools csq --samples - --ncsq 40 --gff-annot {input.ensembl_gff3} --local-csq --fasta-ref {input.genome_fa} -O z -o {output.bcftools_csq}
#        bcftools index --tbi --threads {threads} {output.bcftools_csq}
#
#
#        """
#


rule make_BED_from_genes:
    input:
        gtf=config['ref']['annotation'],
        target_genes=config['Funcotator']['genes']
    output:
        bed="analysis/final/07_annot/Funcotator/target_genes_only/make_BED_from_genes/target_genes.bed"
    benchmark:
        "benchmarks/final/07_annot/make_BED_from_gene/bench.txt"
    params:
    threads: 1
    envmodules:
        config['modules']['R']
    resources:
        mem_gb=64,
        log_prefix="make_BED_from_genes"
    shell:
        """
        Rscript --vanilla -e "library(rtracklayer); genes <- readLines('{input.target_genes}'); gtf <- import('{input.gtf}'); score(gtf) <- 0; outs <- gtf[gtf\$gene_name %in% genes & gtf\$type=='gene']; names(outs) <- outs\$gene_name; export(outs, '{output.bed}')"
        """

rule filter_variants_to_targets:
    """
    Add FORMAT/VAF tag, filter VCF for just target regions, implement final filters for either a specific sample or all of them.
    """
    input:
        vcf="analysis/final/07_annot/Funcotator/all.merged.filt.PASS.Funcotator.vcf.gz",
        regions="analysis/final/07_annot/Funcotator/target_genes_only/make_BED_from_genes/target_genes.bed"
    output:
        vcf="analysis/final/07_annot/Funcotator/target_genes_only/vcf/{sample}.merged.filt.PASS.Funcotator.targetsOnly.filt.vcf.gz",
        funcotator_columns=temp("analysis/final/07_annot/Funcotator/target_genes_only/vcf/{sample}.funcotator_column_names.txt"),
    benchmark:
        "benchmarks/final/07_annot/filter_variants_to_targets/{sample}.txt"
    envmodules:
        config['modules']['bcftools'],
    params:
        sample=lambda wildcards: "" if wildcards.sample == "all" else "-s " + wildcards.sample
    resources:
        mem_gb=100,
        log_prefix=lambda wildcards: "_".join(wildcards)
    threads: 8
    shell:
        """
        # Filter for just regions of interest, add VAF tag, do final filters
        bcftools view --threads {threads} -R {input.regions} -O u {params.sample} {input.vcf} | \
        bcftools +fill-tags - -O u -- -t FORMAT/VAF | \
        bcftools view --threads {threads} -O z -o {output.vcf} -i 'FMT/DP>9 & FMT/VAF>0.3' -

        bcftools index --tbi --threads {threads} {output.vcf}
 
        # grep the FUNCOTATION column names from the VCF header
        zcat {output.vcf} | grep '##INFO=<ID=FUNCOTATION' | grep -Po '([^"\s]+\|[^"\s]+)+' > {output.funcotator_columns}

        """

rule make_variant_table:
    """
    Final formatting.
    """
    input:
        vcf="analysis/final/07_annot/Funcotator/target_genes_only/vcf/{sample}.merged.filt.PASS.Funcotator.targetsOnly.filt.vcf.gz",
        funcotator_columns="analysis/final/07_annot/Funcotator/target_genes_only/vcf/{sample}.funcotator_column_names.txt",
    output:
        table="analysis/final/07_annot/Funcotator/target_genes_only/table/{sample}.merged.filt.PASS.Funcotator.targetsOnly.filt.tsv",
        final_table="analysis/final/07_annot/Funcotator/target_genes_only/table/{sample}.merged.filt.PASS.Funcotator.targetsOnly.filt.xlsx"
    benchmark:
        "benchmarks/final/07_annot/make_variant_table/{sample}.txt"
    envmodules:
        config['modules']['R'],
        config['modules']['gatk']
    params:
    resources:
        mem_gb=100,
        log_prefix=lambda wildcards: "_".join(wildcards)
    threads: 8
    shell:
        """
        # convert vcf to table format
        gatk VariantsToTable -V {input.vcf} \
        -F CHROM -F POS -F TYPE -F REF -F ALT -F FUNCOTATION \
        -GF GT -GF AD -GF VAF \
        -O {output.table}

        # split the FUNCOTATION column
        Rscript --vanilla -e 'library(dplyr); library(readr); library(stringr); library(tidyr); library(openxlsx);\
        colnms <- str_split_1(readLines("{input.funcotator_columns}"), "\\\|");\
        tab <- read_tsv("{output.table}") %>% \
        dplyr::select(!one_of("FUNCOTATION"), FUNCOTATION) %>% \
        dplyr::mutate(FUNCOTATION=str_remove_all(FUNCOTATION, "^\\\[|\\\]$")) %>% \
        separate(col=FUNCOTATION, into=colnms, sep="\\\|", remove=TRUE);\
        write.xlsx(tab, "{output.final_table}", overwrite=TRUE)' 

        """

rule custom_filters_and_formating_on_final_table:
    """
    Custom filters and formating on final table.
    """
    input:
        "analysis/final/07_annot/Funcotator/target_genes_only/table/{sample}.merged.filt.PASS.Funcotator.targetsOnly.filt.xlsx"
    output:
        v1="analysis/final/07_annot/Funcotator/target_genes_only/custom_filters/{sample}.merged.filt.PASS.Funcotator.targetsOnly.filt.v1.xlsx",
        v2="analysis/final/07_annot/Funcotator/target_genes_only/custom_filters/{sample}.merged.filt.PASS.Funcotator.targetsOnly.filt.v2.xlsx",
        v3="analysis/final/07_annot/Funcotator/target_genes_only/custom_filters/{sample}.merged.filt.PASS.Funcotator.targetsOnly.filt.v3.xlsx",
    benchmark:
        "benchmarks/final/07_annot/custom_filters_and_formating_on_final_table/{sample}.txt"
    envmodules:
        config['modules']['R'],
    params:
    resources:
        mem_gb=100,
        log_prefix=lambda wildcards: "_".join(wildcards)
    threads: 1
    script:
        "scripts/custom_filters_on_final_table.R"

rule multiqc:
    input:
        expand("analysis/fastqc/{sample.sample}_R1_fastqc.html", sample=samples.itertuples()),
        expand("analysis/fastqc/{sample.sample}_R2_fastqc.html", sample=samples[samples['se_or_pe']=="PE"].itertuples()),
        expand("analysis/trim_galore/{sample.sample}_R{read}_val_{read}_fastqc.html", sample=samples[samples['se_or_pe']=="PE"].itertuples(), read=["1","2"]),
        expand("analysis/trim_galore/{sample.sample}_R1_trimmed_fastqc.html", sample=samples[samples['se_or_pe']=="SE"].itertuples()),
        expand("analysis/bwamem/flagstat/{sample.sample}.flagstat", sample=samples.itertuples()),
        expand("analysis/bwamem/CollectInsertSizeMetrics/{sample.sample}.insert_size_metrics.txt", sample=samples[samples['se_or_pe']=="PE"].itertuples()),
        expand("analysis/bwamem/{sample.sample}.bam.idxstat", sample=samples.itertuples()),
        expand("analysis/bwamem/CollectAlignmentSummaryMetrics/{sample.sample}.aln_metrics.txt", sample=samples.itertuples()),
        expand("analysis/bwamem/{sample.sample}.samblaster.e", sample=samples.itertuples()),
        expand("analysis/fastq_screen/{sample.sample}_R1_screen.html", sample=samples.itertuples()),
        expand("analysis/fastq_screen/{sample.sample}_R2_screen.html", sample=samples[samples['se_or_pe']=="PE"].itertuples()),
        expand("analysis/qualimap/{sample.sample}/done", sample=samples.itertuples()),
    output:
        "analysis/multiqc/multiqc_report.html",
    benchmark:
        "benchmarks/multiqc/multiqc.txt"
    params:
        workdir="analysis/multiqc/",
        dirs=["analysis/trim_galore/",
        "analysis/fastqc/",
        "analysis/qualimap/",
        "analysis/bwamem/flagstat/",
        "analysis/bwamem/",
        "analysis/fastq_screen/",
        "analysis/bwamem/CollectAlignmentSummaryMetrics/"],
        outfile="multiqc_report"
    envmodules:
        config['modules']['multiqc']
    threads: 4
    resources:
        mem_gb=100,
        log_prefix="MultiQC"
    shell:
        """
        multiqc \
        --force \
        --outdir {params.workdir} \
        --filename {params.outfile} \
        {params.dirs} 

        """

