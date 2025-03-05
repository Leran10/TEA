"""
Snakemake pipeline for Transposable Element analysis
Tools included: STAR, RepeatMasker, TEtranscripts
"""

import os

# Configuration
configfile: "config.yaml"

# Define global variables
SAMPLES = config["samples"]
REFERENCE = config["reference_genome"]
TE_ANNOTATION = config["te_annotation"]
GTF = config["gene_annotation"]

# Get reference genome basename for outputs
ref_basename = os.path.basename(config["reference_genome"])

# Final output rule - this defines what files we want at the end
rule all:
    input:        
        # RepeatMasker outputs
        f"resources/repeatmasker/{ref_basename}.out",
        "resources/repeatmasker/TE.gtf",
        
        # Per-sample outputs
        expand("results/star/{sample}/Aligned.out.bam", sample=SAMPLES),
        expand("results/star/{sample}/Aligned.sortedByCoord.out.bam", sample=SAMPLES),
        expand("results/star/{sample}/Aligned.sortedByCoord.out.bam.bai", sample=SAMPLES),
        
        # Differential expression results
        "results/tetranscripts/differential",
        
        # QC report
        "results/multiqc/multiqc_report.html"

# STAR alignment
rule star_align:
    input:
        r1 = "data/fastq/{sample}_R1.fastq.gz",
        r2 = "data/fastq/{sample}_R2.fastq.gz",
        index = directory(config["star_index"])
    output:
        bam = "results/star/{sample}/Aligned.out.bam",
        log = "results/star/{sample}/Log.final.out"
    threads: 
        config["threads"]["star"]
    resources:
        mem_mb = config["memory"]["star"]
    shell:
        """
        STAR --runThreadN {threads} \
            --genomeDir {input.index} \
            --readFilesIn {input.r1} {input.r2} \
            --readFilesCommand zcat \
            --outFilterMultimapNmax 100 \
            --winAnchorMultimapNmax 100 \
            --outMultimapperOrder Random \
            --runRNGseed 777 \
            --outSAMmultNmax 1 \
            --outSAMtype BAM Unsorted \
            --outFilterType BySJout \
            --alignSJoverhangMin 8 \
            --alignSJDBoverhangMin 1 \
            --outFilterMismatchNmax 999 \
            --alignIntronMin 20 \
            --alignIntronMax 1000000 \
            --alignMatesGapMax 1000000 \
            --outFileNamePrefix results/star/{wildcards.sample}/ \
            --outReadsUnmapped Fastx \
            --quantMode GeneCounts
        """

# Sort BAM files
rule samtools_sort:
    input:
        "results/star/{sample}/Aligned.out.bam"
    output:
        "results/star/{sample}/Aligned.sortedByCoord.out.bam"
    threads: 4
    resources:
        mem_mb = 8000
    shell:
        """
        samtools sort -@ {threads} -m 2G -o {output} {input}
        """

# Index BAM files
rule samtools_index:
    input:
        "results/star/{sample}/Aligned.sortedByCoord.out.bam"
    output:
        "results/star/{sample}/Aligned.sortedByCoord.out.bam.bai"
    threads: 4
    shell:
        """
        samtools index {input}
        """

# Run RepeatMasker on the reference genome to generate masked regions and TE annotations
rule repeatmasker:
    input:
        genome = config["reference_genome"]
    output:
        masked = "resources/repeatmasker/{genome_basename}.masked",
        out = "resources/repeatmasker/{genome_basename}.out",
        gff = "resources/repeatmasker/{genome_basename}.out.gff",
        gtf = "resources/repeatmasker/TE.gtf"  # This is the file needed for TEtranscripts
    threads: 
        config["threads"]["repeatmasker"]
    resources:
        mem_mb = config["memory"]["repeatmasker"]
    params:
        species = config["repeatmasker_species"],
        genome_basename = lambda wildcards, input: os.path.basename(input.genome)
    shell:
        """
        # Create output directory
        mkdir -p resources/repeatmasker/
        
        # Run RepeatMasker with your specific parameters
        RepeatMasker --species {params.species} -pa {threads} -gff -dir resources/repeatmasker/ {input.genome}
        
        # Rename the output files to match expected names if needed
        if [ ! -f resources/repeatmasker/{params.genome_basename}.masked ]; then
            ln -s {params.genome_basename} resources/repeatmasker/{params.genome_basename}.masked
        fi
        
        # Convert RepeatMasker output to GTF format for TEtranscripts
        python scripts/prepare_te_gtf.py resources/repeatmasker/{params.genome_basename}.out resources/repeatmasker/TE.gtf
        """

# This rule has been removed as we only want to run RepeatMasker once on the reference genome

# TEtranscripts analysis for differential expression
rule tetranscripts:
    input:
        # Treatment samples BAMs
        treatment_bams = lambda wildcards: expand("results/star/{sample}/Aligned.sortedByCoord.out.bam", 
                                         sample=config["treatment_samples"]),
        # Control samples BAMs
        control_bams = lambda wildcards: expand("results/star/{sample}/Aligned.sortedByCoord.out.bam", 
                                       sample=config["control_samples"]),
        # Required BAI files for all samples
        treatment_bais = lambda wildcards: expand("results/star/{sample}/Aligned.sortedByCoord.out.bam.bai", 
                                         sample=config["treatment_samples"]),
        control_bais = lambda wildcards: expand("results/star/{sample}/Aligned.sortedByCoord.out.bam.bai", 
                                       sample=config["control_samples"]),
        # Direct paths to annotation files
        gene_gtf = config["gene_annotation"],
        te_gtf = config["te_annotation"]
    output:
        directory("results/tetranscripts/differential")
    threads: 
        config["threads"]["tetranscripts"]
    resources:
        mem_mb = config["memory"]["tetranscripts"]
    params:
        mode = config["tetranscripts_mode"]
    shell:
        """
        # Create output directory
        mkdir -p {output}
        
        # Run TEtranscripts for differential expression
        TEtranscripts \
            -t {input.treatment_bams} \
            -c {input.control_bams} \
            --GTF {input.gene_gtf} \
            --TE {input.te_gtf} \
            --sortByPos \
            --outdir {output} \
            --mode {params.mode} \
            --DESeq
        """

# Optional: MultiQC report generation
rule multiqc:
    input:
        expand("results/star/{sample}/Log.final.out", sample=SAMPLES)
    output:
        "results/multiqc/multiqc_report.html"
    shell:
        """
        multiqc results/ -o results/multiqc/
        """