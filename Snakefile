configfile: "config.yaml"
import os

# Define samples from config
SAMPLES = config["samples"]
TREATMENT_SAMPLES = config["treatment_samples"]
CONTROL_SAMPLES = config["control_samples"]

# Get reference genome and GTF paths from config
GENOME_PATH = config["genome_fasta"]
GTF_PATH = config["gtf_file"]
TE_GTF_PATH = config["te_gtf_file"]  # New parameter for the fixed TE GTF file

rule all:
    input:
        # Final output: TEtranscripts differential expression analysis
        "results/tetranscripts/differential_expression.tsv",
        # Ensure all BAM files are properly sorted and indexed
        expand("results/star/{sample}/Aligned.sortedByCoord.out.bam.bai", sample=SAMPLES)

# No download rules - user provides their own genome and GTF files

rule star_index:
    input:
        genome = GENOME_PATH,
        gtf = GTF_PATH
    output:
        directory("results/star_index")
    threads: 
        config["threads"]["star_index"]
    shell:
        """
        mkdir -p {output}
        STAR --runMode genomeGenerate \
             --genomeDir {output} \
             --genomeFastaFiles {input.genome} \
             --sjdbGTFfile {input.gtf} \
             --runThreadN {threads}
        """

rule star_align:
    input:
        index = "results/star_index",
        r1 = "data/raw_fastq/{sample}_R1.fastq.gz",
        r2 = "data/raw_fastq/{sample}_R2.fastq.gz"
    output:
        bam = "results/star/{sample}/Aligned.out.bam",
        log = "results/star/{sample}/Log.final.out"
    params:
        prefix = "results/star/{sample}/"
    threads: 
        config["threads"]["star_align"]
    shell:
        """
        mkdir -p {params.prefix}
        STAR --genomeDir {input.index} \
             --readFilesIn {input.r1} {input.r2} \
             --readFilesCommand zcat \
             --outFileNamePrefix {params.prefix} \
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
             --runThreadN {threads}
        """

# RepeatMasker and TE GTF creation steps removed - using pre-existing TE GTF file

rule sort_bam:
    input:
        "results/star/{sample}/Aligned.out.bam"
    output:
        bam = "results/star/{sample}/Aligned.sortedByCoord.out.bam",
        bai = "results/star/{sample}/Aligned.sortedByCoord.out.bam.bai"
    threads:
        config["threads"]["samtools"]
    shell:
        """
        samtools sort -@ {threads} -o {output.bam} {input}
        samtools index {output.bam}
        """

rule tetranscripts:
    input:
        treatment_bams = expand("results/star/{sample}/Aligned.sortedByCoord.out.bam", sample=TREATMENT_SAMPLES),
        control_bams = expand("results/star/{sample}/Aligned.sortedByCoord.out.bam", sample=CONTROL_SAMPLES),
        gtf = GTF_PATH,
        te_gtf = TE_GTF_PATH  # Use the fixed TE GTF path from config
    output:
        diff_file = "results/tetranscripts/differential_expression.tsv"
    params:
        outdir = "results/tetranscripts"
    threads:
        config["threads"]["tetx"]
    shell:
        """
        mkdir -p {params.outdir}
        
        # Ensure all input files exist
        for file in {input.treatment_bams} {input.control_bams}; do
            if [ ! -f "$file" ]; then
                echo "ERROR: Input BAM file $file does not exist!"
                exit 1
            fi
        done
        if [ ! -f "{input.gtf}" ]; then
            echo "ERROR: GTF file {input.gtf} does not exist!"
            exit 1
        fi
        if [ ! -f "{input.te_gtf}" ]; then
            echo "ERROR: TE GTF file {input.te_gtf} does not exist!"
            exit 1
        fi
        
        TEtranscripts \
            -t {input.treatment_bams} \
            -c {input.control_bams} \
            --GTF {input.gtf} \
            --TE {input.te_gtf} \
            --sortByPos \
            --outdir {params.outdir} \
            --mode multi \
            --DESeq \
            --project differential_expression
        """

# Quality control rules removed