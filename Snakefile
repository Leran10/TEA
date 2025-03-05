"""
Snakemake pipeline for CRISPR screen analysis.
Tools used: BBDuk.sh (BBTools) for read trimming and MAGeCK for sgRNA enrichment/depletion analysis
"""

# Import required modules
import os
import yaml
from os.path import join, exists
import glob

# Load configuration
configfile: "config.yaml"

# Load library definitions
LIBRARIES_CONFIG = "libraries/libraries.yaml"
if exists(LIBRARIES_CONFIG):
    with open(LIBRARIES_CONFIG, 'r') as f:
        LIBRARIES = yaml.safe_load(f)
else:
    LIBRARIES = {}

# Define variables
SAMPLES = config["samples"]
COMPARISONS = config["comparisons"]
OUTPUT_DIR = config["output_dir"]

# Determine the library file and parameters
if "library_name" in config and config["library_name"] in LIBRARIES:
    # Using a predefined library by name
    library_name = config["library_name"]
    library_info = LIBRARIES[library_name]
    LIBRARY = library_info["path"]
    
    # Use library-specific parameters if not explicitly defined in config
    if "vector_sequence" not in config:
        VECTOR = library_info.get("vector_sequence", "GAAACACCG")
    else:
        VECTOR = config["vector_sequence"]
        
    if "backbone_sequence" not in config:
        BACKBONE = library_info.get("backbone_sequence", "GTTTTAGA")
    else:
        BACKBONE = config["backbone_sequence"]
        
    if "bbduk_min_length" not in config:
        MIN_LENGTH = library_info.get("min_length", 20)
    else:
        MIN_LENGTH = config["bbduk_min_length"]
        
    if "bbduk_max_length" not in config:
        MAX_LENGTH = library_info.get("max_length", 20)
    else:
        MAX_LENGTH = config["bbduk_max_length"]
        
    print(f"Using predefined library: {library_name}")
    print(f"Description: {library_info.get('description', 'No description available')}")
    
elif "library_file" in config:
    # Direct file path provided
    LIBRARY = config["library_file"]
    
    # Set default values for vector and backbone sequences if not specified
    VECTOR = config.get("vector_sequence", "GAAACACCG")
    BACKBONE = config.get("backbone_sequence", "GTTTTAGA")
    MIN_LENGTH = config.get("bbduk_min_length", 20)
    MAX_LENGTH = config.get("bbduk_max_length", 20)
    
    print(f"Using library file: {LIBRARY}")
    
else:
    raise ValueError("Neither 'library_name' nor 'library_file' specified in config.yaml")

# Define output directories
BBDUK_DIR = join(OUTPUT_DIR, "bbduk_output")
MAGECK_DIR = join(OUTPUT_DIR, "mageck_output")

# Ensure output directories exist
for d in [BBDUK_DIR, join(BBDUK_DIR, "logs"), MAGECK_DIR]:
    os.makedirs(d, exist_ok=True)

# All expected output files
rule all:
    input:
        # BBDuk trimmed files
        expand(join(BBDUK_DIR, "{sample}_rtrim.fastq"), sample=SAMPLES),
        # BBDuk summary report
        join(BBDUK_DIR, "bbduk_summary_report.html"),
        join(BBDUK_DIR, "bbduk_summary_stats.tsv"),
        # MAGeCK count file
        join(MAGECK_DIR, "mageck_counts.count.txt"),
        # MAGeCK test results for each comparison
        expand(join(MAGECK_DIR, "{comparison}.gene_summary.txt"), comparison=COMPARISONS)

# BBDuk trimming rule
rule bbduk_trim:
    input:
        fastq = "{sample}.fastq"
    output:
        lfilter = join(BBDUK_DIR, "{sample}_lfilter.fastq"),
        ltrim = join(BBDUK_DIR, "{sample}_ltrim.fastq"),
        rfilter = join(BBDUK_DIR, "{sample}_rfilter.fastq"),
        rtrim = join(BBDUK_DIR, "{sample}_rtrim.fastq"),
        rtrim_miss = join(BBDUK_DIR, "{sample}_rtrim_miss.fastq")
    params:
        # Use parameters from config if available, otherwise use defaults
        vector = VECTOR,
        backbone = BACKBONE,
        k = config.get("bbduk_kmer_length", 8),
        rcomp = config.get("bbduk_reverse_complement", "f"),
        mm = config.get("bbduk_mismatch", "f"),
        hdist = config.get("bbduk_hamming_distance", 0),
        min_length = MIN_LENGTH,
        max_length = MAX_LENGTH,
        log_dir = join(BBDUK_DIR, "logs")
    threads: config["threads"]["bbduk"]
    resources:
        mem_mb = config["memory"]["bbduk"]
    log:
        lfilter = join(BBDUK_DIR, "logs", "{sample}_lfilter_log.txt"),
        ltrim = join(BBDUK_DIR, "logs", "{sample}_ltrim_log.txt"),
        rfilter = join(BBDUK_DIR, "logs", "{sample}_rfilter_log.txt"),
        rtrim = join(BBDUK_DIR, "logs", "{sample}_rtrim_log.txt")
    shell:
        """
        # Left filter - find reads with vector sequence
        bbduk.sh in={input.fastq} outm={output.lfilter} literal={params.vector} \
            k={params.k} rcomp={params.rcomp} mm={params.mm} hdist={params.hdist} 2> {log.lfilter}
        
        # Left trim - trim off vector sequence
        bbduk.sh in={output.lfilter} out={output.ltrim} literal={params.vector} \
            ktrim=l k={params.k} rcomp={params.rcomp} hdist={params.hdist} 2> {log.ltrim}
        
        # Right filter - find reads with backbone sequence
        bbduk.sh in={output.ltrim} outm={output.rfilter} literal={params.backbone} \
            k={params.k} rcomp={params.rcomp} mm={params.mm} hdist={params.hdist} 2> {log.rfilter}
        
        # Right trim - trim backbone and keep only reads with specific length
        bbduk.sh in={output.rfilter} out={output.rtrim} outm={output.rtrim_miss} \
            literal={params.backbone} ktrim=r k={params.k} rcomp={params.rcomp} hdist={params.hdist} \
            ml={params.min_length} maxLength={params.max_length} 2> {log.rtrim}
        """

# Generate BBDuk processing summary
rule bbduk_summary:
    input:
        raw_fastqs = expand("{sample}.fastq", sample=SAMPLES),
        lfilter_fastqs = expand(join(BBDUK_DIR, "{sample}_lfilter.fastq"), sample=SAMPLES),
        ltrim_fastqs = expand(join(BBDUK_DIR, "{sample}_ltrim.fastq"), sample=SAMPLES),
        rfilter_fastqs = expand(join(BBDUK_DIR, "{sample}_rfilter.fastq"), sample=SAMPLES),
        rtrim_fastqs = expand(join(BBDUK_DIR, "{sample}_rtrim.fastq"), sample=SAMPLES),
        rtrim_miss_fastqs = expand(join(BBDUK_DIR, "{sample}_rtrim_miss.fastq"), sample=SAMPLES),
        logs = expand(join(BBDUK_DIR, "logs", "{sample}_{step}_log.txt"), 
                    sample=SAMPLES, step=["lfilter", "ltrim", "rfilter", "rtrim"])
    output:
        html = join(BBDUK_DIR, "bbduk_summary_report.html"),
        stats = join(BBDUK_DIR, "bbduk_summary_stats.tsv")
    params:
        vector = VECTOR,
        backbone = BACKBONE,
        samples = SAMPLES
    script:
        "scripts/generate_bbduk_summary.py"

# MAGeCK count rule
rule mageck_count:
    input:
        fastqs = expand(join(BBDUK_DIR, "{sample}_rtrim.fastq"), sample=SAMPLES),
        library = LIBRARY,
        # Add dependency on the bbduk summary to ensure it runs before mageck
        bbduk_summary = join(BBDUK_DIR, "bbduk_summary_stats.tsv")
    output:
        count_file = join(MAGECK_DIR, "mageck_counts.count.txt")
    params:
        prefix = join(MAGECK_DIR, "mageck_counts"),
        sample_labels = ",".join(SAMPLES)
    threads: config["threads"]["mageck"]
    resources:
        mem_mb = config["memory"]["mageck"]
    shell:
        """
        mageck count -l {input.library} \
            -n {params.prefix} \
            --sample-label {params.sample_labels} \
            --fastq {input.fastqs}
        """

# MAGeCK test rule for each comparison
rule mageck_test:
    input:
        count_file = join(MAGECK_DIR, "mageck_counts.count.txt")
    output:
        gene_summary = join(MAGECK_DIR, "{comparison}.gene_summary.txt"),
        sgrna_summary = join(MAGECK_DIR, "{comparison}.sgrna_summary.txt")
    params:
        prefix = join(MAGECK_DIR, "{comparison}"),
        # Support both string and list formats for test and control samples
        test_sample = lambda wildcards: (
            ",".join(config["comparisons"][wildcards.comparison]["test"]) 
            if isinstance(config["comparisons"][wildcards.comparison]["test"], list)
            else config["comparisons"][wildcards.comparison]["test"]
        ),
        control_sample = lambda wildcards: (
            ",".join(config["comparisons"][wildcards.comparison]["control"])
            if isinstance(config["comparisons"][wildcards.comparison]["control"], list)
            else config["comparisons"][wildcards.comparison]["control"]
        ),
        # Get additional MAGeCK parameters if provided
        fdr = config.get("mageck_fdr", 0.05),
        normalization_method = config.get("mageck_normalization_method", "median")
    threads: config["threads"]["mageck"]
    resources:
        mem_mb = config["memory"]["mageck"]
    shell:
        """
        mageck test -k {input.count_file} \
            -t {params.test_sample} \
            -c {params.control_sample} \
            --normchoice {params.normalization_method} \
            --fdr {params.fdr} \
            -n {params.prefix}
        """