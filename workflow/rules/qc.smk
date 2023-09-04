################################################################################
# Performing fastqc before trimming according to quality.
# WARNING: This code is made to work with TEbench
################################################################################

rule fastqc_single_file:
  input:
    rules.gzip_fastq_single.output.singleFastq
  output:
    "../results/qc/fastqc/raw_fastq/single/{singleEndName}_fastqc.html",
    "../results/qc/fastqc/raw_fastq/single/{singleEndName}_fastqc.zip"
  threads: 1
  conda: "../envs/fastqc.yaml"
  benchmark: "benchmark/fastqc_single_file/{singleEndName}.tsv"
  shell:
    "fastqc -o ../results/qc/fastqc/raw_fastq/single/ {input}"

rule fastqc_paired_file:
  input:
    fq1 = "../results/data/paired/{pairedEndName}_1.fastq.gz",
    fq2 = "../results/data/paired/{pairedEndName}_2.fastq.gz" 
  output:
    "../results/qc/fastqc/raw_fastq/paired/{pairedEndName}_1_fastqc.html",
    "../results/qc/fastqc/raw_fastq/paired/{pairedEndName}_1_fastqc.zip",
    "../results/qc/fastqc/raw_fastq/paired/{pairedEndName}_2_fastqc.html",
    "../results/qc/fastqc/raw_fastq/paired/{pairedEndName}_2_fastqc.zip"
  threads: 1
  conda: "../envs/fastqc.yaml"
  benchmark: "benchmark/fastqc_paired_file/{pairedEndName}.tsv"
  shell:
    """
    fastqc -o ../results/qc/fastqc/raw_fastq/paired/ {input.fq1}
    fastqc -o ../results/qc/fastqc/raw_fastq/paired/ {input.fq2}
    """

################################################################################
# Performing fastqc after trimming according to quality.
# WARNING: This code is made to work with TEbench
################################################################################

rule fastqc_after_trim_single:
  input:
    rules.trimming_single.output.trimmedFq
  output:  
    "../results/qc/fastqc/trimmed_fastq/single/{singleEndName}_trimmed_fastqc.html",
    "../results/qc/fastqc/trimmed_fastq/single/{singleEndName}_trimmed_fastqc.zip"
  threads: 1
  conda: "../envs/fastqc.yaml"
  benchmark: "benchmark/fastqc_after_trim_single/{singleEndName}.tsv"
  shell:
    "fastqc -o ../results/qc/fastqc/trimmed_fastq/single/ {input}"
    
    
rule fastqc_after_trim_paired:
  input:
    fq1 = "../results/data/paired/trimmed/{pairedEndName}_1_val_1.fq.gz",
    fq2 = "../results/data/paired/trimmed/{pairedEndName}_2_val_2.fq.gz"
  output:
    "../results/qc/fastqc/trimmed_fastq/paired/{pairedEndName}_1_val_1_fastqc.html",
    "../results/qc/fastqc/trimmed_fastq/paired/{pairedEndName}_1_val_1_fastqc.zip",
    "../results/qc/fastqc/trimmed_fastq/paired/{pairedEndName}_2_val_2_fastqc.html",
    "../results/qc/fastqc/trimmed_fastq/paired/{pairedEndName}_2_val_2_fastqc.zip"
  threads: 1
  conda: "../envs/fastqc.yaml"
  benchmark: "benchmark/fastqc_after_trim_paired/{pairedEndName}.tsv"
  shell:
    """
    fastqc -o ../results/qc/fastqc/trimmed_fastq/paired/ {input.fq1}
    fastqc -o ../results/qc/fastqc/trimmed_fastq/paired/ {input.fq2}
    """
