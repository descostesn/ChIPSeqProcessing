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
  conda: "envs/fastqc.yaml"
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
  conda: "envs/fastqc.yaml"
  benchmark: "benchmark/fastqc_paired_file/{pairedEndName}.tsv"
  shell:
    """
    fastqc -o ../results/qc/fastqc/raw_fastq/paired/ {input.fq1}
    fastqc -o ../results/qc/fastqc/raw_fastq/paired/ {input.fq2}
    """