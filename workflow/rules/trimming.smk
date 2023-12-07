################################################################################
# Trimming the quality of the fastq file for a minimum value of 20.
# WARNING: This code is made to work with another workflow
################################################################################

rule trimming_single:
  input:
    rules.gzip_fastq_single.output.singleFastq
  output:
    trimmedFq = "../results/data/single/trimmed/{singleEndName}_trimmed.fq.gz",
    report = "../results/data/single/trimmed/{singleEndName}.fastq.gz_trimming_report.txt"
  threads: 6
  conda: "../envs/trimgalore.yaml"
  benchmark: "benchmark/trimming_single/{singleEndName}.tsv"
  shell:
    "trim_galore --quality 20 -o ../results/data/single/trimmed/ --cores 4 {input}"
  
  
rule trimming_paired:
  input:
    fq1 = "../results/data/paired/{pairedEndName}_1.fastq.gz",
    fq2 = "../results/data/paired/{pairedEndName}_2.fastq.gz"
  output:
    trimmedFq1 = "../results/data/paired/trimmed/{pairedEndName}_1_val_1.fq.gz",
    report1 = "../results/data/paired/trimmed/{pairedEndName}_1.fastq.gz_trimming_report.txt",
    trimmedFq2 = "../results/data/paired/trimmed/{pairedEndName}_2_val_2.fq.gz",
    report2 = "../results/data/paired/trimmed/{pairedEndName}_2.fastq.gz_trimming_report.txt"
  threads: 6
  conda: "../envs/trimgalore.yaml"
  benchmark: "benchmark/trimming_paired/{pairedEndName}.tsv"
  shell:
    "trim_galore --quality 20 -o ../results/data/paired/trimmed/ --cores 4 --paired {input.fq1} {input.fq2}"
