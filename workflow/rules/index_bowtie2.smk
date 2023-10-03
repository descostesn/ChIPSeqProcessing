################################################################################
# Performing indexing of the sorted bam files obtained with bowtie2
# WARNING: This code is made to work with TEbench
################################################################################

rule index_bam_bowtie_single:
  input:
    rules.sort_bowtie_single.output
  output:
    "../results/bam/single/bowtie2_results/{genome}/{singlebestmulti}_sorted.bam.bai"
  threads: 1
  conda: "../envs/samtools.yaml"
  benchmark: "benchmark/index_bam_bowtie_single/{genome}/{singlebestmulti}.tsv"
  shell:
    "samtools index {input}" 

rule index_bam_bowtie_paired:
  input:
    rules.sort_bowtie_paired.output
  output:
    "../results/bam/paired/bowtie2_results/{genome}/{pairedbestmulti}_sorted.bam.bai"
  threads: 1
  conda: "../envs/samtools.yaml"
  benchmark: "benchmark/index_bam_bowtie_paired/{genome}/{pairedbestmulti}.tsv"
  shell:
    "samtools index {input}"
