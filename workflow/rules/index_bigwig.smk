################################################################################
# Performing indexing of the sorted and duplicate removed bam files to 
# generate the bigwig files.
# WARNING: This code is made to work with another workflow
################################################################################

rule index_bam_bigwig_single:
  input:
    rules.remove_duplicates_bowtie_single.output.bamFile
  output:
    "../results/bam/single/bowtie2_results/{genome}/{singlebestmulti}_sorted_nodups.bam.bai"
  threads: 1
  conda: "../envs/samtools.yaml"
  benchmark: "benchmark/index_bam_bigwig_single/{genome}/{singlebestmulti}.tsv"
  shell:
    "samtools index {input}"


rule index_bam_bigwig_paired:
  input:
    rules.remove_duplicates_bowtie_paired.output.bamFile
  output:
    "../results/bam/paired/bowtie2_results/{genome}/{pairedbestmulti}_sorted_noDups.bam.bai"
  threads: 1
  conda: "../envs/samtools.yaml"
  benchmark: "benchmark/index_bam_bigwig_paired/{genome}/{pairedbestmulti}.tsv"
  shell:
    "samtools index {input}"    