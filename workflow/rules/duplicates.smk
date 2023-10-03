################################################################################
# Removing duplicates tags from the bam files.
# WARNING: This code is made to work with TEbench
################################################################################

rule remove_duplicates_bowtie_single:
  input:
    rules.sort_bowtie_single.output
  output:
  # put back protected for bamfile
    bamFile = "../results/bam/single/bowtie2_results/{genome}/{singlebestmulti}_sorted_nodups.bam",
    report = "../results/bam/single/bowtie2_results/{genome}/{singlebestmulti}_sorted_nodups.txt"
  threads: 1
  conda: "../envs/picard.yaml"
  benchmark: "benchmark/remove_duplicates_bowtie_single/{genome}/{singlebestmulti}.tsv"
  shell:
    "picard MarkDuplicates I={input} O={output.bamFile} M={output.report} TAGGING_POLICY=All REMOVE_DUPLICATES=true ASSUME_SORT_ORDER=coordinate"
  
  
rule remove_duplicates_bowtie_paired:
  input:
    rules.sort_bowtie_paired.output
  output:
  # put back protected for bamfile
    bamFile = "../results/bam/paired/bowtie2_results/{genome}/{pairedbestmulti}_sorted_noDups.bam",
    report = "../results/bam/paired/bowtie2_results/{genome}/{pairedbestmulti}_sorted_noDups.txt"
  threads: 1
  conda: "../envs/picard.yaml"
  benchmark: "benchmark/remove_duplicates_bowtie_paired/{genome}/{pairedbestmulti}.tsv"
  shell:
    "picard MarkDuplicates I={input} O={output.bamFile} M={output.report} TAGGING_POLICY=All REMOVE_DUPLICATES=true ASSUME_SORT_ORDER=coordinate"
