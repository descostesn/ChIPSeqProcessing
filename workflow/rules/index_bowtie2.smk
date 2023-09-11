rule index_bam_bowtie_single:
  input:
    rules.sort_bowtie_single.output
  output:
    "../results/bam/single/bowtie2_results/{genome}/{singlebestmultiall}_sorted.bam.bai"
  threads: 1
  conda: "../envs/samtools.yaml"
  benchmark: "benchmark/index_bam_bowtie_single/{genome}/{singlebestmultiall}.tsv"
  shell:
    "samtools index {input}" 

rule index_bam_bowtie_paired:
  input:
    rules.sort_bowtie_paired.output
  output:
    "../results/bam/paired/bowtie2_results/{genome}/{pairedbestmultiall}_sorted.bam.bai"
  threads: 1
  conda: "../envs/samtools.yaml"
  benchmark: "benchmark/index_bam_bowtie_paired/{genome}/{pairedbestmultiall}.tsv"
  shell:
    "samtools index {input}"
