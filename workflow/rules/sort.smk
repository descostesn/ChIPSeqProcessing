################################################################################
# Sorting the bam files obtained with bowtie2 before indexing them
# WARNING: This code is made to work with TEbench
################################################################################

rule sort_bowtie_single:
  input:
    "../results/bam/single/bowtie2_results/{genome}/{singlebestmulti}.bam"
  output:
    temp("../results/bam/single/bowtie2_results/{genome}/{singlebestmulti}_sorted.bam")
  threads: 24
  conda: "../envs/samtools.yaml"
  benchmark: "benchmark/sort_bowtie_single/{genome}/{singlebestmulti}.tsv"
  params:
    memory="2G"
  shell:
    """
    samtools sort -m {params.memory} --threads {threads} -O bam -o {output} {input}
    """
    
rule sort_bowtie_paired:
  input:
    "../results/bam/paired/bowtie2_results/{genome}/{pairedbestmulti}.bam"
  output:
    temp("../results/bam/paired/bowtie2_results/{genome}/{pairedbestmulti}_sorted.bam")
  threads: 24
  conda: "../envs/samtools.yaml"
  benchmark: "benchmark/sort_bowtie_paired/{genome}/{pairedbestmulti}.tsv"
  params:
    memory="2G"
  shell:
    """
    samtools sort -m {params.memory} --threads {threads} -O bam -o {output} {input}
    """
  