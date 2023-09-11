rule sort_bowtie_single:
  input:
    bamFile = "../results/bam/single/bowtie2_results/{genome}/{singlebestmultiall}.bam"
  output:
    bam = "../results/bam/single/bowtie2_results/{genome}/{singlebestmultiall}_sorted.bam"
  threads: 24
  conda: "../envs/samtools.yaml"
  benchmark: "benchmark/sort_bowtie_single/{genome}/{singlebestmultiall}.tsv"
  params:
    memory="2G"
  shell:
    """
    samtools sort -m {params.memory} --threads {threads} -O bam -o {output} {input}
    """
    
rule sort_bowtie_paired:
  input:
    bamFile = "../results/bam/paired/bowtie2_results/{genome}/{pairedbestmultiall}.bam"
  output:
    bam = "../results/bam/paired/bowtie2_results/{genome}/{pairedbestmultiall}_sorted.bam"
  threads: 24
  conda: "../envs/samtools.yaml"
  benchmark: "benchmark/sort_bowtie_paired/{genome}/{pairedbestmultiall}.tsv"
  params:
    memory="2G"
  shell:
    """
    samtools sort -m {params.memory} --threads {threads} -O bam -o {output} {input}
    """
  