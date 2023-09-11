################################################################################
# Converting sam files obtained with Bowtie2 to Bam and sorting.
# WARNING: This code is made to work with TEbench
################################################################################

rule sam2bam_best_single:
  input:
    rules.bowtie2_best_single.output.sam
  output:
    temp("../results/bam/single/bowtie2_results/{genome}/{singleEndName}_trimmed_best.bam"),
  threads: 1
  conda: "../envs/samtools.yaml"
  benchmark: "benchmark/sam2bam/{genome}/{singleEndName}_trimmed_best.tsv"
  shell:
    """
    samtools view -h -b {input} > {output}
    """

rule sam2bam_MR_single:
  input:
    rules.bowtie2_MR_single.output.sam
  output:
    temp("../results/bam/single/bowtie2_results/{genome}/{singleEndName}_trimmed_k{multi, [0-9]+}.bam")
  threads: 1
  conda: "../envs/samtools.yaml"
  benchmark: "benchmark/sam2bam/{genome}/{singleEndName}_trimmed_k{multi}.tsv"
  shell:
    """
    samtools view -h -b {input} > {output}
    """

rule sam2bam_all_single:
  input:
    rules.bowtie2_all_single.output.sam
  output:
    temp("../results/bam/single/bowtie2_results/{genome}/{singleEndName}_trimmed_all.bam")
  threads: 1
  conda: "../envs/samtools.yaml"
  benchmark: "benchmark/sam2bam/{genome}/{singleEndName}_trimmed_all.tsv"
  shell:
    """
    samtools view -h -b {input} > {output}
    """

rule sam2bam_best_paired:
  input:
    rules.bowtie2_best_paired.output.sam
  output:
    temp("../results/bam/paired/bowtie2_results/{genome}/{pairedEndName}_trimmed_best.bam")
  threads: 1
  conda: "../envs/samtools.yaml"
  benchmark: "benchmark/sam2bam/{genome}/{pairedEndName}_trimmed_best.tsv"
  shell:
    """
    samtools view -h -b {input} > {output}
    """

rule sam2bam_MR_paired:
  input:
    rules.bowtie2_MR_paired.output.sam
  output:
    temp("../results/bam/paired/bowtie2_results/{genome}/{pairedEndName}_trimmed_k{multi, [0-9]+}.bam")
  threads: 1
  conda: "../envs/samtools.yaml"
  benchmark: "benchmark/sam2bam/{genome}/{pairedEndName}_trimmed_k{multi}.tsv"
  shell:
    """
    samtools view -h -b {input} > {output}
    """

rule sam2bam_all_paired:
  input:
    rules.bowtie2_all_paired.output.sam
  output:
    temp("../results/bam/paired/bowtie2_results/{genome}/{pairedEndName}_trimmed_all.bam")
  threads: 1
  conda: "../envs/samtools.yaml"
  benchmark: "benchmark/sam2bam/{genome}/{pairedEndName}_trimmed_all.tsv"
  shell:
    """
    samtools view -h -b {input} > {output}
    """
