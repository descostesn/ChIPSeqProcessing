################################################################################
# Building bowtie2 index.
# WARNING: This code is made to work with TEbench
################################################################################

rule build_bowtie_index:
  input:
    rules.download_genome_fasta.output
  output:
    "../results/data/bowtie2_index/{genome}/{prefix}.1.bt2",
    "../results/data/bowtie2_index/{genome}/{prefix}.2.bt2",
    "../results/data/bowtie2_index/{genome}/{prefix}.3.bt2",
    "../results/data/bowtie2_index/{genome}/{prefix}.4.bt2",
    "../results/data/bowtie2_index/{genome}/{prefix}.rev.1.bt2",
    "../results/data/bowtie2_index/{genome}/{prefix}.rev.2.bt2"
  threads: 8
  conda: "../envs/bowtie2.yaml"
  benchmark: "benchmark/build_bowtie_index/{genome}.{prefix}.tsv"
  params:
    outputFolder="../results/data/bowtie2_index"
  shell:
    """
    mkdir -p {params.outputFolder}/{wildcards.genome}
    bowtie2-build -f --threads {threads} {input} {params.outputFolder}/{wildcards.genome}/{wildcards.prefix}
    """


################################################################################
# Performing alignement of the reads with bowtie2 in three modes:
# 1) Best (best alignment is kept)
# 2) -k (at most k alignments are kept, several values for k)
# 3) -a (all aligments are kept)
# WARNING: This code is made to work with TEbench
################################################################################


rule bowtie2_best_single:
  input:
    rules.trimming_single.output.trimmedFq
  output:
    sam = temp("../results/bam/single/bowtie2_results/{genome}/{singleEndName}_trimmed_best.sam"),
    log = "../results/bam/single/bowtie2_results/{genome}/{singleEndName}_trimmed_best.log"
  threads: 20
  conda: "../envs/bowtie2.yaml"
  benchmark: "benchmark/bowtie2/single/best/{genome}/{singleEndName}_best.tsv"
  params:
    indexPath = lambda wildcards: expand("../results/data/bowtie2_index/{genome}/{prefix}", genome = GENOMEID, prefix = PREFIXFASTAGTF)
  shell:
    """
    bowtie2 -q -p {threads} -x {params.indexPath} -U {input} -S {output.sam} 2> {output.log}
    """

rule bowtie2_MR_single:
  input:
    rules.trimming_single.output.trimmedFq
  output:
    sam = temp("../results/bam/single/bowtie2_results/{genome}/{singleEndName}_trimmed_k{multi}.sam"),
    log = "../results/bam/single/bowtie2_results/{genome}/{singleEndName}_trimmed_k{multi}.log"
  threads: 20
  conda: "../envs/bowtie2.yaml"
  benchmark: "benchmark/bowtie2/single/thresholdK/{genome}/{singleEndName}_trimmed_k{multi}.tsv"
  params:
    indexPath = lambda wildcards: expand("../results/data/bowtie2_index/{genome}/{prefix}", genome = GENOMEID, prefix = PREFIXFASTAGTF),
    multiThreshold = "{multi}"
  shell:
    """
    bowtie2 -q -p {threads} -x {params.indexPath} -k {params.multiThreshold} -U {input} -S {output.sam} 2> {output.log}
    """

rule bowtie2_best_paired:
  input:
    trimmedFq1 = "../results/data/paired/trimmed/{pairedEndName}_1_val_1.fq.gz",
    trimmedFq2 = "../results/data/paired/trimmed/{pairedEndName}_2_val_2.fq.gz"
  output:
    sam = temp("../results/bam/paired/bowtie2_results/{genome}/{pairedEndName}_trimmed_best.sam"),
    log = "../results/bam/paired/bowtie2_results/{genome}/{pairedEndName}_trimmed_best.log"
  threads: 20
  conda: "../envs/bowtie2.yaml"
  benchmark: "benchmark/bowtie2/paired/best/{genome}/{pairedEndName}_trimmed_best.tsv"
  params:
    indexPath = lambda wildcards: expand("../results/data/bowtie2_index/{genome}/{prefix}", genome = GENOMEID, prefix = PREFIXFASTAGTF)
  shell:
    """
    bowtie2 -q -p {threads} -x {params.indexPath} --no-mixed --no-discordant --dovetail -1 {input.trimmedFq1} -2 {input.trimmedFq2} -S {output.sam} 2> {output.log}
    """

rule bowtie2_MR_paired:
  input:
    trimmedFq1 = "../results/data/paired/trimmed/{pairedEndName}_1_val_1.fq.gz",
    trimmedFq2 = "../results/data/paired/trimmed/{pairedEndName}_2_val_2.fq.gz"
  output:
    sam = temp("../results/bam/paired/bowtie2_results/{genome}/{pairedEndName}_trimmed_k{multi}.sam"),
    log = "../results/bam/paired/bowtie2_results/{genome}/{pairedEndName}_trimmed_k{multi}.log"
  threads: 24
  conda: "../envs/bowtie2.yaml"
  benchmark: "benchmark/bowtie2/paired/thresholdK/{genome}/{pairedEndName}_trimmed_k{multi}.tsv"
  params:
    indexPath = lambda wildcards: expand("../results/data/bowtie2_index/{genome}/{prefix}", genome = GENOMEID, prefix = PREFIXFASTAGTF),
    multiThreshold = "{multi}"
  shell:
    """
    bowtie2 -q -p {threads} -x {params.indexPath} -k {params.multiThreshold} --no-mixed --no-discordant --dovetail -1 {input.trimmedFq1} -2 {input.trimmedFq2} -S {output.sam} 2> {output.log}
    """
