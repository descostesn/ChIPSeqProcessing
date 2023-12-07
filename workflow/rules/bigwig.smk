################################################################################
# Generating the bigwigs for all experiments with and without normalization.
# WARNING: This code is made to work with another workflow
################################################################################


rule bigwig_single:
  input:
    bamFile = "../results/bam/single/bowtie2_results/{genome}/{singlebestmulti}_sorted_nodups.bam",
    baiFile = "../results/bam/single/bowtie2_results/{genome}/{singlebestmulti}_sorted_nodups.bam.bai",
    sizeFile = "../results/qc/elongation_size_single_csaw/{genome}/{singlebestmulti}.txt"
  output:
    "../results/bigwig/single/{genome}/{singlebestmulti}.bw"
  threads: 1
  conda: "../envs/deeptools.yaml"
  benchmark: "benchmark/bigwig_single/{genome}/{singlebestmulti}.tsv"
  shell:
    """
    SIZE=`cat {input.sizeFile}`
    bamCoverage --bam {input.bamFile} --outFileName {output} --outFileFormat bigwig --binSize 50 --numberOfProcessors {threads} --effectiveGenomeSize 2652783500 --extendReads $SIZE
    """

rule bigwig_paired:
  input:
    bamFile = "../results/bam/paired/bowtie2_results/{genome}/{pairedbestmulti}_sorted_noDups.bam",
    baiFile = "../results/bam/paired/bowtie2_results/{genome}/{pairedbestmulti}_sorted_noDups.bam.bai"
  output:
    "../results/bigwig/paired/{genome}/{pairedbestmulti}.bw"
  threads: 1
  conda: "../envs/deeptools.yaml"
  benchmark: "benchmark/bigwig_paired/{genome}/{pairedbestmulti}.tsv"
  shell:
    "bamCoverage --bam {input.bamFile} --outFileName {output} --outFileFormat bigwig --binSize 50 --numberOfProcessors {threads} --effectiveGenomeSize 2652783500 --extendReads"


rule bigwig_norm_single:
  input:
    bamFile = "../results/bam/single/bowtie2_results/{genome}/{singlebestmulti}_sorted_nodups.bam",
    baiFile = "../results/bam/single/bowtie2_results/{genome}/{singlebestmulti}_sorted_nodups.bam.bai",
    sizeFile = "../results/qc/elongation_size_single_csaw/{genome}/{singlebestmulti}.txt"
  output:
    "../results/bigwig/single/{genome}/{singlebestmulti}_norm.bw"
  threads: 1
  conda: "../envs/deeptools.yaml"
  benchmark: "benchmark/bigwig_norm_single/{genome}/{singlebestmulti}.tsv"
  shell:
    """
    SIZE=`cat {input.sizeFile}`
    bamCoverage --bam {input.bamFile} --outFileName {output} --outFileFormat bigwig --binSize 50 --numberOfProcessors {threads} --effectiveGenomeSize 2652783500 --extendReads $SIZE --normalizeUsing RPGC --ignoreForNormalization chrX
    """

rule bigwig_norm_paired:
  input:
    bamFile = "../results/bam/paired/bowtie2_results/{genome}/{pairedbestmulti}_sorted_noDups.bam",
    baiFile = "../results/bam/paired/bowtie2_results/{genome}/{pairedbestmulti}_sorted_noDups.bam.bai"
  output:
    "../results/bigwig/paired/{genome}/{pairedbestmulti}_norm.bw"
  threads: 1
  conda: "../envs/deeptools.yaml"
  benchmark: "benchmark/bigwig_norm_paired/{genome}/{pairedbestmulti}.tsv"
  shell:
    "bamCoverage --bam {input.bamFile} --outFileName {output} --outFileFormat bigwig --binSize 50 --numberOfProcessors {threads} --effectiveGenomeSize 2652783500 --extendReads --normalizeUsing RPGC --ignoreForNormalization chrX"
    