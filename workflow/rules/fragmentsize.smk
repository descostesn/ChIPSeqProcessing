################################################################################
# Estimate the real fragment size for the single-ended experiment by cross-
# correlation.
# WARNING: This code is made to work with TEbench
################################################################################


rule retrieve_elongationSize_single:
  input:
    bamFile = "../results/bam/single/bowtie2_results/{genome}/{singlebestmulti}_sorted_nodups.bam",
    baiindex = "../results/bam/single/bowtie2_results/{genome}/{singlebestmulti}_sorted_nodups.bam.bai"
  output:
    "../results/qc/elongation_size_single_csaw/{genome}/{singlebestmulti}.txt"
  threads: 1
  conda: "../envs/rCoreAndLibraries.yaml"
  benchmark: "benchmark/retrieve_elongationSize_single/{genome}/{singlebestmulti}.tsv"
  script:
    "../scripts/elongationSingleEnd.R"
