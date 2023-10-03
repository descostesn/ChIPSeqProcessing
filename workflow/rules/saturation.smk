################################################################################
# Performing plots of the number of multireads according to the mode of
# alignment
# WARNING: This code is made to work with TEbench
################################################################################

rule saturation_plot_bowtie_single:
  input:
    logVecBest = "../results/bam/single/bowtie2_results/{genome}/{singleEndName}_trimmed_best.log",
    logVecMulti = expand("../results/bam/single/bowtie2_results/{{genome}}/{{singleEndName}}_trimmed_k{multi}.log", multi=MULTITHRESHOLD)
  output:
    pngSingle = "../results/qc/bowtie2_saturation_multireads/{genome}/single/{singleEndName}.png",
    reportSingle = "../results/qc/bowtie2_saturation_multireads/{genome}/single/{singleEndName}.txt"
  threads: 1
  conda: "../envs/rCoreAndLibraries.yaml"
  benchmark: "benchmark/saturation_plot_bowtie_single/{singleEndName}.tsv"
  params:
    type = "single"
  script:
    "../scripts/saturation_plot.R" 


rule saturation_plot_bowtie_paired:
  input:
    logVecBest = "../results/bam/paired/bowtie2_results/{genome}/{pairedEndName}_trimmed_best.log",
    logVecMulti = expand("../results/bam/paired/bowtie2_results/{{genome}}/{{pairedEndName}}_trimmed_k{multi}.log", multi=MULTITHRESHOLD)
  output:
    pngPaired = "../results/qc/bowtie2_saturation_multireads/{genome}/paired/{pairedEndName}.png",
    reportPaired = "../results/qc/bowtie2_saturation_multireads/{genome}/paired/{pairedEndName}.txt"
  threads: 1
  conda: "../envs/rCoreAndLibraries.yaml"
  benchmark: "benchmark/saturation_plot_bowtie_paired/{pairedEndName}.tsv"
  params:
    type = "paired"
  script:
    "../scripts/saturation_plot.R" 
