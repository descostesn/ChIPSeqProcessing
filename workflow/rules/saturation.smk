################################################################################
# Performing plots of the number of multireads according to the mode of
# alignment
# WARNING: This code is made to work with another workflow
################################################################################

rule saturation_plot_bowtie_single:
  input:
    logVecBest = "../results/bam/single/bowtie2_results/{genome}/{singleEndName}_trimmed_best.log",
    logVecMulti = expand("../results/bam/single/bowtie2_results/{{genome}}/{{singleEndName}}_trimmed_k{multi}.log", multi=MULTITHRESHOLD)
  output:
    png = "../results/qc/bowtie2_saturation_percentmultireads/{genome}/single/{singleEndName}.png",
    report = "../results/qc/bowtie2_saturation_percentmultireads/{genome}/single/{singleEndName}.txt"
  threads: 1
  conda: "../envs/rCoreAndLibraries.yaml"
  benchmark: "benchmark/saturation_plot_bowtie_single/{genome}/{singleEndName}.tsv"
  params:
    type = "single"
  script:
    "../scripts/saturation_plot.R" 


rule saturation_plot_bowtie_paired:
  input:
    logVecBest = "../results/bam/paired/bowtie2_results/{genome}/{pairedEndName}_trimmed_best.log",
    logVecMulti = expand("../results/bam/paired/bowtie2_results/{{genome}}/{{pairedEndName}}_trimmed_k{multi}.log", multi=MULTITHRESHOLD)
  output:
    png = "../results/qc/bowtie2_saturation_percentmultireads/{genome}/paired/{pairedEndName}.png",
    report = "../results/qc/bowtie2_saturation_percentmultireads/{genome}/paired/{pairedEndName}.txt"
  threads: 1
  conda: "../envs/rCoreAndLibraries.yaml"
  benchmark: "benchmark/saturation_plot_bowtie_paired/{genome}/{pairedEndName}.tsv"
  params:
    type = "paired"
  script:
    "../scripts/saturation_plot.R" 


#!!!!!!!!!!!!!!!!!!!!
#Do rule to retrieve a table of nb of alignment per read
#!!!!!!!!!!!!!!!!!!!!!!