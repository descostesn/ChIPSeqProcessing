################################################################################
# Performing peak detection with macs2
# WARNING: This code is made to work with another workflow
################################################################################

rule macs2_narrow_single:
  input:
    rules.remove_duplicates_bowtie_single.output.bamFile
  output:
    "../results/peak_detection/single/macs2/{genome}/{qvalthres}/{modeltype}/{singlebestmulti}_peaks.narrowPeak"


no_model
model_based




rule macs2_broad_single:
  input:
    rules.remove_duplicates_bowtie_single.output.bamFile
  output:
    "../results/peak_detection/single/macs2/{genome}/{qvalthres}/no_model_broad/{singlebestmulti}_peaks.broadPeak"

