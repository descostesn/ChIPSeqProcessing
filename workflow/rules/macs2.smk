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

    drwxr-sr-x 5 descoste boulard 4096 Sep 15  2023 0.001
drwxr-sr-x 5 descoste boulard 4096 Sep 15  2023 0.01
drwxr-sr-x 5 descoste boulard 4096 Sep 15  2023 0.02
drwxr-sr-x 5 descoste boulard 4096 Sep 15  2023 0.03
drwxr-sr-x 5 descoste boulard 4096 Sep 15  2023 0.04
drwxr-sr-x 5 descoste boulard 4096 Sep 15  2023 1e-04
drwxr-sr-x 5 descoste boulard 4096 Sep 15  2023 1e-05


rule macs2_broad_single:
  input:
    rules.remove_duplicates_bowtie_single.output.bamFile
  output:
    "../results/peak_detection/single/macs2/{genome}/{qvalthres}/no_model_broad/{singlebestmulti}_peaks.broadPeak"

