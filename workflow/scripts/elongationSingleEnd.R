##########
# This script uses cross-correlation plots to estimate the fragment length of 
# single-end experiments.
# Descostes December 2019
##########

library(csaw)

#############
## PARAMS
#############


bamfile <- snakemake@input$bamFile
save(bamfile, file="/g/romebioinfo/Projects/TEbench/workflow/jobs/ChIPSeqProcessing_retrieve_elongationSize_single/bamfile.Rdat")
stop("saving for testing")
load("/g/romebioinfo/Projects/TEbench/workflow/jobs/ChIPSeqProcessing_retrieve_elongationSize_single/bamfile.Rdat")


#############
## MAIN
#############


maxdelay <- 500
param <- readParam(minq = 20)

dedupon <- reform(param, dedup=TRUE)
x <- correlateReads(bamfile, maxdelay, param=dedupon)
plot(0:maxdelay, x, type = "l", ylab = "CCF", xlab = "Delay (bp)")

extensionsize <- maximizeCcf(x)
outputfile <- paste0("results/elongation_size_single_csaw/",
    strsplit(basename(bamfile), "\\.bam")[[1]], ".txt")

cat("Writing output to: ", outputfile, "\n")

write(extensionsize, file = outputfile, ncolumns = 1)

Sys.sleep(10)
