################
# This script counts the number of reads having a number of matches equal
# to a threshold k. It also creates an histogram showing the number of sequences
# having x matches for each experiment.
#
# Descostes - Nov 2023 -R 4.3.2 (mamba env)
################

library(Rsamtools)
library(parallel)
library(collapse)
library(GenomicFiles)
library(GenomicRanges)
library(IRanges)

##################
# PARAMETERS
##################

bamvec <- c(
    "/g/romebioinfo/Projects/TEbench/results/bam/paired/bowtie2_results/mm39/H3K9me3_SRR2136759_GSM1841035_GSE71589_ESC_bio2_trimmed_best_sorted_noDups.bam",
    "/g/romebioinfo/Projects/TEbench/results/bam/paired/bowtie2_results/mm39/H3K9me3_SRR2136759_GSM1841035_GSE71589_ESC_bio2_trimmed_k10_sorted_noDups.bam",
    "/g/romebioinfo/Projects/TEbench/results/bam/paired/bowtie2_results/mm39/H3K9me3_SRR2136759_GSM1841035_GSE71589_ESC_bio2_trimmed_k50_sorted_noDups.bam",
    "/g/romebioinfo/Projects/TEbench/results/bam/paired/bowtie2_results/mm39/H3K9me3_SRR2136759_GSM1841035_GSE71589_ESC_bio2_trimmed_k100_sorted_noDups.bam",
    "/g/romebioinfo/Projects/TEbench/results/bam/paired/bowtie2_results/mm39/H3K9me3_SRR2136759_GSM1841035_GSE71589_ESC_bio2_trimmed_k150_sorted_noDups.bam",
    "/g/romebioinfo/Projects/TEbench/results/bam/paired/bowtie2_results/mm39/H3K9me3_SRR2136759_GSM1841035_GSE71589_ESC_bio2_trimmed_k200_sorted_noDups.bam",
    "/g/romebioinfo/Projects/TEbench/results/bam/paired/bowtie2_results/mm39/H3K9me3_SRR2136759_GSM1841035_GSE71589_ESC_bio2_trimmed_k250_sorted_noDups.bam",
    "/g/romebioinfo/Projects/TEbench/results/bam/paired/bowtie2_results/mm39/H3K9me3_SRR2136759_GSM1841035_GSE71589_ESC_bio2_trimmed_k300_sorted_noDups.bam",
    "/g/romebioinfo/Projects/TEbench/results/bam/paired/bowtie2_results/mm39/H3K9me3_SRR2136759_GSM1841035_GSE71589_ESC_bio2_trimmed_k350_sorted_noDups.bam",
    "/g/romebioinfo/Projects/TEbench/results/bam/paired/bowtie2_results/mm39/H3K9me3_SRR2136759_GSM1841035_GSE71589_ESC_bio2_trimmed_k400_sorted_noDups.bam"
)

expname <- "H3K9me3_SRR2136759_GSM1841035_GSE71589_ESC_bio2"

namesbamvec <- c(
    "kbest", "k10", "k50", "k100", "k150", "k200",
    "k250", "k300", "k350", "k400"
)

ncores <- 8

outputfold <- "/g/romebioinfo/Projects/TEbench/results/tmp"


##################
# FUNCTIONS
##################

readingonchrom <- function(bamfile) {

    message("\t\t Reading chunk of sequences")

    ## Retrieving the sequence name (rname)
    params <- Rsamtools::ScanBamParam(
        what = c("qname"),
        flag = Rsamtools::scanBamFlag(isUnmappedQuery = FALSE)
    )
    Rsamtools::scanBam(bamfile, param = params)[[1]][["qname"]]
}


fasttable <- function(value) {

    message("\t\t Counting frequencies")

    xf <- collapse::qF(value)
    levelvec <- levels(xf)
    levels(xf) <- seq_len(length(levelvec))
    res <- tabulate(xf)
    names(res) <- levelvec
    return(res)
}

combinefreqtable <- function(x, y) {
    message("Merging chunk to previous one")
     df <- data.frame(id = c(names(x), names(y)), counts = c(x, y))
     res <- tapply(df$counts, df$id, FUN = sum)
     return(res)
}

plotcounts <- function(outputfold1, tablabeloccupancy, nbmatchseqvec) {
    message("\t\t Plotting occupancy and mismatches")
    png(file = file.path(outputfold1, "nbchrommatchesbyseq.png"))
    barplot(tablabeloccupancy,
        xlab = "Number of chromosomes with match per sequence",
        ylab = "Number of sequences"
    )
    dev.off()

    png(file = file.path(outputfold1, "nbseqmatches.png"))
    hist(nbmatchseqvec,
        breaks = 1000, xlab = "Number of matches",
        ylab = "Number of sequences"
    )
    dev.off()
}

##################
# MAIN
##################

register(MulticoreParam(ncores))


freqseqlist <- mapply(function(bamfile, bamname, expname) {

    message("Reading ", expname, "-", bamname)
    bamfile <- Rsamtools::BamFile(bamvec[1], yieldSize = 1000000)
    freqseq <- GenomicFiles::reduceByYield(bamfile, YIELD = readingonchrom,
    MAP = fasttable, REDUCE = combinefreqtable, parallel = TRUE, iterate = TRUE)

    gc(verbose = FALSE)
    return(freqseq)
}, bamvec, namesbamvec, MoreArgs = list(expname))
save(freqseqlist,
    file = "/g/romebioinfo/Projects/TEbench/results/tmp/freqseqlist.Rdat")
