################
# This script counts the number of reads having a number of matches equal
# to a threshold k. It also creates an histogram showing the number of sequences
# having x matches for each experiment.
#
# Descostes - Nov 2023 - R 4.2.2-foss-2022b
################

library(Rsamtools)
library(parallel)



##################
# PARAMETERS
##################

bamvec <- c("/g/romebioinfo/Projects/TEbench/results/bam/paired/bowtie2_results/mm39/H3K9me3_SRR2136759_GSM1841035_GSE71589_ESC_bio2_trimmed_best_sorted_noDups.bam",
"/g/romebioinfo/Projects/TEbench/results/bam/paired/bowtie2_results/mm39/H3K9me3_SRR2136759_GSM1841035_GSE71589_ESC_bio2_trimmed_k10_sorted_noDups.bam",
"/g/romebioinfo/Projects/TEbench/results/bam/paired/bowtie2_results/mm39/H3K9me3_SRR2136759_GSM1841035_GSE71589_ESC_bio2_trimmed_k50_sorted_noDups.bam",
"/g/romebioinfo/Projects/TEbench/results/bam/paired/bowtie2_results/mm39/H3K9me3_SRR2136759_GSM1841035_GSE71589_ESC_bio2_trimmed_k100_sorted_noDups.bam",
"/g/romebioinfo/Projects/TEbench/results/bam/paired/bowtie2_results/mm39/H3K9me3_SRR2136759_GSM1841035_GSE71589_ESC_bio2_trimmed_k150_sorted_noDups.bam",
"/g/romebioinfo/Projects/TEbench/results/bam/paired/bowtie2_results/mm39/H3K9me3_SRR2136759_GSM1841035_GSE71589_ESC_bio2_trimmed_k200_sorted_noDups.bam",
"/g/romebioinfo/Projects/TEbench/results/bam/paired/bowtie2_results/mm39/H3K9me3_SRR2136759_GSM1841035_GSE71589_ESC_bio2_trimmed_k250_sorted_noDups.bam",
"/g/romebioinfo/Projects/TEbench/results/bam/paired/bowtie2_results/mm39/H3K9me3_SRR2136759_GSM1841035_GSE71589_ESC_bio2_trimmed_k300_sorted_noDups.bam",
"/g/romebioinfo/Projects/TEbench/results/bam/paired/bowtie2_results/mm39/H3K9me3_SRR2136759_GSM1841035_GSE71589_ESC_bio2_trimmed_k350_sorted_noDups.bam",
"/g/romebioinfo/Projects/TEbench/results/bam/paired/bowtie2_results/mm39/H3K9me3_SRR2136759_GSM1841035_GSE71589_ESC_bio2_trimmed_k400_sorted_noDups.bam")

expname <- "H3K9me3_SRR2136759_GSM1841035_GSE71589_ESC_bio2"

namesbamvec <- c("kbest", "k10","k50", "k100", "k150", "k200",
    "k250", "k300", "k350", "k400")

chromvec <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12",
"13", "14", "15", "16", "17", "18", "19", "X", "Y")

ncores <- 21


##################
# MAIN
##################

# bamFile <- bamvec[10]
# bamName <- namesbamvec[10]
mapply(function(bamFile, bamName, expname, chromvec, ncores) {

    message("Reading ", expname, "-", bamName)

    ## Importing bam reads per chromosome
    #chrom <- chromvec[1]
    freqseqlist <- parallel::mclapply(chromvec, function(chrom) {

        message("\t chr", chrom)
        ## Retrieving the sequence name (rname) on the defined chrom (qname)
        params <- Rsamtools::ScanBamParam(what = c('rname', 'qname'),
            which = GRanges(chrom, IRanges(1, 1e8)),
            flag = scanBamFlag(isUnmappedQuery = FALSE))
        message("\t\t Reading")
        x <- Rsamtools::scanBam(bamFile, index = bamFile, param = params)[[1]]
        chromread <- as.character(unique(x[[2]]))

        if (!isTRUE(all.equal(length(chromread), 1)) ||
            !isTRUE(all.equal(chromread, chrom)))
            stop("Problem in reading the bam file per chromosome")

        message("\t\t Counting")
        freqseq <- table(x[[1]])

        return(freqseq)
    }, mc.cores = ncores)

}, bamvec, namesbamvec, MoreArgs = list(expname, chromvec, ncores))