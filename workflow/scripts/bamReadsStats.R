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

chromvec <- c(
    "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12",
    "13", "14", "15", "16", "17", "18", "19", "X", "Y"
)

ncores <- 21

outputfold <- "/g/romebioinfo/Projects/TEbench/results/tmp"


##################
# FUNCTIONS
##################

uniqueormultichrom <- function(currentseq) {
    idx <- which(currentseq != 0)

    if (isTRUE(all.equal(length(idx), 1))) {
        return("unique")
    } else if (isTRUE(all.equal(length(idx), 0))) {
        return("none") # This should not happen
    } else {
        return("several")
    }
}

##################
# MAIN
##################

# bamFile <- bamvec[3]
# bamName <- namesbamvec[3]
mapply(function(bamFile, bamName, expname, chromvec, ncores, outputfold) {
    message("Reading ", expname, "-", bamName)
    outputfold1 <- file.path(outputfold, expname, bamName)

    ## Importing bam reads per chromosome
    # chrom <- chromvec[1]
    freqseqlist <- parallel::mclapply(chromvec, function(chrom, outfold) {
        message("\t chr", chrom)
        ## Retrieving the sequence name (rname) on the defined chrom (qname)
        params <- Rsamtools::ScanBamParam(
            what = c("rname", "qname"),
            which = GRanges(chrom, IRanges(1, 1e8)),
            flag = scanBamFlag(isUnmappedQuery = FALSE)
        )
        message("\t\t Reading")
        x <- Rsamtools::scanBam(bamFile, index = bamFile, param = params)[[1]]
        chromread <- as.character(unique(x[[2]]))

        if (!isTRUE(all.equal(length(chromread), 1)) ||
            !isTRUE(all.equal(chromread, chrom))) {
            stop("Problem in reading the bam file per chromosome")
        }

        message("\t\t Counting")
        freqseq <- table(x[[1]])
        mat <- cbind(
            ID = rownames(freqseq),
            Freq = as.numeric(freqseq)
        )
        valuesvec <- as.numeric(mat[, "Freq"])
        sumstat <- summary(valuesvec, digits = 4)
        titlehist <- paste(names(sumstat), sumstat, collapse = "-")

        if (!file.exists(outfold)) {
            dir.create(outfold, recursive = TRUE)
        }

        message("\t\t Plotting")
        png(filename = file.path(outfold, paste0("chr", chrom, ".png")))
        hist(valuesvec,
            main = titlehist,
            xlab = paste0("Nb of matches on chr", chrom),
            ylab = "Nb of sequences", breaks = 1000
        )
        dev.off()

        png(filename = file.path(
            outfold,
            paste0("chr", chrom, "-limitedQuart.png")
        ))
        hist(valuesvec,
            main = titlehist,
            xlab = paste0("Nb of matches on chr", chrom),
            ylab = "Nb of sequences", breaks = 1000,
            xlim = c(0, sumstat[5] + 40)
        )
        dev.off()

        return(freqseq)
    }, outputfold1, mc.cores = ncores)

    names(freqseqlist) <- chromvec

    ## Retrieving the names of all sequences
    seqnamesvec <- unique(unlist(lapply(freqseqlist, names)))
    ## Retrieving the count for each sequence on each chromosome
    reslist <- parallel::mclapply(freqseqlist, function(x, allnames) {
        res <- x[allnames]
        idxNA <- which(is.na(res))
        if (!isTRUE(all.equal(length(idxNA), 0))) {
            res[idxNA] <- 0
        }
        return(res)
    }, seqnamesvec, mc.cores = ncores)
    ## Verify that all elements of list have the same number of counts
    if (!isTRUE(all.equal(length(unique(lengths(reslist))), 1))) {
        stop("Problem in retrieving the counts for each chromosome")
    }

    ## Building a matrix with number of matches per chromosomes
    matfreqperchrom <- do.call("cbind", reslist)

    ## Several operations:
    ## 1) Count nb of sequences having matches on several chromosomes
    ## 2) Building a vector of counts for each sequence
    countsreslist <- apply(matfreqperchrom, 1, function(currentseqfreq) {
        oneormorechrom <- uniqueormultichrom(currentseqfreq)
        allmatches <- sum(currentseqfreq)
    }, simplify = FALSE)

}, bamvec, namesbamvec, MoreArgs = list(expname, chromvec, ncores, outputfold))
