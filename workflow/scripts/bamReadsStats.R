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

readingseqonchrom <- function(chrom, bamfile) {
    message("\t\t Reading on chr ", chrom)

    ## Retrieving the sequence name (rname) on the defined chrom (qname)
    params <- Rsamtools::ScanBamParam(
        what = c("rname", "qname"),
        which = GRanges(chrom, IRanges(1, 1e8)),
        flag = scanBamFlag(isUnmappedQuery = FALSE)
    )
    x <- Rsamtools::scanBam(bamfile, index = bamfile, param = params)[[1]]

    ## Sanity check
    chromread <- as.character(unique(x[[2]]))
    if (!isTRUE(all.equal(length(chromread), 1)) ||
        !isTRUE(all.equal(chromread, chrom))) {
        stop("Problem in reading the bam file per chromosome")
    }

    return(x)
}

computefreqonchrom <- function(x) {
    message("\t\t Counting nb of sequence matches on the chromosome")
    freqseq <- table(x[[1]])
    return(freqseq)
}

plotnbmatchesperchrom <- function(freqseq, outfold, chrom) {
    message("\t\t Plotting")

    if (!file.exists(outfold)) {
        dir.create(outfold, recursive = TRUE)
    }

    mat <- cbind(ID = rownames(freqseq), Freq = as.numeric(freqseq))
    valuesvec <- as.numeric(mat[, "Freq"])
    sumstat <- summary(valuesvec, digits = 4)
    titlehist <- paste(names(sumstat), sumstat, collapse = "-")


    png(filename = file.path(outfold, paste0("chr", chrom, ".png")))
    hist(valuesvec,
        main = titlehist,
        xlab = paste0("Nb of matches on chr", chrom),
        ylab = "Nb of sequences", breaks = 1000
    )
    dev.off()

    png(filename = file.path(outfold, paste0("chr", chrom, "-limitedQuart.png"))) # nolint
    hist(valuesvec,
        main = titlehist,
        xlab = paste0("Nb of matches on chr", chrom),
        ylab = "Nb of sequences", breaks = 1000,
        xlim = c(0, sumstat[5] + 40)
    )
    dev.off()
}

buildingmatmatchesperchrom <- function(freqseqlist, ncores) {
    message("\t\t Building matrix")

    ## Retrieving the names of all sequences
    seqnamesvec <- unique(unlist(lapply(freqseqlist, names)))

    ## Retrieving the count for each sequence on each chromosome
    reslist <- parallel::mclapply(freqseqlist, function(x, allnames) {
        res <- x[allnames]
        idxna <- which(is.na(res))
        if (!isTRUE(all.equal(length(idxna), 0))) {
            res[idxna] <- 0
        }
        return(res)
    }, seqnamesvec, mc.cores = ncores)

    ## Verify that all elements of list have the same number of counts
    if (!isTRUE(all.equal(length(unique(lengths(reslist))), 1))) {
        stop("Problem in retrieving the counts for each chromosome")
    }

    ## Building a matrix with number of matches per chromosomes
    matfreqperchrom <- do.call("cbind", reslist)

    return(matfreqperchrom)
}

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

labelsandmismatches <- function(ncores, matfreqperchrom) {
    message("\t\t Computing occupancy and mismatches")
    cl <- makeCluster(ncores)
    labeloccupancyvec <- parRapply(cl, matfreqperchrom, uniqueormultichrom)
    nbmatchseqvec <- parRapply(cl, matfreqperchrom, sum)
    stopCluster(cl)
    return(list(labeloccupancyvec, nbmatchseqvec))
}

plotcounts <- function(outputfold1, tablabeloccupancy, nbmatchseqvec) {
    message("\t\t Plotting occupancy and mismatches")
    png(file = file.path(outputfold1, "nbchrommatchesbyseq.png"))
    barplot(tablabeloccupancy,
        xlab = "Number of chromosomes with match per sequence",
        ylab = "Number of sequences")
    dev.off()

    png(file = file.path(outputfold1, "nbseqmatches.png"))
    hist(nbmatchseqvec, breaks = 1000, xlab = "Number of matches",
        ylab = "Number of sequences")
    dev.off()
}

##################
# MAIN
##################

nbmatchseqlist <- mapply(function(bamfile, bamname, expname, chromvec,
    ncores, outputfold) {

    message("Reading ", expname, "-", bamname)
    outputfold1 <- file.path(outputfold, expname, bamname)

    ## Counting number of matches per chromosome
    freqseqlist <- parallel::mclapply(chromvec, function(chrom, outfold) {
        x <- readingseqonchrom(chrom, bamfile)
        freqseq <- computefreqonchrom(x)
        plotnbmatchesperchrom(freqseq, outfold, chrom)
        return(freqseq)
    }, outputfold1, mc.cores = ncores)
    names(freqseqlist) <- chromvec

    ## Building a matrix of the number of matches of each sequence on each
    ## chromosome
    matfreqperchrom <- buildingmatmatchesperchrom(freqseqlist, ncores)

    ## Several operations:
    ## 1) Count nb of sequences having matches on several chromosomes
    ## 2) Building a vector of counts for each sequence
    countsreslist <- labelsandmismatches(ncores, matfreqperchrom)
    tablabeloccupancy <- table(countsreslist[[1]])
    nbmatchseqvec <- countsreslist[[2]]
    plotcounts(outputfold1, tablabeloccupancy, nbmatchseqvec)
    gc(verbose = FALSE)
    return(nbmatchseqvec)
}, bamvec, namesbamvec, MoreArgs = list(expname, chromvec, ncores, outputfold))


message("Retrieving nb of sequences with max nb of matches for all bam files")
maxnbmatchvec <- mapply(function(nbmatchseq, bamname) {

    message("\t Computing frequencies for ", bamname)
    tabnbmatchseq <- table(nbmatchseq)
    return(tabnbmatchseq[which.max(tabnbmatchseq)])
}, nbmatchseqlist, namesbamvec, SIMPLIFY = TRUE)
names(maxnbmatchvec) <- namesbamvec

png(file = file.path(output_folder, expname, "nbseqperthreshold.png"))
barplot(maxnbmatchvec, xlab = "Bam threshold", ylab = "Nb of sequences")
dev.off()
