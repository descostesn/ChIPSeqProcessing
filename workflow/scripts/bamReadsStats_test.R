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


##################
# PARAMETERS
##################

bamvec <- c("/g/romebioinfo/Projects/TEbench/results/bam/paired/bowtie2_results/mm39/H3K9me3_SRR2136759_GSM1841035_GSE71589_ESC_bio2_trimmed_k150_sorted_noDups.bam")
expname <- "H3K9me3_SRR2136759_GSM1841035_GSE71589_ESC_bio2"
namesbamvec <- c("k150")
chromvec <- "9"

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

fasttable <- function(x) {
    xf <- collapse::qF(x)
    levelvec <- levels(xf)
    levels(xf) <- seq_len(length(levelvec))
    res <- tabulate(xf)
    names(res) <- levelvec
    return(res)
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

bamfile <- bamvec
bamname <- namesbamvec
chrom <- chromvec
message("Reading ", expname, "-", bamname)
x <- readingseqonchrom(chrom, bamfile)
message("Done")
