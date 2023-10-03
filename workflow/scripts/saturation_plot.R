###########
# Script that generates the saturation plot of multireads for single end exp.
# Descostes Dec 2019
###########

library(stringr)


#############
## PARAMS
#############

logbest <- snakemake@input$logVecBest
logmulti <- snakemake@input$logVecMulti
theoretical_suffixes <- c("best", "k10",  "k50",  "k100", "k150", "k200",
    "k250", "k300", "k350", "k400")


#############
## FUNCTIONS
#############

saturation_study <- function(logbest, logmulti, output_folder) {

    # Grouping log files per experiment
    alllog <- c(logbest, logmulti)
    factorexp <- as.factor(basename(unlist(lapply(strsplit(alllog, "_trimmed"),
        "[", 1))))
    listperexp <- split(alllog, factorexp)

    # Creation of the list of multiread percentages for each experiment
    percentlist <- lapply(listperexp, function(expfilevec) {

        percentvec <- sapply(expfilevec, function(expfile) {
            stats <- readLines(expfile)
            if (!isTRUE(all.equal(length(stats), 6)))
                stop("In the script saturation_plot.R which is in the rule ",
                    "saturation_plot_*, the file ", expfile, " does not ",
                    "contain 6 lines")

            if (isTRUE(all.equal(snakemake@params$type, "single")))
                patternmulti <- "aligned >1 times"
            else
                patternmulti <- "aligned concordantly >1 times"

            if (!isTRUE(all.equal(length(grep(patternmulti, stats[5])), 1)))
                stop("In the script saturation_plot.R which is in the rule ",
                    "saturation_plot_*, the file ", expfile, " does not ",
                    "contain multi nb on line 5")

            result <- strsplit(str_extract(stats[5], "[0-9]+\\.[0-9]+\\%"),
                "%")[[1]]
            return(result)})

        return(as.numeric(percentvec))})

    # Plotting the saturation curve for each experiment
    invisible(mapply(function(percentvec, outputname) {
        png(filename = paste0(output_folder, outputname, ".png"))
        barplot(percentvec,
                names.arg=c("best", "10", "50", "100", "150", "all"),
                xlab="Maximum reportable alignments",
                ylab = "Proportions of multi-reads")
        dev.off()
    }, percentlist, names(percentlist)))

    # Recording percent numbers for each exp
    suffixlist <- lapply(listperexp, function(expvec) {
        return(unlist(strsplit(unlist(lapply(strsplit(expvec, "trimmed_"),
            "[", 2)), "\\.log")))})

    suffixvec <- unique(unlist(suffixlist))

    if (!isTRUE(all.equal(theoretical_suffixes, suffixvec)))
        stop("In the script saturation_plot.R which is in the rule ",
            "saturation_plot_*, there is a problem with suffixes retrieval")

    invisible(mapply(function(percentvec, outputname, knames) {
        percentmat <- t(as.matrix(percentvec))
        colnames(percentmat) <- knames
        outfold <- paste0(output_folder, outputname, ".txt")
        write.table(percentmat, file = outfold, sep="\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)}, percentlist,
            names(percentlist), MoreArgs = list(suffixvec)))
}



#############
## MAIN
#############

cat("Generating saturation plots for ", snakemake@params$type, " end\n")
saturation_study(logbest, logmulti,
    paste0("../results/qc/bowtie2_saturation_multireads/",
        snakemake@params$type, "/"))
