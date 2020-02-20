#' @title Differential circRNA expression analysis adapted from DESeq2
#'
#' @description The helper functions getDeseqRes() identifies
#' differentially expressed circRNAs. The latter uses respectively the R
#' Bioconductor packages DESeq2 which implements a beta-binomial model to
#' model changes in circRNA expression.
#'
#' @param backSplicedJunctions A data frame containing the back-spliced junction
#' coordinates and counts in each analyzed sample.
#' See \code{\link{getBackSplicedJunctions}} and \code{\link{mergeBSJunctions}}
#' (to group circRNA detected by multiple detection tools) on how to generated
#' this data frame.
#'
#' @param condition A string specifying which conditions to compare. Only 2
#' conditions at the time can be analyzed. Separate the 2 conditions with a
#' dash, e.g. A-B. Use the same name used in column condition in experiment.txt.
#' log2FC calculation is perfomed by comparing the condition positioned
#' forward against the condition positioned backward in the alphabet.
#' E.g. if there are 2 conditions A and B then a negative log2FC means that
#' in condition B there is a downregulation of the corresponding circRNA.
#' If a positive log2FC is found means that there is an upregulation in
#' condition B of that circRNA.
#'
#' @param pAdjustMethod A character string stating the method used to adjust
#' p-values for multiple testing. See \code{\link[stats]{p.adjust}}.
#' Deafult value is "BH".
#'
#' @param ... Arguments to be passed to the DESeq function used internally from
#' DESeq2 package. If nothing is specified the default
#' values of the function \code{\link[DESeq2]{DESeq}} are used.
#'
#' @param pathToExperiment A string containing the path to the experiment.txt
#' file. The file experiment.txt contains the experiment design information.
#' It must have at least 3 columns with headers:
#' \describe{
#' \item{label:}{(1st column) - unique names of the samples (short but informative).}
#' \item{fileName:}{(2nd column) - name of the input files - e.g. circRNAs_X.txt, where
#' x can be can be 001, 002 etc.}
#' \item{group:}{ (3rd column) - biological conditions - e.g. A or B; healthy or diseased
#' if you have only 2 conditions.}
#' }
#'
#' By default pathToExperiment is set to NULL and the file it is searched in
#' the working directory. If experiment.txt is located in a different directory
#' then the path needs to be specified.
#'
#' @return A data frame.
#'
#' @examples
#' # Load a data frame containing detected back-spliced junctions
#' data("mergedBSJunctions")
#'
#' pathToExperiment <- system.file("extdata", "experiment.txt",
#'     package ="circRNAprofiler")
#'
#' # Filter circRNAs
#' filteredCirc <- filterCirc(
#'     mergedBSJunctions,
#'     allSamples = FALSE,
#'     min = 5,
#'     pathToExperiment)
#'
#' # Find differentially expressed circRNAs
#' deseqResBvsA <- getDeseqRes(
#'     filteredCirc,
#'     condition = "A-B",
#'     pAdjustMethod = "BH",
#'     pathToExperiment)
#'
#' @import dplyr
#' @import DESeq2
#' @importFrom utils read.table
#' @importFrom rlang .data
#'@export
getDeseqRes <-
    function(backSplicedJunctions,
        condition,
        pAdjustMethod = "BH",
        pathToExperiment = NULL,
        ...) {
        # Read experiment.txt
        experiment <- .readExperiment(pathToExperiment)
        if (nrow(experiment)) {
            cond <- strsplit(condition, "-")[[1]]
            experiment <-
                experiment[experiment$condition %in% cond,]
            experiment$condition <-
                as.factor(experiment$condition)

            # Analysis with DESeq2
            dds <-
                DESeq2::DESeqDataSetFromMatrix(
                    countData = backSplicedJunctions[, experiment$label],
                    colData =
                        experiment[,-which(names(experiment) %in% "fileName")],
                    design = ~ condition
                )

            # Differential expression analysis - standard analysis
            dds <- DESeq2::DESeq(dds, ...)
            statistics <- DESeq2::results(dds, pAdjustMethod = pAdjustMethod)

            # Get deseqRes data frame
            deseqRes <- .getDeseqResDF(backSplicedJunctions,
                statistics,
                dds)
        } else{
            deseqRes <- data.frame()
            cat("experiment.txt not found in wd (or empty). Differential expression analysis can
                not be run. Type ?getDeseqRes and see pathToExperiment param.\n")
        }
        return(deseqRes)
    }


# get deseqRes data frame
.getDeseqResDF <-
    function(backSplicedJunctions,
        statistics, dds) {
        deseqRes <-
            dplyr::bind_cols(
                data.frame(
                    backSplicedJunctions$id,
                    backSplicedJunctions$gene,
                    statistics,
                    counts(dds, norm = TRUE)
                )
            ) %>%
            dplyr::select(-c(.data$baseMean, .data$lfcSE, .data$stat)) %>%
            dplyr::rename(
                log2FC = .data$log2FoldChange,
                gene = .data$backSplicedJunctions.gene,
                id = .data$backSplicedJunctions.id
            ) %>%
            dplyr::mutate(id = as.character(.data$id),
                gene = as.character(.data$gene))
        return(deseqRes)
    }

# If the function you are looking for is not here check supportFunction.R
# Functions in supportFunction.R are used by multiple functions.
