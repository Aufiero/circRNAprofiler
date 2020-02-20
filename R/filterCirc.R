#' @title Filter circRNAs
#'
#' @description The functions filterCirc() filters circRNAs on
#' different criteria: condition and read counts. The info reported in
#' experiment.txt file are needed for filtering step.
#'
#' @param backSplicedJunctions A data frame containing back-spliced junction
#' coordinates and counts. See \code{\link{getBackSplicedJunctions}} and
#' \code{\link{mergeBSJunctions}} (to group circRNA detected by multiple
#' detection tools) on how to generated this data frame.
#'
#' @param allSamples A string specifying whether to apply the filter to all
#' samples. Default valu is FALSE.
#'
#' @param min An integer specifying the read counts cut-off.
#' If allSamples = TRUE and min = 0 all circRNAs are kept.
#' If allSamples = TRUE and min = 3, a circRNA is kept if all samples have at
#' least 3 counts. If allSamples = FALSE and min = 2 the filter is applied to
#' the samples of each condition separately meaning that a circRNA is kept if
#' at least 2 counts are present in all sample of 1 of the conditions.
#' Default value is 3.
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
#'
#' @import dplyr
#' @importFrom utils read.table
#' @importFrom magrittr %>%
#'
#'@export
filterCirc <- function(backSplicedJunctions,
    allSamples = FALSE,
    min = 3,
    pathToExperiment = NULL) {

    # Read experiment.txt
    experiment <- .readExperiment(pathToExperiment)
    if (nrow(experiment)) {

        # Get colum names
        colNames <- c(.getBasicColNames(), experiment$label)
        # Create the data frame that will be filled with the circRNA prediction
        # perfomed by the prediction tools used.
        filteredCirc <-
            data.frame(matrix(nrow = 0, ncol = length(colNames)))
        colnames(filteredCirc) <- colNames

        # The filter is applied to all samples
        if (allSamples) {
            filteredCirc <- backSplicedJunctions %>%
                dplyr::filter_at(vars(experiment$label), all_vars(. >= min))

            # The filter is applied to the samples of each condition separately
        } else {
            conditions <- unique(experiment$condition)

            for (i in seq_along(conditions)) {
                cond <- experiment[experiment$condition == conditions[i], "label"]

                backSplicedJunctionsCond <- backSplicedJunctions %>%
                    dplyr::filter_at(vars(cond), all_vars(. >= min))

                filteredCirc <-
                    rbind(filteredCirc, backSplicedJunctionsCond)
            }

            filteredCirc <-
                filteredCirc[!duplicated(filteredCirc),]
        }
    } else{
        filteredCirc <- backSplicedJunctions
        cat("experiment.txt not found in wd (or empty), data frame can not be filtered.
            Type ?filterCirc and see pathToExperiment param.\n")
    }
    # Return a filtered data frame
    return(filteredCirc)
}


# If the function you are looking for is not here check supportFunction.R
# Functions in supportFunction.R are used by multiple functions.
