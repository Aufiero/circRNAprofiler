#' @title Return column names
#'
#' @description The function getBasicColNames() returns the basic column names.
#'
#' @return A character vector
#'
#' @examples
#' # Inner function
#' getBasicColNames()
#'
#' @export
getBasicColNames <- function() {
    basicColumns <- c("id",
                      "gene",
                      "strand",
                      "chrom",
                      "startUpBSE", # back-spliced junction
                      "endDownBSE") # back-spliced junction
    return(basicColumns)

}


#' @title Check backSplicedJunctions data frame
#'
#' @description The function checkBSJsDF() verifies that the functions to import
#' the detected circRNAs generate the correct data structure and content.
#'
#' @param backSplicedJunctions A data frame containing back-spliced junction
#' coordinates.
#'
#' @param addColNames A string vector containing the columns that the data frame
#' must contains in addition to the basic columns. By default addColNames is set
#' to NULL.
#'
#' @return A data frame
#'
#' @importFrom magrittr %>%
#' @import dplyr
#'
#' @examples
#' # Load short version of the gencode v19 annotation file
#' data("gtf")
#'
#' # Path to an example file containing circRNA detected by MapSplice2
#' pathToFile <- system.file("extdata", "mapsplice/circRNAs_001.txt",
#'     package="circRNAprofiler")
#'
#' # Inner function.
#' # Import circRNAs.
#' backSplicedJunctions <- importMapSplice(pathToFile, gtf)
#'
#' # Inner function
#' # check table
#' checkBSJsDF(backSplicedJunctions, addColNames = "coverage")
#'
#' @import magrittr
#' @export
checkBSJsDF <- function(backSplicedJunctions, addColNames = NULL) {
    colNames <- c(getBasicColNames(), addColNames)

    if (!all(colNames  %in% colnames(backSplicedJunctions))) {
        missingNamesId <- which(!colnames(backSplicedJunctions) %in% colNames)
        stop("missing or wrong column names: ",
             paste(colnames(backSplicedJunctions)[missingNamesId],
                   collapse = " \t"))
    }

    for (i in seq_along(backSplicedJunctions$strand)) {
        # For gene located on the positive strand the coordinates of an
        # acceptor site of an exon must be smaller than the coordinates of
        # the donor site
        if (backSplicedJunctions$strand[i] == "+" &
            backSplicedJunctions$startUpBSE[i] > backSplicedJunctions$endDownBSE[i]) {
            stop(
                "found startUpBSE > endDownBSE - for gene on the positive
                strand the coordinates of an startUpBSE (acceptor site)
                of an exon must be smaller than the coordinates
                of the endDownBSE(donor site)"
            )
        }

        # For gene located on the negative strand the coordinates of an
        # acceptor site of an exon must be greater than the coordinates
        # of the donor site
        if (backSplicedJunctions$strand[i] == "-" &
            backSplicedJunctions$startUpBSE[i] < backSplicedJunctions$endDownBSE[i]) {
            stop(
                "found startUpBSE < endDownBSE - for gene on the negative
                strand the coordinates of an startUpBSE (acceptor site)
                of an exon must be greater than the coordinates
                of the endDownBSE (donor site)"
            )
        }
    }
    return(backSplicedJunctions)
}


#' @title Return column names
#'
#' @description The function getTargetsColNames() returns the column names.
#'
#' @return A character vector.
#'
#' @examples
#' # Inner function
#' getTargetsColNames()
#'
#' @export
getTargetsColNames <- function() {
    targetsColumns <- c(
        "id",
        "gene",
        "transcript",
        "strand",
        "chrom",
        "startGR",
        "endGR",
        "length",
        "seq",
        "type"
    )
    return(targetsColumns)

}
