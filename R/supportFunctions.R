# The function getBasicColNames() returns the basic column names.
getBasicColNames <- function() {
    basicColumns <- c("id",
        "gene",
        "strand",
        "chrom",
        "startUpBSE", # back-spliced junction
        "endDownBSE") # back-spliced junction
    return(basicColumns)

}


# The function checkBSJsDF() verifies that the functions to import
# the detected circRNAs generate the correct data structure and content.
checkBSJsDF <- function(backSplicedJunctions, addColNames = NULL) {
    colNames <- c(getBasicColNames(), addColNames)

    if (!all(colNames  %in% colnames(backSplicedJunctions))) {
        missingNamesId <-
            which(!colnames(backSplicedJunctions) %in% colNames)
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



# The function getTargetsColNames() returns the column names.
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
