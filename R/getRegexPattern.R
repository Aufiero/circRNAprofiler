#' @title Convert IUPAC sequence to an regular expression
#'
#' @description The function getRegexPattern() converts a nucleotide sequence
#' with IUPAC codes to an regular expression.
#'
#' @param iupacSeq A character string with IUPAC codes.
#'
#' @param isDNA A logical specifying whether the iupacSeq is DNA or RNA.
#' Deafult value is FALSE.
#'
#' @return A character string.
#'
#' @examples
#' regextPattern <- getRegexPattern("CGUKMBVNN", isDNA = FALSE)
#'
#' @importFrom utils data
#' @export
getRegexPattern <- function(iupacSeq, isDNA = FALSE) {
    # Convert the sequence to UPPER CASE
    iupacSeq <- toupper(iupacSeq)

    if (isDNA) {
        col <- paste0("regex", "DNA")
    } else{
        col <- paste0("regex", "RNA")
    }

    # get IUPAC code and the corresponding regular expressions
    iupac <- NULL
    data("iupac", package= "circRNAprofiler", envir = environment())
    mt <-
        match(base::strsplit(iupacSeq, "")[[1]], iupac$code)
    regexPattern <- paste(iupac[mt, col], collapse = "")

    return(regexPattern)
}


# If the function you are looking for is not here check supportFunction.R
# Functions in supportFunction.R are used by multiple functions.
