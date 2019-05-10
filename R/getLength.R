#' @title Retieve feature length
#'
#' @description The function getLength() calculates the length (nt) of the
#' back-spliced exons and the corresponding flanking introns.
#'
#' @param annotatedBSJs A data frame with the annotated back-spliced junctions.
#'
#' @return A data frame.
#'
#' @examples
#' # Inner function
#' annotatedBSJs <- data.frame(
#'     matrix(nrow = 1, ncol = length(getAnnotatedBSJsColNames())))
#'
#' colnames(annotatedBSJs) <- getAnnotatedBSJsColNames()
#' annotatedBSJs$id <- "SYCP2:-:chr20:58497514:58497445"
#' annotatedBSJs$startUpIntron <- 58507116
#' annotatedBSJs$endUpIntron <- 58497515
#' annotatedBSJs$startUpBSE <- 58497514
#' annotatedBSJs$endUpBSE <- 58497445
#' annotatedBSJs$startDownBSE <- 58497514
#' annotatedBSJs$endDownBSE <- 58497445
#' annotatedBSJs$startDownIntron <- 58497444
#' annotatedBSJs$endDownIntron <- 58496509
#'
#' # Retrieve length
#' getLength(annotatedBSJs)
#'
#' @import dplyr
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
getLength <- function(annotatedBSJs) {
    colNames <- c(
        "id",
        "lenUpIntron",
        "lenUpBSE",
        "lenDownBSE",
        "lenDownIntron",
        "meanLengthBSEs",
        "meanLengthIntrons"
    )

    # Create an empty dataframe
    exInLength <-
        data.frame(matrix(
            nrow = nrow(annotatedBSJs),
            ncol = length(colNames)
        ))
    colnames(exInLength) <- colNames


    # Calculate the back-spliced exons and introns length and select the needed
    # columns
    exInLength <- annotatedBSJs %>%
        dplyr::mutate(
            lenUpIntron = abs(.data$endUpIntron - .data$startUpIntron),
            lenUpBSE = abs(.data$endUpBSE - .data$startUpBSE),
            lenDownBSE = abs(.data$endDownBSE - .data$startDownBSE),
            lenDownIntron = abs(.data$endDownIntron - .data$startDownIntron),
            meanLengthBSEs = base::rowMeans(
                base::cbind(.data$lenUpBSE, .data$lenDownBSE), na.rm =
                    TRUE),
            meanLengthIntrons = base::rowMeans(
                base::cbind(.data$lenUpIntron, .data$lenDownIntron),
                na.rm =
                    TRUE
            )
        ) %>%
        dplyr::select(
            .data$id,
            .data$lenUpIntron,
            .data$lenUpBSE,
            .data$lenDownBSE,
            .data$lenDownIntron,
            .data$meanLengthBSEs,
            .data$meanLengthIntrons
        )
    # Return a data frame containing the length (bp) of the back-spliced exons
    # and the corresponding falnking introns.

    # Repalce NaN values with NA
    exInLength[is.na(exInLength)] <- NA_character_
    return(exInLength)

}
