# The function getLength() calculates the length (nt) of the
# back-spliced exons and the corresponding flanking introns.
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
            meanLengthBSEs = base::rowMeans(base::cbind(.data$lenUpBSE, .data$lenDownBSE), na.rm =
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
