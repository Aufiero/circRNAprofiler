#' @title Retrieve random back-spliced junctions
#'
#' @description The function getRandomBSJunctions() retrieves random
#' back-spliced junctions from the user genome annotation.
#'
#' @param gtf A dataframe containing genome annotation information This can be
#' generated with \code{\link{formatGTF}}.
#'
#' @param n Integer specifying the number of randomly selected transcripts
#' from which random back-spliced junctions are extrated. Default value = 100.
#'
#' @param f An integer specifying the fraction of single exon circRNAs that
#' have to be present in the output data frame. Default value is 10.
#'
#' @return A data frame.
#'
#' @examples
#' # Load short version of the gencode v19 annotation file
#' data("gtf")
#'
#' # Get 10 random back-spliced junctions
#' randomBSJunctions <- getRandomBSJunctions(gtf, n = 10, f = 10)
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @import dplyr
#' @export
getRandomBSJunctions <- function(gtf, n = 100, f = 10) {
    basicColumns <- getBasicColNames()

    # Create an empty data frame
    randomBSJunctions <-
        data.frame(matrix(nrow = n, ncol = length(basicColumns)))
    colnames(randomBSJunctions) <- basicColumns

    # calculate the percentage of back-spliced junctions from sigle exons
    c <- round( (n / 100) * f, 0)

    # Select one random exon from n randomly selected transcript
    bsExons1 <- gtf %>%
        dplyr::filter(.data$type == "exon") %>%
        dplyr::group_by(.data$transcript_id) %>%
        dplyr::filter(n() >= 3) %>%
        dplyr::ungroup() %>%
        dplyr::filter(.data$transcript_id %in% sample(
            unique(.data$transcript_id), c)) %>%
        dplyr::group_by(.data$transcript_id) %>%
        dplyr::sample_n(1) %>%
        dplyr::arrange(.data$exon_number)

    # If only one exon is picked then the row is duplicated. This step is only
    # made to simply the code below without introducing an additional if
    # statment when a single exon is present in the isoform sampled in
    # the above step.
    bsExons1 <- dplyr::bind_rows(bsExons1, bsExons1)

    # Select two random back-spliced exons (up and down) from n randomly
    # selected transcript
    bsExons2 <- gtf %>%
        dplyr::filter(.data$type == "exon") %>%
        dplyr::group_by(.data$transcript_id) %>%
        dplyr::filter(n() >= 3) %>%
        dplyr::ungroup() %>%
        dplyr::filter(.data$transcript_id %in%
                sample(unique(.data$transcript_id), n - c)) %>%
        dplyr::group_by(.data$transcript_id) %>%
        dplyr::sample_n(2) %>%
        dplyr::arrange(.data$exon_number)

    # Join the 2 data frame bsExons1 and bsExons2
    bsExons12 <-  dplyr::bind_rows(bsExons1, bsExons2)

    # Align duplicates on the same row
    indexDup <- which(duplicated(bsExons12$transcript_id))
    duplicates <- bsExons12[indexDup, ]
    cleanedDF <- bsExons12[-indexDup, ]
    mt <- match(cleanedDF$transcript_id, duplicates$transcript_id)

    allBSEs <- dplyr::bind_cols(cleanedDF, duplicates[mt, ])

    # For negative strand
    allBSEsNeg <- allBSEs[allBSEs$strand == "-", ]
    if (nrow(allBSEsNeg) > 0) {
        randomBSJunctions$gene[seq_along(allBSEsNeg$exon_number)] <-
            allBSEsNeg$gene_name
        randomBSJunctions$strand[seq_along(allBSEsNeg$exon_number)] <-
            allBSEsNeg$strand
        randomBSJunctions$chrom[seq_along(allBSEsNeg$exon_number)] <-
            allBSEsNeg$chrom
        randomBSJunctions$startUpBSE[seq_along(allBSEsNeg$exon_number)] <-
            allBSEsNeg$end
        randomBSJunctions$endDownBSE[seq_along(allBSEsNeg$exon_number)] <-
            allBSEsNeg$start1
    }

    # For positive strand
    allBSEsPos <- allBSEs[allBSEs$strand == "+", ]
    if (nrow(allBSEsPos) > 0) {
        randomBSJunctions$gene[(nrow(allBSEsNeg) + 1):n] <-
            allBSEsPos$gene_name
        randomBSJunctions$strand[(nrow(allBSEsNeg) + 1):n] <-
            allBSEsPos$strand
        randomBSJunctions$chrom[(nrow(allBSEsNeg) + 1):n] <-
            allBSEsPos$chrom
        randomBSJunctions$startUpBSE[(nrow(allBSEsNeg) + 1):n] <-
            allBSEsPos$start
        randomBSJunctions$endDownBSE[(nrow(allBSEsNeg) + 1):n] <-
            allBSEsPos$end1
    }

    # Generate a unique identifier by combining the values of the following
    # columns: gene, strand, chrom, endDownBSE and startUpBSE.
    # The values are separated by a semicolumns (:)
    randomBSJunctions$id <- paste(
        randomBSJunctions$gene,
        randomBSJunctions$strand,
        randomBSJunctions$chrom,
        randomBSJunctions$startUpBSE,
        randomBSJunctions$endDownBSE,
        sep = ":"
    )


    # Return the data frame with the genimic coordinates of randomly selected
    # back-spliced exon juntions
    return(randomBSJunctions)

}
