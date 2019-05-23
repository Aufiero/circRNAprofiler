# Satisfy R CMD check
if (getRversion() >= "3.1.0")
    utils::globalVariables(c("ebicat38", "ebicat37", "gwascat", "."))


#' @title Format annotation file
#'
#' @description The function formatGTF() formats the given annotation file.
#'
#' @param pathToGTF A string containing the path to the the GTF file.
#' Use the same annotation file used during the RNA-seq mapping procedure.
#'
#' @return A data frame.
#'
#' @examples
#' gtf <- formatGTF()
#'
#' @importFrom stringr str_extract
#' @importFrom rtracklayer import
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @import dplyr
#'
#' @export
formatGTF <- function(pathToGTF = NULL) {
    if (!is.null(pathToGTF)) {
        # suppress warning for closing unused connection
        gtf <- suppressWarnings(rtracklayer::import(pathToGTF)) %>%
            as.data.frame()

        # Keep only the needed rows
        # Rename seqnames to chrom
        gtf <- gtf %>%
            dplyr::filter(.data$type == "exon") %>%
            dplyr::rename(chrom = .data$seqnames) %>%
            dplyr::mutate(strand = as.character(.data$strand),
                chrom = as.character(.data$chrom))

        needColumns <-
            c(
                "chrom",
                "start",
                "end",
                "width",
                "strand",
                "type",
                "gene_name",
                "transcript_id"
            )

        # For UCSC and Ncbi the exon_number column needs to be added
        if (!("exon_number" %in% colnames(gtf))) {
            pos <-  gtf  %>%
                dplyr::filter(.data$strand == "+") %>%
                dplyr::group_by(.data$transcript_id) %>%
                dplyr::arrange(.data$start, .by_group = TRUE) %>%
                dplyr::mutate (exon_number = row_number())

            neg <-  gtf  %>%
                dplyr::filter(.data$strand == "-") %>%
                dplyr::group_by(.data$transcript_id) %>%
                dplyr::arrange(.data$start, .by_group = TRUE) %>%
                dplyr::mutate (exon_number = rev(row_number()))

            formattedGTF <- rbind(pos, neg)
            formattedGTF <- formattedGTF %>%
                dplyr::select(c(needColumns), "exon_number") %>%
                dplyr::mutate(exon_number = as.numeric(.data$exon_number)) %>%
                as.data.frame()

        } else {
            formattedGTF <- gtf %>%
                dplyr::select(c(needColumns), "exon_number") %>%
                dplyr::mutate(exon_number = as.numeric(.data$exon_number)) %>%
                as.data.frame()
        }

        # For ensamble and ncbi the "chr" needs to be added to the chromosome number
        if (is.na(stringr::str_extract(formattedGTF$chrom[1], "chr"))) {
            formattedGTF$chrom <- paste0("chr", formattedGTF$chrom)
        }
    } else{
        formattedGTF <- data.frame()
        cat("Specify path to GTF file.")
    }

    return(formattedGTF)
}
