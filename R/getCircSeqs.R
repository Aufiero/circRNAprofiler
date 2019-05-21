#' @title Retrieve circRNA sequences
#'
#' @description The function getCircSeqs() retrieves the circRNA
#' sequences. The circRNA sequence is given by the sequences of the exons
#' in between the back-spliced-junctions.
#'
#' @param annotatedBSJs A data frame with the annotated back-spliced junctions.
#' It can be generated with \code{\link{annotateBSJs}}.
#'
#' @param gtf A data frame containing genome annotation information. It can be
#' generated with \code{\link{formatGTF}}.
#'
#' @param bsGenome A BSgenome object from which to retrieve the sequences.
#' It can be generated with \code{\link[BSgenome]{getBSgenome}}.
#' See available.genomes() to see the BSgenome package currently available.
#'
#' @return A list.
#'
#' @examples
#' # Load a data frame containing detected back-spliced junctions
#' data("mergedBSJunctions")
#'
#' # Load short version of the gencode v19 annotation file
#' data("gtf")
#'
#' # Annotate the first 10 back-spliced junctions
#' annotatedBSJs <- annotateBSJs(mergedBSJunctions[1:10, ],gtf)
#'
#' # Load BSgenome object
#' bsGenome <- BSgenome::getBSgenome("BSgenome.Hsapiens.UCSC.hg19")
#'
#' # Retrieve target sequences
#' targets <- getCircSeqs(
#'     annotatedBSJs,
#'     gtf,
#'     bsGenome)
#'
#' @importFrom BSgenome getSeq
#' @importFrom Biostrings reverseComplement
#' @importFrom Biostrings RNAString
#' @importFrom rlang .data
#' @import dplyr
#' @export
getCircSeqs <-
    function(annotatedBSJs,
        gtf,
        bsGenome) {
        # Create a enmpty list of 1 elements
        targets <- vector("list", 1)
        names(targets)[1] <- "circ"

        # Create an empty data frame
        targets[[1]] <-
            data.frame(matrix(
                nrow = nrow(annotatedBSJs),
                ncol = length(getTargetsColNames())
            ))
        colnames(targets[[1]]) <- getTargetsColNames()

        # Fill the data frame with the needed information
        targets[[1]]$id <- annotatedBSJs$id
        targets[[1]]$gene <- annotatedBSJs$gene
        targets[[1]]$transcript <- annotatedBSJs$transcript
        targets[[1]]$chrom <- annotatedBSJs$chrom
        targets[[1]]$length <- annotatedBSJs$lenCircRNA
        targets[[1]]$strand <- annotatedBSJs$strand
        targets[[1]]$type <- rep("circ", nrow(targets[[1]]))

        for (k in seq_along(targets[[1]]$id)) {
            if (targets[[1]]$strand[k] == "+") {
                targets[[1]]$startGR[k] <- annotatedBSJs$startUpBSE[k]
                targets[[1]]$endGR[k] <- annotatedBSJs$endDownBSE[k]
            } else if (targets[[1]]$strand[k] == "-") {
                targets[[1]]$startGR[k] <- annotatedBSJs$endDownBSE[k]
                targets[[1]]$endGR[k] <- annotatedBSJs$startUpBSE[k]
            }
        }

        # Retrieve the sequences
        for (i in seq_along(annotatedBSJs$id)) {
            # If the transcript is not present we do not have the coordinates to
            # retrieve the seq
            if (!is.na(annotatedBSJs$transcript[i])) {
                exonsToSelect <-
                    annotatedBSJs$exNumUpBSE[i]:annotatedBSJs$exNumDownBSE[i]

                transcript <- gtf %>%
                    dplyr::filter(
                        .data$type == "exon",
                        .data$gene_name == annotatedBSJs$gene[i],
                        .data$transcript_id == annotatedBSJs$transcript[i],
                        .data$exon_number %in% exonsToSelect
                    )

                exonSeqs <- as.character()

                if (transcript$strand[1] == "-") {
                    transcript <- transcript  %>%
                        dplyr::arrange(desc(.data$exon_number))

                    # Note: the sequences retrieved by getSeq function
                    # from BSgenome correspond to the positive strand of
                    # the DNA (genome reference).For a circRNA that arises
                    # from a gene transcribed from negative strand we need to
                    # take the reverse complement of the sequence. Complement
                    # because the sense strand is the negative strand of the
                    # DNA. Reverse because we always report the sequences from
                    # 5' to 3'in the final data frame.

                    for (j in seq_along(transcript$exon_number)) {
                        exonSeqs[j] <- gsub("T",
                            "U",
                            as.character(
                                BSgenome::getSeq(
                                    bsGenome,
                                    names = transcript$chrom[j],
                                    transcript$start[j],
                                    transcript$end[j]
                                )
                            ))

                    }

                    joinedExonSeqs <- paste(exonSeqs, collapse = "")

                    # For the negative strand the last 40 nucleotide of the
                    # upstream back-spliced exons are attached at the end to
                    # the downstream back-spliced exon to recreate the
                    # back-spliced sequence. In this way the mir analysis
                    # can be perfomed also across the back-spliced
                    # junctions.

                    joinedExonSeqsWithBSJ <-
                        paste0(substr(
                            joinedExonSeqs,
                            nchar(joinedExonSeqs) - 49,
                            nchar(joinedExonSeqs)
                        ),
                            joinedExonSeqs)


                    targets[[1]]$seq[i] <-
                        as.character(reverseComplement(Biostrings::RNAString(joinedExonSeqsWithBSJ)))

                } else if (transcript$strand[1] == "+") {
                    transcript <- transcript  %>%
                        arrange(.data$exon_number)
                    # For the positive strand no modification is needed because
                    # the sense strand corresponds to the positive strand of
                    # the DNA that is the genome reference.
                    for (j in seq_along(transcript$exon_number)) {
                        exonSeqs[j] <- gsub("T",
                            "U",
                            as.character(
                                BSgenome::getSeq(
                                    bsGenome,
                                    names = transcript$chrom[j],
                                    transcript$start[j],
                                    transcript$end[j]
                                )
                            ))
                    }

                    joinedExonSeqs <- paste(exonSeqs, collapse = "")

                    # For the positive strand the first 50 nucleotide of the
                    # upstream back-spliced exon are taken and attached at
                    # the end of the downstream back-spliced exon.
                    joinedExonSeqsWithBSJ <-
                        paste0(joinedExonSeqs, substr(joinedExonSeqs, 1, 50))


                    targets[[1]]$seq[i] <- joinedExonSeqsWithBSJ

                }
            }

        }

        return(targets)

        }
