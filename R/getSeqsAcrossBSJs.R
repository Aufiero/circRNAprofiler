#' @title Retrieve back-spliced junction sequences
#'
#' @description The function getSeqsAcrossBSJs() retrieves
#' the sequences across the back-spliced junctions. A total of 11 nucleotides
#' from each side of the back-spliced junction are taken and concatenated
#' together.
#'
#' @param annotatedBSJs A data frame with the annotated back-spliced junctions.
#' It can be generated with \code{\link{annotateBSJs}}.
#'
#' @param gtf A dataframe containing genome annotation information. It can be
#' generated with \code{\link{formatGTF}}.
#'
#' @param species A string specifying the species from which to retrieve the
#' sequences. For Example: Mmusculus or Hsapiens.
#' See available.genomes() from BSgenome package to see all the species.
#' Default value is "Hsapiens".
#'
#' @param genome A string specifying the genome assembly from which to retrieve
#' the sequences. For Example: mm10 or mm9 for Mouse, hg38 or hg19 for Human.
#' See available.genomes() from BSgenome package to see all genomes.
#' Default value is "hg19".
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
#' annotatedBSJs <- annotateBSJs(mergedBSJunctions[1:10, ], gtf)
#'
#' # Retrieve target sequences
#' targets <- getSeqsAcrossBSJs(
#'     annotatedBSJs,
#'     gtf,
#'     species = "Hsapiens",
#'     genome = "hg19")
#'
#'
#' @importFrom BSgenome getBSgenome
#' @importFrom BSgenome getSeq
#' @importFrom Biostrings RNAString
#' @importFrom Biostrings reverseComplement
#' @importFrom BiocManager install
#' @importFrom rlang .data
#' @import dplyr
#'
#' @export
getSeqsAcrossBSJs <-
    function(annotatedBSJs,
        gtf,
        species = "Hsapiens",
        genome = "hg19") {
        # Create a enmpty list of 1 elements
        targets <- vector("list", 1)
        names(targets)[1] <- "bsj"

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
        targets[[1]]$strand <- annotatedBSJs$strand
        targets[[1]]$type <- rep("bsj", nrow(targets[[1]]))

        for (k in seq_along(targets[[1]]$id)) {
            if (targets[[1]]$strand[k] == "+") {
                targets[[1]]$startGR[k] <- annotatedBSJs$startUpBSE[k]
                targets[[1]]$endGR[k] <- annotatedBSJs$endDownBSE[k]
            } else if (targets[[1]]$strand[k] == "-") {
                targets[[1]]$startGR[k] <- annotatedBSJs$endDownBSE[k]
                targets[[1]]$endGR[k] <- annotatedBSJs$startUpBSE[k]
            }
        }


        # Load required genome from which to retrieve the target sequences
        requiredGenome <-
            paste0("BSgenome.", species, ".", "UCSC", ".", genome)

        # Check if it a correct genome
        if (is.element(requiredGenome, BSgenome::available.genomes())) {
            # It returns FALSE if the package does not exist
            if (!requireNamespace(requiredGenome, quietly = TRUE)) {
                BiocManager::install(requiredGenome)
            }
        } else{
            stop(
                "species name or genome is not correct:
                see available.genomes() from BSgenome package"
            )
        }

        genome <- BSgenome::getBSgenome(requiredGenome)

        # Retrieve the sequences
        for (i in seq_along(annotatedBSJs$id)) {
            # If the transcript is not present we do not have the coordinates
            # to retrieve the seq
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

                    # Note: the sequences retrieved by getSeq function from
                    # BSgenome correspond to the positive strand of the DNA
                    # (genome reference).
                    # For a circRNA that arises from a gene transcribed from
                    # the negative strand we need to take the reverse
                    # complement of the sequence. Complement because the
                    # sense strand is the negative strand of the DNA.
                    # Reverse because we always report the sequences from
                    # 5' to 3' in the final data frame.

                    for (j in seq_along(transcript$exon_number)) {
                        exonSeqs[j] <- gsub("T",
                            "U",
                            as.character(
                                BSgenome::getSeq(
                                    genome,
                                    names = transcript$chrom[j],
                                    transcript$start[j],
                                    transcript$end[j]
                                )
                            ))

                    }

                    joinedExonSeqs <- paste(exonSeqs, collapse = "")

                    # The first 11 nucleotides of the upstream back-spliced
                    # exon are taken and attached at the end of the downstream
                    # back-spliced exon.

                    bsjSeq <-
                        paste0(
                            base::substr(
                                joinedExonSeqs,
                                nchar(joinedExonSeqs) - 10,
                                nchar(joinedExonSeqs)
                            ),
                            base::substr(joinedExonSeqs, 1, 11)
                        )

                    targets[[1]]$length[i] <- nchar(bsjSeq)
                    targets[[1]]$seq[i] <-
                        as.character(Biostrings::reverseComplement(Biostrings::RNAString(bsjSeq)))


                } else if (transcript$strand[1] == "+") {
                    transcript <- transcript  %>%
                        dplyr::arrange(.data$exon_number)
                    # For the positive strand no modification is needed because
                    # the sense strand corresponds to the positive strand of
                    # the DNA that is the genome reference.
                    for (j in seq_along(transcript$exon_number)) {
                        exonSeqs[j] <- gsub("T",
                            "U",
                            as.character(
                                BSgenome::getSeq(
                                    genome,
                                    names = transcript$chrom[j],
                                    transcript$start[j],
                                    transcript$end[j]
                                )
                            ))
                    }

                    joinedExonSeqs <- paste(exonSeqs, collapse = "")

                    # The first 11 nucleotides of the upstream back-spliced
                    # exon are taken and attached at the end of the downstream
                    # back-spliced exon.
                    bsjSeq <- paste0(
                        base::substr(
                            joinedExonSeqs,
                            nchar(joinedExonSeqs) - 10,
                            nchar(joinedExonSeqs)
                        ),
                        base::substr(joinedExonSeqs, 1, 11)
                    )


                    targets[[1]]$length[i] <- nchar(bsjSeq)
                    targets[[1]]$seq[i] <- bsjSeq

                }
            }

        }

        return(targets)

        }
