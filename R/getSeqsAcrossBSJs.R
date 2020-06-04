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
#' @param genome A BSgenome object containing the genome sequences.
#' It can be generated with \code{\link{getBSgenome}}.
#' See \code{\link[BSgenome]{available.genomes}} to see the BSgenome package
#' currently available.
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
#' # Annotate the first back-spliced junctions
#' annotatedBSJs <- annotateBSJs(mergedBSJunctions[1, ], gtf)
#'
#' # Get genome
#' if (requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)){
#' 
#' genome <- BSgenome::getBSgenome("BSgenome.Hsapiens.UCSC.hg19")
#'
#' # Retrieve target sequences
#' targets <- getSeqsAcrossBSJs(
#'     annotatedBSJs,
#'     gtf,
#'     genome)
#' }
#' 
#'
#' @importFrom BSgenome getSeq
#' @import BSgenome.Hsapiens.UCSC.hg19
#' @importFrom Biostrings RNAString
#' @importFrom Biostrings reverseComplement
#' @importFrom rlang .data
#' @import dplyr
#' @export
getSeqsAcrossBSJs <-
    function(annotatedBSJs,
        gtf,
        genome) {
        # Create a enmpty list of 1 elements
        targets <- vector("list", 1)
        names(targets)[1] <- "bsj"

        # Create an empty data frame
        targets[[1]] <- .getTargetsDF(annotatedBSJs)

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

                bsjSeq <- .recreateBSJseq(transcript, genome)
                targets[[1]]$length[i] <- nchar(bsjSeq)
                targets[[1]]$seq[i] <- bsjSeq

            }
        }
        return(targets)
    }


# recreate BSJ sequence
.recreateBSJseq <- function(transcript, genome) {
    if (transcript$strand[1] == "-") {
        transcript <- transcript  %>%
            dplyr::arrange(desc(.data$exon_number))
        # Get exons sequences
        exonSeqs <- .getExonSeqs(transcript, genome)
        # Join sequences
        joinedExonSeqs <- paste(exonSeqs, collapse = "")
        # Recreate bsj seq
        # For negative strand we take the reverse complement
        bsjSeqToReverse <-
            paste0(
                base::substr(
                    joinedExonSeqs,
                    nchar(joinedExonSeqs) - 10,
                    nchar(joinedExonSeqs)
                ),
                base::substr(joinedExonSeqs, 1, 11)
            )
        bsjSeq <-
            as.character(Biostrings::reverseComplement(Biostrings::RNAString(bsjSeqToReverse)))

    } else if (transcript$strand[1] == "+") {
        transcript <- transcript  %>%
            dplyr::arrange(.data$exon_number)
        # Get exons sequences
        exonSeqs <- .getExonSeqs(transcript, genome)
        # Join sequences
        joinedExonSeqs <- paste(exonSeqs, collapse = "")
        # Recreate bsj seq
        # For the positive strand no modification is needed
        bsjSeq <- paste0(
            base::substr(
                joinedExonSeqs,
                nchar(joinedExonSeqs) - 10,
                nchar(joinedExonSeqs)
            ),
            base::substr(joinedExonSeqs, 1, 11)
        )
    }
    return(bsjSeq)
}


# If the function you are looking for is not here check supportFunction.R
# Functions in supportFunction.R are used by multiple functions.
