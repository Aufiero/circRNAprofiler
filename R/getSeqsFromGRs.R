#' @title Retrieve sequences flanking back-spliced junctions
#'
#' @description The function getSeqsFromGRs() includes 3 modules to retrieve
#' 3 types of sequences. Sequences of the introns flanking back-spliced
#' junctions, sequences from a defined genomic window surrounding the
#' back-spliced junctions and sequences of the back-spliced exons.
#'
#' @param annotatedBSJs A data frame with the annotated back-spliced junctions.
#' This data frame can be generated with \code{\link{annotateBSJs}}.
#'
#' @param genome A BSgenome object containing the genome sequences.
#' It can be generated with in \code{\link{getBSgenome}}.
#' See \code{\link[BSgenome]{available.genomes}} to see the BSgenome package
#' currently available
#'
#' @param lIntron An integer indicating how many nucleotides are taken from
#' the introns flanking the back-spliced junctions. This number must be
#' positive. Default value is 100.
#'
#' @param lExon An integer indicating how many nucleotides are taken from
#' the back-spliced exons starting from the back-spliced junctions. This number
#' must be positive. Default value is 10.
#'
#' @param type A character string specifying the sequences to retrieve.
#' If type = "ie" the sequences are retrieved from the the genomic ranges
#' defined by using the lIntron and lExon given in input. If type = "bse"
#' the sequences of the back-spliced exons are retrieved. If type = "fi"
#' the sequences of the introns flanking the back-spliced exons are retrieved.
#' Default value is "ie".
#'
#' @return A list.
#'
#' @examples
#' # Load data frame containing predicted back-spliced junctions
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
#' genome <- BSgenome::getBSgenome("BSgenome.Hsapiens.UCSC.hg19")
#'
#' # Retrieve target sequences
#' targets <- getSeqsFromGRs(
#'     annotatedBSJs,
#'     genome,
#'     lIntron = 200,
#'     lExon = 10,
#'     type = "ie"
#'     )
#' }
#' 
#' @importFrom Biostrings reverseComplement
#' @importFrom BSgenome getSeq
#' @import BSgenome.Hsapiens.UCSC.hg19
#' @export
getSeqsFromGRs <-
    function(annotatedBSJs,
        genome,
        lIntron = 100,
        lExon = 10,
        type = "ie") {
        # only 3 options are possible for the argument type
        match.arg(type, c("ie", "bse", "fi"))

        # Retrieve the coordinates from whcih to extract the sequences.
        grCoords <- switch(
            type,
            bse = .getBSEcoords(annotatedBSJs),
            fi = .getFIcoords(annotatedBSJs),
            ie = .getIECoords(annotatedBSJs, lIntron, lExon)
        )
        # Create a list with 2 elements to store the sequences extracted from
        # the upstream and downstream genomic ranges
        targets <- vector("list", 2)
        names(targets)[1] <- "upGR"
        names(targets)[2] <- "downGR"

        for (i in seq_along(targets)) {
            # Create an empty list of data frames
            targets[[i]] <- .getTargetsDF(annotatedBSJs)

            if (i == 1) {
                indexStart <- which(colnames(grCoords) == "startUpGR")
                indexEnd <- which(colnames(grCoords) == "endUpGR")
            } else{
                indexStart <- which(colnames(grCoords) == "startDownGR")
                indexEnd <- which(colnames(grCoords) == "endDownGR")
            }
            # Fill the data frame with the extracted genomic range coordinates
            targets[[i]]$id <- grCoords$id
            targets[[i]]$gene <- grCoords$gene
            targets[[i]]$transcript <- grCoords$transcript
            targets[[i]]$strand <- grCoords$strand
            targets[[i]]$chrom <- grCoords$chrom
            targets[[i]]$startGR <- grCoords[, indexStart]
            targets[[i]]$endGR <- grCoords[, indexEnd]
            targets[[i]]$type <- rep(type, nrow(targets[[i]]))

            for (j in seq_along(targets[[1]]$id)) {
                chrom <- targets[[i]]$chrom[j]
                startGR <- targets[[i]]$startGR[j]
                endGR <- targets[[i]]$endGR[j]
                strand <- targets[[i]]$strand[j]
                seq <-
                    .getSeqFromBSgenome(genome, chrom, startGR, endGR)

                # Check coordinates and retrieve the seq
                if (strand == "-" &
                        !is.na(startGR) & !is.na(endGR)) {
                    # For a circRNA that arises from a gene transcribed from
                    # negative strand we take the reverse complement
                    targets[[i]]$seq[j] <-
                        gsub("T",
                            "U",
                            as.character(Biostrings::reverseComplement(seq)))
                    targets[[i]]$length[j] <-
                        nchar(targets[[i]]$seq[j])

                } else if (strand == "+" &
                        !is.na(startGR) & !is.na(endGR)) {
                    # For the positive strand no modification is needed because
                    # the sense strand corresponds to the positive strand that
                    # is the genome reference.
                    targets[[i]]$seq[j] <-
                        gsub("T", "U", as.character(seq))
                    targets[[i]]$length[j] <-
                        nchar(targets[[i]]$seq[j])
                }
            }
        }
        return(targets)
    }


# The function getIECoords() retrieves the coordinates of the
# genomic ranges defined by the values of the arguments lIntron and lExon.
# The start point are the back-spliced junction coordinates reported in the
# annotateBSJs data frame.
.getIECoords <- function(annotatedBSJs,
    lIntron = 100,
    lExon = 10) {
    # lIntron and lExon must be positive numbers
    if (lIntron < 0) {
        stop("lIntron must be a positive number")
    }
    if (lExon < 0) {
        stop("lExon must be a positive number")
    }
    # Create an empty dataframe
    grCoords <- .getGRcoordsDF(annotatedBSJs)

    for (i in seq_along(annotatedBSJs$id)) {
        grCoords$id[i] <- annotatedBSJs$id[i]
        grCoords$gene[i] <- annotatedBSJs$gene[i]
        grCoords$transcript[i] <- annotatedBSJs$transcript[i]
        grCoords$strand[i] <- annotatedBSJs$strand[i]
        grCoords$chrom[i] <- annotatedBSJs$chrom[i]

        if (!is.na(annotatedBSJs$transcript[i]) &
                !is.na(annotatedBSJs$strand[i])) {
            annotatedBSJsToAnalyze <- annotatedBSJs[i,]
            # Get upstream and dowstream genomic ranges coordinates
            # for POSITIVE STRAND
            if (annotatedBSJsToAnalyze$strand == "+") {
                grCoordsToAnalyze <-
                    .grCoordsForPositive(annotatedBSJsToAnalyze, lIntron, lExon)
                grCoords$startUpGR[i] <- grCoordsToAnalyze$startUpGR
                grCoords$endUpGR[i] <- grCoordsToAnalyze$endUpGR
                grCoords$startDownGR[i] <-
                    grCoordsToAnalyze$startDownGR
                grCoords$endDownGR[i] <- grCoordsToAnalyze$endDownGR

            }
            # Get upstream and dowstream genomic ranges coordinates
            # for NEGATIVE STRAND
            if (annotatedBSJsToAnalyze$strand == "-") {
                grCoordsToAnalyze <-
                    .grCoordsForNegative(annotatedBSJsToAnalyze, lIntron, lExon)
                grCoords$startUpGR[i] <- grCoordsToAnalyze$startUpGR
                grCoords$endUpGR[i] <- grCoordsToAnalyze$endUpGR
                grCoords$startDownGR[i] <-
                    grCoordsToAnalyze$startDownGR
                grCoords$endDownGR[i] <- grCoordsToAnalyze$endDownGR
            }
        }
    }
    return (grCoords)
}

# Get upstream and dowstream genomic ranges coordinates for positive strand
.grCoordsForPositive <-
    function(annotatedBSJsToAnalyze,
        lIntron = 100,
        lExon = 10) {
        # Create an empty dataframe
        grCoords <- .getGRcoordsDF(annotatedBSJsToAnalyze)
        # When the BSEs are NOT the FIRST and LAST of the transcript
        if (!is.na(annotatedBSJsToAnalyze$endUpIntron) &
                !is.na(annotatedBSJsToAnalyze$startDownIntron)) {
            # For the upstream genomic range
            grCoords$startUpGR <-
                (annotatedBSJsToAnalyze$startUpBSE) - lIntron
            grCoords$endUpGR <-
                annotatedBSJsToAnalyze$startUpBSE + lExon

            # For the downstream genomic range
            grCoords$startDownGR <-
                annotatedBSJsToAnalyze$endDownBSE - lExon
            grCoords$endDownGR <-
                (annotatedBSJsToAnalyze$endDownBSE) + lIntron

            # When the upstream back-spliced exon IS the FIRST of the
            # transcript
        } else if (is.na(annotatedBSJsToAnalyze$endUpIntron) &
                !is.na(annotatedBSJsToAnalyze$startDownIntron)) {
            # For the upstream genomic range
            grCoords$startUpGR <-
                (annotatedBSJsToAnalyze$startUpBSE)
            grCoords$endUpGR <-
                annotatedBSJsToAnalyze$startUpBSE + lExon

            # For the downstream genomic range
            grCoords$startDownGR <-
                annotatedBSJsToAnalyze$endDownBSE - lExon
            grCoords$endDownGR <-
                (annotatedBSJsToAnalyze$endDownBSE) + lIntron

            # When the downstream back-spliced exon IS the LAST of the
            # transcript
        } else if (!is.na(annotatedBSJsToAnalyze$endUpIntron) &
                is.na(annotatedBSJsToAnalyze$startDownIntron)) {
            # For the upstream genomic range
            grCoords$startUpGR <-
                (annotatedBSJsToAnalyze$startUpBSE) - lIntron
            grCoords$endUpGR <-
                annotatedBSJsToAnalyze$startUpBSE + lExon

            # For the downstream genomic range
            grCoords$startDownGR <-
                annotatedBSJsToAnalyze$endDownBSE - lExon
            grCoords$endDownGR <-
                (annotatedBSJsToAnalyze$endDownBSE)

            # When the back-spliced exons ARE the FIRST and LAST of the
            # transcript
        } else if (is.na(annotatedBSJsToAnalyze$endUpIntron) &
                is.na(annotatedBSJsToAnalyze$startDownIntron)) {
            # For the upstream genomic range
            grCoords$startUpGR <-
                (annotatedBSJsToAnalyze$startUpBSE)
            grCoords$endUpGR <-
                annotatedBSJsToAnalyze$startUpBSE + lExon

            # For the downstream genomic range
            grCoords$startDownGR <-
                annotatedBSJsToAnalyze$endDownBSE - lExon
            grCoords$endDownGR <-
                (annotatedBSJsToAnalyze$endDownBSE)
        }

        grCoords <- grCoords %>%
            dplyr::select(startUpGR,
                endUpGR,
                startDownGR,
                endDownGR)
        return(grCoords)
    }


# Get startUpGR and endUpGR for negative strand
.grCoordsForNegative <-
    function(annotatedBSJsToAnalyze,
        lIntron = 100,
        lExon = 10) {
        # Create an empty dataframe
        grCoords <- .getGRcoordsDF(annotatedBSJsToAnalyze)
        # When the back-spliced exons are NOT the FIRST and LAST of the
        # transcript
        if (!is.na(annotatedBSJsToAnalyze$endUpIntron) &
                !is.na(annotatedBSJsToAnalyze$startDownIntron)) {
            # For the upstream genomic range
            grCoords$startUpGR <-
                annotatedBSJsToAnalyze$startUpBSE - lExon
            grCoords$endUpGR <-
                annotatedBSJsToAnalyze$startUpBSE + lIntron

            # For the downstream genomic range
            grCoords$startDownGR <-
                annotatedBSJsToAnalyze$endDownBSE - lIntron
            grCoords$endDownGR <-
                annotatedBSJsToAnalyze$endDownBSE + lExon

            # When the upstream back-spliced exon IS the FIRST of the
            # transcript
        } else if (is.na(annotatedBSJsToAnalyze$endUpIntron) &
                !is.na(annotatedBSJsToAnalyze$startDownIntron)) {
            # For the upstream genomic range
            grCoords$startUpGR <-
                annotatedBSJsToAnalyze$startUpBSE - lExon
            grCoords$endUpGR <- annotatedBSJsToAnalyze$startUpBSE

            # For the downstream genomic range
            grCoords$startDownGR <-
                annotatedBSJsToAnalyze$endDownBSE - lIntron
            grCoords$endDownGR <-
                annotatedBSJsToAnalyze$endDownBSE + lExon

            # When the downstream the back-spliced exon IS the LAST of
            # the transcript
        } else if (!is.na(annotatedBSJsToAnalyze$endUpIntron) &
                is.na(annotatedBSJsToAnalyze$startDownIntron)) {
            # For the upstream genomic range
            grCoords$startUpGR <-
                annotatedBSJsToAnalyze$startUpBSE - lExon
            grCoords$endUpGR <-
                annotatedBSJsToAnalyze$startUpBSE + lIntron

            # For the downstream genomic range
            grCoords$startDownGR <-
                annotatedBSJsToAnalyze$endDownBSE
            grCoords$endDownGR <-
                annotatedBSJsToAnalyze$endDownBSE + lExon

            # When the back-spliced exons ARE the FIRST and LAST of the
            # transcript
        } else if (is.na(annotatedBSJsToAnalyze$endUpIntron) &
                is.na(annotatedBSJsToAnalyze$startDownIntron)) {
            # For the upstream genomic range
            grCoords$startUpGR <-
                annotatedBSJsToAnalyze$startUpBSE - lExon
            grCoords$endUpGR <- annotatedBSJsToAnalyze$startUpBSE

            # For the downstream genomic range
            grCoords$startDownGR <-
                annotatedBSJsToAnalyze$endDownBSE
            grCoords$endDownGR <-
                annotatedBSJsToAnalyze$endDownBSE + lExon
        }

        grCoords <- grCoords %>%
            dplyr::select(startUpGR,
                endUpGR,
                startDownGR,
                endDownGR)

        return(grCoords)
    }



# The function getBSEcoords() retrieves the coordinates of the bse
.getBSEcoords <- function(annotatedBSJs) {
    # Create an empty dataframe
    grCoords <- .getGRcoordsDF(annotatedBSJs)

    for (i in seq_along(annotatedBSJs$id)) {
        grCoords$id[i] <- annotatedBSJs$id[i]
        grCoords$gene[i] <- annotatedBSJs$gene[i]
        grCoords$transcript[i] <- annotatedBSJs$transcript[i]
        grCoords$strand[i] <- annotatedBSJs$strand[i]
        grCoords$chrom[i] <- annotatedBSJs$chrom[i]

        if (!is.na(grCoords$transcript[i])) {
            # POSITIVE STRAND
            if (annotatedBSJs$strand[i] == "+") {
                # For the upstream genomic range
                grCoords$startUpGR[i] <- annotatedBSJs$startUpBSE[i]
                grCoords$endUpGR[i] <- annotatedBSJs$endUpBSE[i]
                # For the downstream genomic range
                grCoords$startDownGR[i] <-
                    annotatedBSJs$startDownBSE[i]
                grCoords$endDownGR[i] <- annotatedBSJs$endDownBSE[i]

                # NEGATIVE STRAND
            } else if (annotatedBSJs$strand[i] == "-") {
                # For the upstream genomic range
                grCoords$startUpGR[i] <- annotatedBSJs$endUpBSE[i]
                grCoords$endUpGR[i] <- annotatedBSJs$startUpBSE[i]
                # For the downstream genomic range
                grCoords$startDownGR[i] <-
                    annotatedBSJs$endDownBSE[i]
                grCoords$endDownGR[i] <-
                    annotatedBSJs$startDownBSE[i]

            }
        }
    }
    return (grCoords)

}


#The function getFIcoords() retrieves the coordinates of
# the introns flanking the back-spliced exons from the annotateBSJs data frame.
.getFIcoords <- function(annotatedBSJs) {
    # Create an empty dataframe
    grCoords <- .getGRcoordsDF(annotatedBSJs)

    for (i in seq_along(annotatedBSJs$id)) {
        grCoords$id[i] <- annotatedBSJs$id[i]
        grCoords$gene[i] <- annotatedBSJs$gene[i]
        grCoords$transcript[i] <- annotatedBSJs$transcript[i]
        grCoords$strand[i] <- annotatedBSJs$strand[i]
        grCoords$chrom[i] <- annotatedBSJs$chrom[i]

        if (!is.na(grCoords$transcript[i])) {
            # POSITIVE STRAND
            if (annotatedBSJs$strand[i] == "+") {
                # For the upstream genomic range
                grCoords$startUpGR[i] <-
                    annotatedBSJs$startUpIntron[i]
                grCoords$endUpGR[i] <- annotatedBSJs$endUpIntron[i]

                # For the downstream genomic range
                grCoords$startDownGR[i] <-
                    annotatedBSJs$startDownIntron[i]
                grCoords$endDownGR[i] <-
                    annotatedBSJs$endDownIntron[i]

                # For NEGATIVE strand
            } else if (annotatedBSJs$strand[i] == "-") {
                # For the upstream genomic range
                grCoords$startUpGR[i] <-
                    annotatedBSJs$endUpIntron[i]
                grCoords$endUpGR[i] <-
                    annotatedBSJs$startUpIntron[i]

                # For the downstream genomic range
                grCoords$startDownGR[i] <-
                    annotatedBSJs$endDownIntron[i]
                grCoords$endDownGR[i] <-
                    annotatedBSJs$startDownIntron[i]
            }
        }
    }
    return (grCoords)
}


# The function getGRcolNames() returns the column names
.getGRcolNames <- function() {
    grColNames <- c(
        "id",
        "gene",
        "transcript",
        "strand",
        "chrom",
        "startUpGR",
        "endUpGR",
        "startDownGR",
        "endDownGR"
    )
    return(grColNames)
}

# Create grCoords data frame
.getGRcoordsDF <- function(annotatedBSJs) {
    grCoords <-
        data.frame(matrix(
            nrow = nrow(annotatedBSJs),
            ncol = length(.getGRcolNames())
        ))
    colnames(grCoords) <- .getGRcolNames()
    return(grCoords)
}


# get sequences from BS genome
.getSeqFromBSgenome <- function(genome, chrom, startGR, endGR) {
    seq <- BSgenome::getSeq(genome,
        names = chrom,
        startGR,
        endGR)
    return(seq)
}


# If the function you are looking for is not here check supportFunction.R
# Functions in supportFunction.R are used by multiple functions.
