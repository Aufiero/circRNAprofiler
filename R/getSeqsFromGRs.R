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
#' # Load data frame containing predicted back-spliced junctions
#' data("mergedBSJunctions")
#'
#' # Load short version of the gencode v19 annotation file
#' data("gtf")
#'
#' # Annotate the first 10 back-spliced junctions
#' annotatedBSJs <- annotateBSJs(mergedBSJunctions[1:10, ], gtf)
#'
#' # Retrieve target sequences
#' targets <- getSeqsFromGRs(
#'     annotatedBSJs,
#'     lIntron = 200,
#'     lExon = 10,
#'     type = "ie",
#'     species = "Hsapiens",
#'     genome = "hg19")
#'
#'
#' @importFrom Biostrings reverseComplement
#' @importFrom BSgenome getBSgenome
#' @importFrom BSgenome getSeq
#' @importFrom BiocManager install
#'
#' @export
getSeqsFromGRs <-
    function(annotatedBSJs,
        lIntron = 100,
        lExon = 10,
        type = "ie",
        species = "Hsapiens",
        genome = "hg19") {
        requiredGenome <-
            paste0("BSgenome.", species, ".", "UCSC", ".", genome)

        # only 3 options are possible for the argument type
        match.arg(type, c("ie", "bse", "fi"))

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


        # Retrieve the coordinates from whcih to extract the sequences.
        grCoords <- switch(
            type,
            bse = getBSEcoords(annotatedBSJs),
            fi = getFIcoords(annotatedBSJs),
            ie = getIECoords(annotatedBSJs, lIntron, lExon)
        )

        # Create a list with 2 elements to store the sequences extracted from
        # the upstream and downstream genomic ranges
        targets <- vector("list", 2)
        names(targets)[1] <- "upGR"
        names(targets)[2] <- "downGR"


        for (i in seq_along(targets)) {
            # Create an empty list of 9 data frames
            targets[[i]] <-
                data.frame(matrix(
                    nrow = nrow(annotatedBSJs),
                    ncol = length(getTargetsColNames())
                ))
            colnames(targets[[i]]) <- getTargetsColNames()


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

            # Retrieve the sequences from the defined genomic range coordinates
            genome <- BSgenome::getBSgenome(requiredGenome)

            for (j in seq_along(targets[[1]]$id)) {
                # Here check the coordinates instead of the trancript becaseu
                # even if we have the trancript the back-spliced exons can
                # first or the last of the transcript so in that case we do
                # not have the seq.
                if (targets[[i]]$strand[j] == "-" &
                        !is.na(targets[[i]]$startGR[j]) &
                        !is.na(targets[[i]]$endGR[j])) {
                    # Note: the sequences retrieved by getSeq function from
                    # BSgenome correspond to the positive strand of the DNA
                    # (genome reference). For a circRNA that arises from a gene
                    # trinscribed from negative strand we need to take the
                    # reverse complement of such a sequences. Complement
                    # because the sense strand is the negative strand of the
                    # DNA. Reverse because we always report the sequences
                    # from 5' to 3' in the final data frame.
                    targets[[i]]$seq[j] <-
                        gsub("T",
                            "U",
                            as.character(
                                Biostrings::reverseComplement(
                                    BSgenome::getSeq(
                                        genome,
                                        names = targets[[i]]$chrom[j],
                                        targets[[i]]$startGR[j],
                                        targets[[i]]$endGR[j]
                                    )
                                )
                            ))
                    targets[[i]]$length[j] <-
                        nchar(targets[[i]]$seq[j])


                } else if (targets[[i]]$strand[j] == "+" &
                        !is.na(targets[[i]]$startGR[j]) &
                        !is.na(targets[[i]]$endGR[j])) {
                    # For the positive strand no modification is needed because
                    # the sense strand corresponds to the posiive strand that
                    # is the genome reference.
                    targets[[i]]$seq[j] <-
                        gsub("T", "U", as.character(
                            BSgenome::getSeq(
                                genome,
                                targets[[i]]$chrom[j],
                                targets[[i]]$startGR[j],
                                targets[[i]]$endGR[j]
                            )
                        ))
                    targets[[i]]$length[j] <-
                        nchar(targets[[i]]$seq[j])
                } else{

                }
            }
        }


        return(targets)

        }


#' @title Retrive coordinates of genomic ranges surrounding back-spliced
#' junction coordinates
#'
#' @description The function getIECoords() retrieves the coordinates of the
#' genomic ranges defined by the values of the arguments lIntron and lExon.
#' The start point are the back-spliced junction coordinaes reported in the
#' annotateBSJs data frame.
#'
#' @param annotatedBSJs A data frame with the annotated back-spliced junctions.
#' This data frame can be generated with \code{\link{annotateBSJs}}.
#'
#' @param lIntron An integer indicating how many nucleotides are taken from
#' the introns flanking the back-spliced junctions. This number must be
#' positive.
#'
#' @param lExon An integer indicating how many nucleotides are taken from
#' the back-spliced exons starting from the back-spliced junctions. This number
#' must be positive.
#'
#' @return A data frame.
#'
#' @keywords internal
#'
#' @examples
#' # Load data frame containing predicted back-spliced junctions
#' data("mergedBSJunctions")
#'
#' # Load short version of the gencode v19 annotation file
#' data("gtf")
#'
#' # Example with the first back-spliced junction.
#' # Multiple back-spliced junctions can also be analyzed at the same time.
#'
#' # Annotate predicted back-spliced junctions
#' annotatedBSJs <- annotateBSJs(mergedBSJunctions[1, ], gtf)
#'
#' # Inner function
#' # Retrieve genomic range coordinates
#' getIECoords(annotatedBSJs, lIntron = 101, lExon = 10)
#'
#'
#' @export
getIECoords <- function(annotatedBSJs, lIntron, lExon) {
    # lIntron and lExon must be positive numbers
    if (lIntron < 0) {
        stop("lIntron must be a positive number")
    }

    if (lExon < 0) {
        stop("lExon must be a positive number")
    }

    #Column names of the new data frame
    grColNames <- getGRcolNames()

    # Create an empty dataframe
    grCoords <-
        data.frame(matrix(
            nrow = nrow(annotatedBSJs),
            ncol = length(grColNames)
        ))
    colnames(grCoords) <- grColNames


    for (i in seq_along(annotatedBSJs$id)) {
        grCoords$id[i] <- annotatedBSJs$id[i]
        grCoords$gene[i] <- annotatedBSJs$gene[i]
        grCoords$transcript[i] <- annotatedBSJs$transcript[i]
        grCoords$strand[i] <- annotatedBSJs$strand[i]
        grCoords$chrom[i] <- annotatedBSJs$chrom[i]


        if (!is.na(grCoords$transcript[i])) {
            # POSITIVE STRAND

            # When the back-spliced exons are NOT the FIRST and LAST of
            # the transcript
            if (annotatedBSJs$strand[i] == "+" &
                    !is.na(annotatedBSJs$endUpIntron[i]) &
                    !is.na(annotatedBSJs$startDownIntron[i])) {
                # For the upstream genomic range
                grCoords$startUpGR[i] <-
                    (annotatedBSJs$startUpBSE[i]) - lIntron
                grCoords$endUpGR[i] <-
                    annotatedBSJs$startUpBSE[i] + lExon

                # For the downstream genomic range
                grCoords$startDownGR[i] <-
                    annotatedBSJs$endDownBSE[i] - lExon
                grCoords$endDownGR[i] <-
                    (annotatedBSJs$endDownBSE[i]) + lIntron

                # When the upstream back-spliced exon IS the FIRST of the
                # transcript
            } else if (annotatedBSJs$strand[i] == "+" &
                    is.na(annotatedBSJs$endUpIntron[i]) &
                    !is.na(annotatedBSJs$startDownIntron[i])) {
                # For the upstream genomic range
                grCoords$startUpGR[i] <-
                    (annotatedBSJs$startUpBSE[i])
                grCoords$endUpGR[i] <-
                    annotatedBSJs$startUpBSE[i] + lExon

                # For the downstream genomic range
                grCoords$startDownGR[i] <-
                    annotatedBSJs$endDownBSE[i] - lExon
                grCoords$endDownGR[i] <-
                    (annotatedBSJs$endDownBSE[i]) + lIntron

                # When the downstream back-spliced exon IS the LAST of the
                # transcript
            } else if (annotatedBSJs$strand[i] == "+" &
                    !is.na(annotatedBSJs$endUpIntron[i]) &
                    is.na(annotatedBSJs$startDownIntron[i])) {
                # For the upstream genomic range
                grCoords$startUpGR[i] <-
                    (annotatedBSJs$startUpBSE[i]) - lIntron
                grCoords$endUpGR[i] <-
                    annotatedBSJs$startUpBSE[i] + lExon

                # For the downstream genomic range
                grCoords$startDownGR[i] <-
                    annotatedBSJs$endDownBSE[i] - lExon
                grCoords$endDownGR[i] <-
                    (annotatedBSJs$endDownBSE[i])


                # When the back-spliced exons ARE the FIRST and LAST of the
                # transcript

            } else if (annotatedBSJs$strand[i] == "+" &
                    is.na(annotatedBSJs$endUpIntron[i]) &
                    is.na(annotatedBSJs$startDownIntron[i])) {
                # For the upstream genomic range
                grCoords$startUpGR[i] <-
                    (annotatedBSJs$startUpBSE[i])
                grCoords$endUpGR[i] <-
                    annotatedBSJs$startUpBSE[i] + lExon

                # For the downstream genomic range
                grCoords$startDownGR[i] <-
                    annotatedBSJs$endDownBSE[i] - lExon
                grCoords$endDownGR[i] <-
                    (annotatedBSJs$endDownBSE[i])

                # NEGATIVE STRAND

                # When the back-spliced exons are NOT the FIRST and LAST of the
                # transcript
            } else if (!is.na(annotatedBSJs$strand[i]) &
                    annotatedBSJs$strand[i] == "-" &
                    !is.na(annotatedBSJs$endUpIntron[i]) &
                    !is.na(annotatedBSJs$startDownIntron[i])) {
                # For the upstream genomic range
                grCoords$startUpGR[i] <-
                    annotatedBSJs$startUpBSE[i] - lExon
                grCoords$endUpGR[i] <-
                    annotatedBSJs$startUpBSE[i] + lIntron

                # For the downstream genomic range
                grCoords$startDownGR[i] <-
                    annotatedBSJs$endDownBSE[i] - lIntron
                grCoords$endDownGR[i] <-
                    annotatedBSJs$endDownBSE[i] + lExon

                # When the upstream back-spliced exon IS the FIRST of the
                # transcript
            } else if (!is.na(annotatedBSJs$strand[i]) &
                    annotatedBSJs$strand[i] == "-" &
                    is.na(annotatedBSJs$endUpIntron[i]) &
                    !is.na(annotatedBSJs$startDownIntron[i])) {
                # For the upstream genomic range
                grCoords$startUpGR[i] <-
                    annotatedBSJs$startUpBSE[i] - lExon
                grCoords$endUpGR[i] <-
                    annotatedBSJs$startUpBSE[i]

                # For the downstream genomic range
                grCoords$startDownGR[i] <-
                    annotatedBSJs$endDownBSE[i] - lIntron
                grCoords$endDownGR[i] <-
                    annotatedBSJs$endDownBSE[i] + lExon


                # When the downstream the back-spliced exon IS the LAST of
                # the transcript
            } else if (!is.na(annotatedBSJs$strand[i]) &
                    annotatedBSJs$strand[i] == "-" &
                    !is.na(annotatedBSJs$endUpIntron[i]) &
                    is.na(annotatedBSJs$startDownIntron[i])) {
                # For the upstream genomic range
                grCoords$startUpGR[i] <-
                    annotatedBSJs$startUpBSE[i] - lExon
                grCoords$endUpGR[i] <-
                    annotatedBSJs$startUpBSE[i] + lIntron

                # For the downstream genomic range
                grCoords$startDownGR[i] <-
                    annotatedBSJs$endDownBSE[i]
                grCoords$endDownGR[i] <-
                    annotatedBSJs$endDownBSE[i] + lExon

                # When the back-spliced exons ARE the FIRST and LAST of the
                # transcript
            } else if (!is.na(annotatedBSJs$strand[i]) &
                    annotatedBSJs$strand[i] == "-" &
                    is.na(annotatedBSJs$endUpIntron[i]) &
                    is.na(annotatedBSJs$startDownIntron[i])) {
                # For the upstream genomic range
                grCoords$startUpGR[i] <-
                    annotatedBSJs$startUpBSE[i] - lExon
                grCoords$endUpGR[i] <-
                    annotatedBSJs$startUpBSE[i]


                # For the downstream genomic range
                grCoords$startDownGR[i] <-
                    annotatedBSJs$endDownBSE[i]
                grCoords$endDownGR[i] <-
                    annotatedBSJs$endDownBSE[i] + lExon

            }
        }

    }

    return (grCoords)

}


#' @title Retrieve back-spliced exon coordinates
#'
#' @description The function getBSEcoords() retrieves the coordinates of the
#' back-spliced exons from the annotateBSJs data frame.
#'
#' @param annotatedBSJs A data frame with the annotated back-spliced junctions.
#' This data frame can be generated with \code{\link{annotateBSJs}}.
#'
#' @return A data frame.
#'
#' @keywords internal
#'
#' @examples
#' # Load data frame containing predicted back-spliced junctions
#' data("mergedBSJunctions")
#'
#' # Load short version of the gencode v19 annotation file
#' data("gtf")
#'
#' # Example with the first back-spliced junction
#' # Multiple back-spliced junctions can also be analyzed at the same time.
#'
#' # Annotate predicted back-spliced junctions
#' annotatedBSJs <- annotateBSJs(mergedBSJunctions[1, ], gtf)
#'
#' # Inner function
#' # Retrieve back-spliced exon coordinates
#' getBSEcoords(annotatedBSJs)
#'
#'
#' @export
getBSEcoords <- function(annotatedBSJs) {
    #Column names of the new data frame
    grColNames <- getGRcolNames()

    # Create an empty dataframe
    grCoords <-
        data.frame(matrix(
            nrow = nrow(annotatedBSJs),
            ncol = length(grColNames)
        ))
    colnames(grCoords) <- grColNames


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


#' @title Retrieve flanking intron coordinates
#'
#' @description The function getFIcoords() retrieves the coordinates of
#' the introns flanking the back-spliced exons from the annotateBSJs data frame.
#'
#' @param annotatedBSJs A data frame with the annotated back-spliced junctions.
#' This data frame can be generated with \code{\link{annotateBSJs}}.
#'
#' @return A data frame.
#'
#' @keywords internal
#'
#' @examples
#' # Load data frame containing predicted back-spliced junctions
#' data("mergedBSJunctions")
#'
#' # Load short version of the gencode v19 annotation file
#' data("gtf")
#'
#' # Example with the first back-spliced junction.
#' # Multiple back-spliced junctions can also be analyzed at the same time.
#'
#' # Annotate predicted back-spliced junctions
#' annotatedBSJs <- annotateBSJs(mergedBSJunctions[1, ], gtf)
#'
#' # Inner function
#' # Retrieve flanking intron coordinates
#' getFIcoords(annotatedBSJs)
#'
#'
#' @export
getFIcoords <- function(annotatedBSJs) {
    #Column names of the new data frame
    grColNames <- getGRcolNames()

    # Create an empty dataframe
    grCoords <-
        data.frame(matrix(
            nrow = nrow(annotatedBSJs),
            ncol = length(grColNames)
        ))
    colnames(grCoords) <- grColNames


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


#' The function getGRcolNames() returns the column names
#'
#' @return A character vector.
#'
#' @keywords internal
#'
#' @examples
#' # Inner function
#' getGRcolNames()
#'
#' @export
getGRcolNames <- function() {
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
