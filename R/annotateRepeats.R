#' @title Annotate repetitive elements
#'
#' @description The function annotateRepeats() annotates repetitive elements
#' located in the region flanking the back-spliced junctions of each circRNA.
#' Repetitive elements are provided by AnnotationHub storage which
#' collected repeats from RepeatMasker database. See \code{\link{AnnotationHub}}
#' and \url{http://www.repeatmasker.org} for more details.
#'
#' @param targets A list containing the target regions to analyze.
#' It can be generated with \code{\link{getSeqsFromGRs}}.
#'
#' @param annotationHubID A string specifying the AnnotationHub id to use.
#' Type data(ahRepeatMasker) to see all possibile options. E.g. if AH5122 is
#' specified, repetitve elements from Homo sapiens, genome hg19 will be
#' downloaded and annotated. Defult value is "AH5122".
#'
#' @param complementary A logical specifying whether to filter and report
#' only back-spliced junctions of circRNAs which flanking introns contain
#' complementary repeats, that is, repeats belonging to a same family but
#' located on opposite strands.
#'
#' @return A list.
#'
#' @examples
#' # Load data frame containing detected back-spliced junctions
#' data("mergedBSJunctions")
#'
#' # Load short version of the gencode v19 annotation file
#' data("gtf")
#'
#' # Annotate the first 10 back-spliced junctions
#' annotatedBSJs <- annotateBSJs(mergedBSJunctions[1:10, ], gtf,
#'      isRandom = FALSE)
#'
#' # Retrieve targets
#' targets <- getSeqsFromGRs(
#'     annotatedBSJs,
#'     lIntron = 200,
#'     lExon = 10,
#'     type = "ie",
#'     species = "Hsapiens",
#'     genome = "hg19")
#'
#' # Annotate repeats
#' repeats <-
#' annotateRepeats(targets, annotationHubID  = "AH5122", complementary = TRUE)
#'
#'
#' @import dplyr
#' @import AnnotationHub
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom GenomicRanges findOverlaps
#' @importFrom magrittr %>%
#' @importFrom S4Vectors subjectHits
#' @importFrom S4Vectors queryHits
#' @importFrom rlang .data
#' @importFrom stats setNames
#' @export
annotateRepeats <-
    function(targets,
        annotationHubID = "AH5122",
        complementary = TRUE) {
        if (length(targets) == 2 &
                names(targets)[[1]] == "upGR") {
            # Create an empty list of 2 elements
            repeats <- vector("list", 2)
            names(repeats)[1] <- "upGR"
            names(repeats)[2] <- "downGR"

        }  else {
            stop("target sequences not valid, only upstream and downtream GRs
                are allowed.")
        }

        ah <- AnnotationHub::AnnotationHub()
        rm <- ah[[annotationHubID]]

        # Clean targets dataframe by removing the rows with NA values to avoid
        # getting errors with the makeGRangesFromDataFrame function that does
        # not handle NA values

        index1 <-
            which(
                !is.na(targets[[1]]$transcript) &
                    !is.na(targets[[1]]$startGR) &
                    !is.na(targets[[1]]$endGR)
            )
        index2 <-
            which(
                !is.na(targets[[2]]$transcript) &
                    !is.na(targets[[2]]$startGR) &
                    !is.na(targets[[2]]$endGR)
            )

        targets[[1]] <- targets[[1]][intersect(index1, index2), ]
        targets[[2]] <- targets[[2]][intersect(index1, index2), ]

        for (i in seq_along(repeats)) {
            # Create an empty list of 2 elements to store the extracted
            # information
            repeats[[i]] <- vector("list", 2)
            names(repeats[[i]])[1] <- "targets"
            names(repeats[[i]])[2] <- "repeats"

            # Fill the data frame with the target sequences
            repeats[[i]]$targets <- targets[[i]]


            # Make GR object for the upstream region of the circRNAs
            genRanges <- GenomicRanges::makeGRangesFromDataFrame(
                repeats[[i]]$targets,
                keep.extra.columns = TRUE,
                ignore.strand = FALSE,
                seqinfo = NULL,
                seqnames.field = c("chrom"),
                start.field = c("startGR"),
                end.field = c("endGR"),
                strand.field = "strand",
                starts.in.df.are.0based = FALSE
            )

            # Find Overlaps
            overlaps <-
                suppressWarnings(GenomicRanges::findOverlaps(rm, genRanges, ignore.strand =
                        TRUE))
            if (length(overlaps) == 0) {
                # no genomic ranges in common
                repeats[[i]]$repeats <-
                    data.frame(matrix(nrow = 0, ncol = 8))
                colnames(repeats[[i]]$repeats) <- getRepeatsColNames()
            } else{
                # # Keep only targets where a hit is found
                repeats[[i]]$targets <-
                    repeats[[i]]$targets[S4Vectors::subjectHits(overlaps), ] %>%
                    dplyr::filter(!duplicated(.))

                repeats[[i]]$repeats <-
                    data.frame(genRanges[S4Vectors::subjectHits(overlaps)],
                        rm[S4Vectors::queryHits(overlaps)])


                repeats[[i]]$repeats <- repeats[[i]]$repeats %>%
                    dplyr::select(
                        .data$id,
                        .data$name,
                        .data$seqnames.1,
                        .data$start.1,
                        .data$end.1,
                        .data$width.1,
                        .data$strand.1,
                        .data$score
                    ) %>%
                    dplyr::rename(
                        chrom = .data$seqnames.1,
                        start = .data$start.1,
                        end = .data$end.1,
                        width = .data$width.1,
                        strand = .data$strand.1,
                        score = .data$score
                    )

                repeats[[i]]$repeats <-
                    repeats[[i]]$repeats[!duplicated(repeats[[i]]$repeats), ]
            }
        }

        # Find repeats of the same family located in the upstream and
        # downstream genomic ranges and located on different strands

        if (complementary) {
            upGRs <-
                base::cbind(repeats$upGR$repeats,
                    rep("up", nrow(repeats$upGR$repeats)))
            colnames(upGRs)[9] <- "gr"
            downGRs <-
                base::cbind(repeats$downGR$repeats,
                    rep("down", nrow(repeats$downGR$repeats)))
            colnames(downGRs)[9] <- "gr"


            # Report only instances where the same repeats (same family) is
            # present in the upstream and downstream genomic ranges and are
            # located on different strands.
            overlaps <-
                S4Vectors::findMatches(upGRs$id, downGRs$id, select = "all")
            df <-
                data.frame(upGRs[S4Vectors::queryHits(overlaps), ],
                    downGRs[S4Vectors::subjectHits(overlaps), ])

            matchingRepeats <- df %>%
                dplyr::group_by(.data$id) %>%
                mutate_if(is.factor, as.character) %>%
                dplyr::filter(
                    .data$name == .data$name.1 &
                        .data$gr != .data$gr.1 &
                        .data$strand != .data$strand.1
                )


            # if the repeats are found
            if (nrow(matchingRepeats) > 0) {
                repeats$upGR$repeats <-
                    matchingRepeats[, c(1,2,4,5,6,7,8)] %>%
                    dplyr::arrange(.data$id)

                repeats$upGR$targets <-
                    repeats$upGR$targets %>%
                    dplyr::filter(.data$id %in% unique(matchingRepeats$id)) %>%
                    dplyr::arrange(.data$id)


                repeats$downGR$repeats <-
                    matchingRepeats[, c(10,11,12,13,14,15,16,17)] %>%
                    stats::setNames(gsub(".1", "", names(.))) %>%
                    dplyr::arrange(.data$id)


                repeats$downGR$targets <-
                    repeats$downGR$targets %>%
                    dplyr::filter(.data$id %in% unique(matchingRepeats$id)) %>%
                    dplyr::arrange(.data$id)

            } else {
                repeats[[1]]$repeats <-
                    data.frame(matrix(nrow = 0, ncol = 8))
                colnames(repeats[[1]]$repeats) <- getRepeatsColNames()

                repeats[[2]]$repeats <-
                    data.frame(matrix(nrow = 0, ncol = 8))
                colnames(repeats[[2]]$repeats) <- getRepeatsColNames()
            }
        }
        return(repeats)
    }



#' @title Return column names
#'
#' @description The function getRepeatsColNames() returns the column names.
#'
#' @return A character vector.
#'
#' @keywords internal
#'
#' @examples
#' # Inner function
#' getRepeatsColNames()
#'
#' @export
getRepeatsColNames <- function() {
    repeatsColumns <- c("id",
        "name",
        "chrom",
        "start",
        "end",
        "width",
        "strand",
        "score")
    return(repeatsColumns)

}
