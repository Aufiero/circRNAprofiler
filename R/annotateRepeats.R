#' @title Annotate repetitive elements
#'
#' @description The function annotateRepeats() annotates repetitive elements
#' located in the region flanking the back-spliced junctions of each circRNA.
#' Repetitive elements are provided by AnnotationHub storage which
#' collected repeats from RepeatMasker database. See \code{\link{AnnotationHub}}
#' and \url{http://www.repeatmasker.org} for more details.
#' An empty list is returned if none overlapping repeats  are found.
#'
#' @param targets A list containing the target regions to analyze.
#' It can be generated with \code{\link{getSeqsFromGRs}}.
#'
#' @param annotationHubID A string specifying the AnnotationHub id to use.
#' Type data(ahRepeatMasker) to see all possible options. E.g. if AH5122 is
#' specified, repetitive elements from Homo sapiens, genome hg19 will be
#' downloaded and annotated. Default value is "AH5122".
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
#' # Annotate the first back-spliced junctions
#' annotatedBSJs <- annotateBSJs(mergedBSJunctions[1, ], gtf)
#'
#' # Get genome
#' if (requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)){
#'
#' genome <- BSgenome::getBSgenome("BSgenome.Hsapiens.UCSC.hg19")
#'
#' # Retrieve targets
#' targets <- getSeqsFromGRs(
#'     annotatedBSJs,
#'     genome,
#'     lIntron = 200,
#'     lExon = 10,
#'     type = "ie"
#'     )
#'
#' # Annotate repeats
#'
#' # repeats <- annotateRepeats(targets, annotationHubID  = "AH5122", complementary = TRUE)
#'
#' }
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
        
        # Clean targets from NA value
        targets <- .cleanTargets(targets)
        
        for (i in seq_along(repeats)) {
            # Create an empty list of 2 elements to store the extracted
            # information
            repeats[[i]] <- vector("list", 2)
            names(repeats[[i]])[1] <- "targets"
            names(repeats[[i]])[2] <- "repeats"
            
            targetsToAnalyze <- targets[[i]]
            overlaps <-
                .findOverlappingRepeats(rm, targetsToAnalyze)
            repeats[[i]]$targets <- overlaps$targets
            repeats[[i]]$repeats <- overlaps$repeats
            
        }
        
        # Find repeats of the same family located in the upstream and
        # downstream genomic ranges and located on different strands
        
        if (complementary) {
            repeats <- .getComplRepeats(repeats)
        }
        
        return(repeats)
    }



# The function getRepeatsColNames() returns the column names.
.getRepeatsColNames <- function() {
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

# The function getRepeatsColNames() returns complementary repeats.
.getComplRepeats <- function(repeats) {
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
    
    if (length(overlaps) > 0) {
        df <-
            data.frame(upGRs[S4Vectors::queryHits(overlaps),],
                       downGRs[S4Vectors::subjectHits(overlaps),])
        
        
        matchingRepeats <- df %>%
            dplyr::mutate_if(is.factor, as.character) %>%
            dplyr::filter(
                name == name.1 &
                    gr != gr.1 &
                    strand != strand.1
            )
        
        repeats$upGR$repeats <-
            matchingRepeats[, c(1, 2, 4, 5, 6, 7, 8)] %>%
            dplyr::arrange(id)
        
        repeats$upGR$targets <-
            repeats$upGR$targets %>%
            dplyr::filter(id %in% unique(matchingRepeats$id)) %>%
            dplyr::arrange(id)
        
        
        repeats$downGR$repeats <-
            matchingRepeats[, c(10, 11, 12, 13, 14, 15, 16, 17)] %>%
            stats::setNames(gsub(".1", "", names(.))) %>%
            dplyr::arrange(id)
        
        
        repeats$downGR$targets <-
            repeats$downGR$targets %>%
            dplyr::filter(id %in% unique(matchingRepeats$id)) %>%
            dplyr::arrange(id)
        
    } else {
        repeats[[1]]$repeats <-
            data.frame(matrix(nrow = 0, ncol = 8))
        colnames(repeats[[1]]$repeats) <-
            .getRepeatsColNames()
        
        repeats[[2]]$repeats <-
            data.frame(matrix(nrow = 0, ncol = 8))
        colnames(repeats[[2]]$repeats) <-
            .getRepeatsColNames()
    }
    
    return(repeats)
}

# Select the needed column and rename repeats data frame
.renameRepeats <- function(repeats) {
    repeats <- repeats %>%
        dplyr::select(
            id,
            name,
            seqnames.1,
            start.1,
            end.1,
            width.1,
            strand.1,
            score
        ) %>%
        dplyr::rename(
            chrom = seqnames.1,
            start = start.1,
            end = end.1,
            width = width.1,
            strand = strand.1,
            score = score
        )
    
    repeats <- repeats[!duplicated(repeats),]
    return(repeats)
}


.findOverlappingRepeats <- function(rm, targetsToAnalyze) {
    # Make GR object for the upstream region of the circRNAs
    genRanges <- GenomicRanges::makeGRangesFromDataFrame(
        targetsToAnalyze,
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
        # No genomic ranges in common
        repeats <- data.frame(matrix(nrow = 0, ncol = 8))
        colnames(repeats) <- .getRepeatsColNames()
        
        targets <- targetsToAnalyze[NULL, ]
        
    } else{
        repeats <- data.frame(genRanges[S4Vectors::subjectHits(overlaps)],
                              rm[S4Vectors::queryHits(overlaps)])
        repeats <- .renameRepeats(repeats)
        
        # Keep only targets where a hit is found
        targets <-
            repeats[S4Vectors::subjectHits(overlaps), ] %>%
            dplyr::filter(!duplicated(.))
    }
    
    overlaps <- vector("list", 2)
    names(overlaps)[1] <- "repeats"
    names(overlaps)[2] <- "targets"
    
    overlaps$repeats <- repeats
    overlaps$targets <- targets
    return(overlaps)
}


# If the function you are looking for is not here check supportFunction.R
# Functions in supportFunction.R are used by multiple functions.
