#' @title LiftOver back-spliced junction coordinates
#'
#' @description The function liftBSJcoords() maps back-spliced junction
#' coordinates between species ad genome assemblies by using the liftOver utility
#' from UCSC. Only back-spliced junction coordinates where the mapping was
#' successful are reported.
#'
#' @param backSplicedJunctions A data frame containing the back-spliced junction
#' coordinates (predicted or randomly selected).
#' See \code{\link{getRandomBSJunctions}}, \code{\link{getBackSplicedJunctions}}
#' and \code{\link{mergeBSJunctions}} (to group circRNAs detected by multiple
#' detection tools), on how to generated this data frame.
#'
#' @param map A character string in the format <db1>To<Db2> (e.g."hg19ToMm9")
#' specifying the reference genome mapping logic associated with a
#' valid .chain file. Default value is "hg19ToMm9".
#'
#' @param annotationHubID A string specifying the AnnotationHub id associated
#' with a valid *.chain file. Type data(ahChainFiles) to see all possibile
#' options. E.g. if AH14155 is specified, the hg19ToMm9.over.chain.gz will
#' be used to convert the hg19 (Human GRCh37) coordinates to mm10
#' (Mouse GRCm38). Default value is "AH14155".
#'
#' @return a data frame.
#'
#' @examples
#' # Load a data frame containing detected back-spliced junctions
#' data("mergedBSJunctions")
#' 
#' # LiftOver the first 10 back-spliced junction coordinates
#' liftedBSJcoords <- liftBSJcoords(mergedBSJunctions[1:10,], map = "hg19ToMm9")
#'
#' @importFrom IRanges IRanges
#' @importFrom rtracklayer import.chain
#' @importFrom rtracklayer liftOver
#' @importFrom GenomicRanges GRanges
#' @importFrom R.utils gunzip
#' @importFrom utils write.table
#' @import AnnotationHub
#'
#' @export
liftBSJcoords <-
    function(backSplicedJunctions,
        map = "hg19ToMm9",
        annotationHubID = "AH14155") {
        species1 <- base::strsplit(map, "To")[[1]][1]
        species2 <- base::strsplit(map, "To")[[1]][2]

        # Create an empty data frame
        liftedBSJCoords <- .createLiftedBSJCoordsDF(backSplicedJunctions, species1)
        # Capitalize gene name if species is human
        liftedBSJCoords$gene <- .upperHuman(species2, backSplicedJunctions)

        options(warn = -1)    # warning is ignored
        chain <- .getChain(annotationHubID)
        # Create GR objects
        genRanges <- .createGRsForLift(backSplicedJunctions)
        grStartUpBSE <- genRanges$grStartUpBSE
        grEndDownBSE <- genRanges$grEndDownBSE
        # Lift Over coordinates
        liftedOverStartUpBSE <-
            data.frame(rtracklayer::liftOver(grStartUpBSE, chain))
        liftedOverEndDownBSE <-
            data.frame(rtracklayer::liftOver(grEndDownBSE, chain))

        options(warn = 0) # default
        mtStart <- match(rownames(backSplicedJunctions), liftedOverStartUpBSE$group)
        mtEnd <- match(rownames(backSplicedJunctions), liftedOverEndDownBSE$group)
        liftedBSJCoords$strand <- unlist(as.character(liftedOverStartUpBSE$strand[mtStart]))
        liftedBSJCoords$chrom <- unlist(as.character(liftedOverStartUpBSE$seqnames[mtStart]))
        liftedBSJCoords$startUpBSE <-liftedOverStartUpBSE$start[mtStart]
        liftedBSJCoords$endDownBSE <- liftedOverEndDownBSE$end[mtEnd]

        # Get basic colum names
        basicColumns <- .getBasicColNames()
        liftedBSJCoords[, paste(basicColumns[1], species1, sep = "_")] <-
            backSplicedJunctions$id

        # Remove non converted bsjs
        liftedBSJCoords <- liftedBSJCoords[!is.na(liftedBSJCoords$strand),]
        # A gene can be located on different strands in different species
        # Fix (if necessary) the coordinates.
        liftedBSJCoords <- .fixCoords(liftedBSJCoords)

        # Generate a unique identifier
        id <- .getID(liftedBSJCoords)
        liftedBSJCoords$id <- id
        return(liftedBSJCoords)
    }


# Create liftedBSJCoords data frame
.createLiftedBSJCoordsDF <- function(backSplicedJunctions, species1){

    # Get basic colum names
    basicColumns <- .getBasicColNames()

    # Create an empty data frame
    liftedBSJCoords <-
        data.frame(matrix(
            nrow = nrow(backSplicedJunctions),
            ncol = length(basicColumns) + 1
        ))
    colnames(liftedBSJCoords) <-
        c(basicColumns, paste(basicColumns[1], species1, sep = "_"))

    return(liftedBSJCoords)

}


# Get chain file from AnnotationHub
.getChain <- function(annotationHubID){
    ah <- AnnotationHub::AnnotationHub()

    # Import chain. file
    chain <- ah[[annotationHubID]]
    return(chain)
}

# Create GR objects
.createGRsForLift <- function(backSplicedJunctions){
    # Create GR objects
    grStartUpBSE <- GenomicRanges::GRanges(
        seqnames = backSplicedJunctions$chrom,
        ranges = IRanges::IRanges(start = backSplicedJunctions$startUpBSE,
            end = backSplicedJunctions$startUpBSE),
        strand = backSplicedJunctions$strand
    )
    grEndDownBSE <- GenomicRanges::GRanges(
        seqnames = backSplicedJunctions$chrom,
        ranges = IRanges::IRanges(start = backSplicedJunctions$endDownBSE,
            end = backSplicedJunctions$endDownBSE),
        strand = backSplicedJunctions$strand
    )

    genRanges <- vector("list", 2)
    names(genRanges)[1] <- "grStartUpBSE"
    names(genRanges)[2] <- "grEndDownBSE"

    genRanges$grStartUpBSE <- grStartUpBSE
    genRanges$grEndDownBSE <- grEndDownBSE

    return(genRanges)
}


# Capitalize gene name if species is human
.upperHuman <- function(species, backSplicedJunctions){
    if (grepl("Hg", species)) {
        gene <- toupper(backSplicedJunctions$gene)
    } else{
        x <- tolower(backSplicedJunctions$gene)
        substr(x, 1, 1) <- toupper(substr(x, 1, 1))
        gene <- x
    }
    return(gene)
}


# If the function you are looking for is not here check supportFunction.R
# Functions in supportFunction.R are used by multiple functions.
