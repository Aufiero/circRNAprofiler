#' @title LiftOver back-spliced junction coordinates
#'
#' @description The function liftBSJCoords() converts back-spliced junction
#' coordinates between species ad genome assemblies by using the liftOver utility
#' from UCSC.
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
#' # LifOver the first 10 back-spliced junction coordinates
#' liftedBSJCoords <- liftBSJCoords(mergedBSJunctions[1:10,], map = "hg19ToMm9")
#'
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
liftBSJCoords <-
    function(backSplicedJunctions,
        map = "hg19ToMm9",
        annotationHubID = "AH14155") {
        species1 <- base::strsplit(map, "To")[[1]][1]
        species2 <- base::strsplit(map, "To")[[1]][2]

        # Get basic colum names
        basicColumns <- getBasicColNames()

        # Create an empty data frame
        liftedBSJCoords <-
            data.frame(matrix(
                nrow = nrow(backSplicedJunctions),
                ncol = length(basicColumns) + 1
            ))
        colnames(liftedBSJCoords) <-
            c(basicColumns, paste(basicColumns[1], species1, sep = "_"))


        ah <- AnnotationHub::AnnotationHub()
        options(warn = -1) # warning is ignored
        # Import chain. file
        chain <- ah[[annotationHubID]]

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

        # Lift Over coordinates
        liftedOverStartUpBSE <-
            data.frame(rtracklayer::liftOver(grStartUpBSE, chain))
        liftedOverEndDownBSE <-
            data.frame(rtracklayer::liftOver(grEndDownBSE, chain))

        options(warn = 0) # default
        mtStart <-
            match(rownames(backSplicedJunctions),
                liftedOverStartUpBSE$group)
        mtEnd <-
            match(rownames(backSplicedJunctions),
                liftedOverEndDownBSE$group)

        # if the species to convert to the coordinates is human the gene names
        # need to be all upper case
        if (grepl("Hg", species2)) {
            liftedBSJCoords$gene <- toupper(backSplicedJunctions$gene)
        } else{
            liftedBSJCoords$gene <- backSplicedJunctions$gene
        }


        liftedBSJCoords$strand <-
            unlist(as.character(liftedOverStartUpBSE$strand[mtStart]))
        liftedBSJCoords$chrom <-
            unlist(as.character(liftedOverStartUpBSE$seqnames[mtStart]))
        liftedBSJCoords$startUpBSE <-
            liftedOverStartUpBSE$start[mtStart]
        liftedBSJCoords$endDownBSE <-
            liftedOverEndDownBSE$end[mtEnd]


        liftedBSJCoords[, paste(basicColumns[1], species1, sep = "_")] <-
            backSplicedJunctions$id


        # It can happen that some genes are located on different strands in
        # different species so a check needs to be done to switch
        # (if necessary) the coordinates.
        for (i in seq_along(liftedBSJCoords$id)) {
            if (!is.na(liftedBSJCoords$strand[i])  &
                    !is.na(liftedBSJCoords$chrom[i]) &
                    !is.na(liftedBSJCoords$startUpBSE[i]) &
                    !is.na(liftedBSJCoords$endDownBSE[i]) &
                    liftedBSJCoords$strand[i] == "+" &
                    liftedBSJCoords$startUpBSE[i] >
                    liftedBSJCoords$endDownBSE[i]) {
                x <- liftedBSJCoords$endDownBSE[i]
                liftedBSJCoords$startUpBSE[i] <-
                    liftedBSJCoords$endDownBSE[i]
                liftedBSJCoords$endDownBSE[i] <- x

            } else if (!is.na(liftedBSJCoords$strand[i]) &
                    !is.na(liftedBSJCoords$chrom[i]) &
                    !is.na(liftedBSJCoords$startUpBSE[i]) &
                    !is.na(liftedBSJCoords$endDownBSE[i]) &
                    liftedBSJCoords$strand[i] == "-" &
                    liftedBSJCoords$startUpBSE[i] <
                    liftedBSJCoords$endDownBSE[i]) {
                x <- liftedBSJCoords$endDownBSE[i]
                liftedBSJCoords$startUpBSE[i] <-
                    liftedBSJCoords$endDownBSE[i]
                liftedBSJCoords$endDownBSE[i] <- x
            }
        }

        # Generate a unique identifier by combining the values of the following
        # columns: gene, strand, chrom, end_downBSExon and start_upBSExon.
        # The values are separated by a semicolumns (:)
        id <- paste(
            liftedBSJCoords$gene,
            liftedBSJCoords$strand,
            liftedBSJCoords$chrom,
            liftedBSJCoords$startUpBSE,
            liftedBSJCoords$endDownBSE,
            sep = ":"
        )

        liftedBSJCoords$id <- id


        return(liftedBSJCoords)

    }
