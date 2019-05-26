#' @title Annotate GWAS SNPs
#'
#' @description The function annotateSNPsGWAS() annotates GWAS SNPs located
#' in the region flanking the back-spliced junctions of each circRNA.
#' SNPs information including the corresponding genomic coordinates are
#' retrived from the GWAs catalog database.
#' The user can restric the analysis to specific traits/diseases. These must go
#' in the file traits.txt. If this file is absent or empty, all traits in the
#' GWAS catalog are considered in the analysis.
#'
#' @param targets A list containing the target regions to analyze.
#' It can be generated with \code{\link{getSeqsFromGRs}}.
#'
#' @param assembly A string specifying the human genome assembly.
#' Possible options are hg38 or hg19. Current image for GWAS SNPS coordinates
#' is hg38. If hg19 is specified SNPs coordinates are realtime liftOver to
#' hg19 coordinates.
#'
#' @param makeCurrent A logical specifying whether to download the
#' current image of the GWAS catalog.  If FALSE is specified, the function
#' \code{\link[gwascat]{makeCurrentGwascat}} from the \code{\link{gwascat}}
#' package is used to get the more recent image. If TRUE is specified the
#' image data(ebicat37) or data(ebicat38) are used. Default value is FALSE.
#'
#' @param pathToTraits A string containing the path to the traits.txt
#' file. contains diseases/traits specified by the user. It must
#' have one column with header id. By default pathToTraits is set to NULL and
#' the file it is searched in the working directory. If traits.txt is located
#' in a different directory then the path needs to be specified. If this file is
#' absent or empty SNPs associated with all diseases/traits in
#' the GWAS catalog are considered in the SNPs analysis.
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
#' isRandom = FALSE)
#'
#' # Get genome
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
#' # Annotate GWAS SNPs
#' snpsGWAS <- annotateSNPsGWAS(targets, assembly = "hg19")
#'
#' @import dplyr
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom GenomicRanges findOverlaps
#' @importFrom magrittr %>%
#' @importFrom utils read.table
#' @importFrom S4Vectors subjectHits
#' @importFrom S4Vectors queryHits
#' @importFrom rlang .data
#' @importFrom GenomeInfoDb seqlevelsStyle
#' @importFrom gwascat makeCurrentGwascat
#' @importFrom gwascat subsetByTraits
#' @importFrom utils data
#' @export
annotateSNPsGWAS <-
    function(targets,
        assembly = "hg19",
        makeCurrent = FALSE,
        pathToTraits = NULL) {
        options(readr.num_columns = 0)

        # Create a snpsGWAS list
        snpsGWAS <- createSNPsGWASlist(targets)
        # Retrieve SNPs from the GWAS catalog
        gwas <- getGWAS(assembly, makeCurrent)
        # Read traits.txt
        traitsFromFile <- readTraits(pathToTraits)
        # Check if there there are traits
        if (nrow(traitsFromFile) > 0) {
            gwas <- gwascat::subsetByTraits(gwas, tr = traitsFromFile$id)

        } else{
            cat("traits.txt is empty or absent. All
                traits in the GWAS catalog will be analyzed")
        }
        # Clean targets from NA value
        targets <- cleanTargets(targets)

        for (i in seq_along(snpsGWAS)) {
            # Create an empty list of 2 elements to store the extracted
            # information
            snpsGWAS[[i]] <- vector("list", 2)
            names(snpsGWAS[[i]])[1] <- "targets"
            names(snpsGWAS[[i]])[2] <- "snps"

            targetsToAnalyze <- targets[[i]]
            overlaps <- findOverlappingSNPs(gwas, targetsToAnalyze)
            snpsGWAS[[i]]$targets <- overlaps$targets
            snpsGWAS[[i]]$snps <- overlaps$snps

        }
        return(snpsGWAS)
    }



# Get GWAS SNPs
getGWAS <- function(assembly = "hg19",
    makeCurrent = FALSE) {
    if (assembly == "hg19") {
        if (makeCurrent) {
            gwas <- suppressWarnings(tryCatch(
                gwascat::makeCurrentGwascat(fixNonASCII = FALSE, genome = "GRCh37"),
                error = function(e)
                    NULL
            ))
        } else{
            gwas <- NULL
        }

        if (!is.null(gwas)) {
            GenomeInfoDb::seqlevelsStyle(gwas) <- "UCSC"

        } else{
            #Load lifted over image
            utils::data(ebicat37)
            gwas <- ebicat37
        }

    } else if (assembly == "hg38") {
        if (makeCurrent) {
            gwas <- suppressWarnings(tryCatch(
                gwascat::makeCurrentGwascat(fixNonASCII = FALSE, genome = "GRCh38"),
                error = function(e)
                    NULL
            ))
        } else{
            gwas <- NULL
        }

        if (!is.null(gwas)) {
            GenomeInfoDb::seqlevelsStyle(gwas) <- "UCSC"

        } else{
            # Load image dated 3 August 2015
            utils::data(ebicat38)
            gwas <- ebicat38
        }

    } else {
        stop("Possible genome assembly: hg19, hg38")
    }
    return(gwas)
}



# Create a snpsGWAS list
createSNPsGWASlist <- function(targets) {
    if (length(targets) == 2 &
            names(targets)[[1]] == "upGR") {
        # Create an empty list of 2 elements
        snpsGWAS <- vector("list", 2)
        names(snpsGWAS)[1] <- "upGR"
        names(snpsGWAS)[2] <- "downGR"

    }  else {
        stop("target sequences not valid, only upstream and downtream
            GRs are allowed.")
    }

    return(snpsGWAS)
}


# Read traits.txt
readTraits <- function(pathToTraits = NULL) {
    if (is.null(pathToTraits)) {
        pathToTraits <- "traits.txt"
    }

    if (file.exists(pathToTraits)) {
        traitsFromFile <-
            utils::read.table(
                pathToTraits,
                stringsAsFactors = FALSE,
                header = TRUE,
                sep = "\t"
            )
        # colnames(traitsFromFile)[1] <- "id"
    } else{
        traitsFromFile <- data.frame()
    }

    return(traitsFromFile)
}


# get SNPsGWAS column names
getSNPsGWASColNames <- function() {
    colNames <- c(
        "id",
        "snp",
        "chrom",
        "coord",
        "mappedGene",
        "diseaseTrait",
        "pvalue",
        "context",
        "strongestSNPriskAllele",
        "pubmedID",
        "study"
    )

    return(colNames)
}

# Select the needed column and rename snps data frame
renameSNPsGWAS <- function(snps){
    snps <- snps %>%
        dplyr::select(
            .data$id,
            .data$SNPS,
            .data$seqnames.1,
            .data$start.1,
            .data$MAPPED_GENE,
            .data$DISEASE.TRAIT,
            .data$P.VALUE,
            .data$CONTEXT,
            .data$STRONGEST.SNP.RISK.ALLELE,
            .data$PUBMEDID,
            .data$STUDY
        ) %>%
        dplyr::rename(
            snp = .data$SNPS,
            chrom = .data$seqnames.1 ,
            coord = .data$start.1,
            mappedGene = .data$MAPPED_GENE,
            diseaseTrait = .data$DISEASE.TRAIT,
            pvalue = .data$P.VALUE,
            context = .data$CONTEXT,
            strongestSNPriskAllele =
                .data$STRONGEST.SNP.RISK.ALLELE,
            pubmedID = .data$PUBMEDID,
            study = .data$STUDY
        ) %>%
        dplyr::arrange(.data$id)
    return(snps)
}



# Find the overlapping gwas snps
findOverlappingSNPs <- function(gwas, targetsToAnalyze) {

     # Make GR objects
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

    overlappingGRs <-
        suppressWarnings(GenomicRanges::findOverlaps(gwas, genRanges, ignore.strand = TRUE))

    if (length(overlappingGRs) == 0) {
        # No genomic ranges in common
        snps <-  data.frame(matrix(nrow = 0, ncol = 11))
        colnames(snps) <- getSNPsGWASColNames()

        targets <- targetsToAnalyze[NULL, ]

    } else{
        snps <- data.frame(genRanges[S4Vectors::subjectHits(overlappingGRs)],
                gwas[S4Vectors::queryHits(overlappingGRs)])
        snps <- renameSNPsGWAS(snps)

        # Keep only targets where a hit is found
        targets <- targetsToAnalyze[S4Vectors::subjectHits(overlappingGRs), ] %>%
            dplyr::filter(!duplicated(.)) %>%
            dplyr::arrange(.data$id)
    }

    overlaps <- vector("list", 2)
    names(overlaps)[1] <- "snps"
    names(overlaps)[2] <- "targets"

    overlaps$snps <- snps
    overlaps$targets <- targets

    return(overlaps)
}



# If the function you are looking for is not here check supportFunction.R
# Functions in supportFunction.R are used by multiple functions.
