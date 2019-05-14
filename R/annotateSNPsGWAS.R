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
#' @param genome A string specifying the human genome assembly the target
#' regions came from. Possible options are hg38 or hg19.
#' Current image for GWAS SNPS coordinates is hg38. If hg19 is specified
#' SNPs coordinates are realtime liftOver to hg19 coordinates.
#' Internally at first the function \code{\link[gwascat]{makeCurrentGwascat}}
#' from the \code{\link{gwascat}} package is used to get the more recent image.
#' If an error occurs due to connection to the db, data(ebicat37) or
#' data(ebicat38) are used. Default value is "hg19".
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
#' # Retrieve targets
#' targets <- getSeqsFromGRs(
#'     annotatedBSJs,
#'     lIntron = 200,
#'     lExon = 10,
#'     type = "ie",
#'     species = "Hsapiens",
#'     genome = "hg19")
#'
#' # Annotate GWAS SNPs
#' snpsGWAS <- annotateSNPsGWAS(targets, genome = "hg19")
#'
#' @import gwascat
#' @import dplyr
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom GenomicRanges findOverlaps
#' @importFrom magrittr %>%
#' @importFrom utils read.table
#' @importFrom S4Vectors subjectHits
#' @importFrom S4Vectors queryHits
#' @importFrom rlang .data
#' @importFrom utils data
#' @importFrom GenomeInfoDb seqlevelsStyle
#' @export
annotateSNPsGWAS <-
    function(targets, genome = "hg19", pathToTraits = NULL) {
        options(readr.num_columns = 0)
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

        #base::library(gwascat)
        if (genome == "hg19") {

            gwas <- suppressWarnings(tryCatch(
                makeCurrentGwascat(fixNonASCII = FALSE, genome = "GRCh37"),
                error = function(e) NULL
            ))

            if (!is.null(gwas)) {
                GenomeInfoDb::seqlevelsStyle(gwas) <- "UCSC"


            }else{
                #Load lifted over image
                utils::data(ebicat37)
                gwas <- ebicat37
            }

        } else if (genome == "hg38") {
            gwas <- suppressWarnings(tryCatch(
                makeCurrentGwascat(fixNonASCII = FALSE, genome = "GRCh38"),
                error = function(e) NULL
            ))

            if (!is.null(gwas)) {
                GenomeInfoDb::seqlevelsStyle(gwas) <- "UCSC"

            }else{
                # Load image dated 3 August 2015
                utils::data(ebicat38)
                gwas <- ebicat38
            }

        } else {
            stop("Possible genome assembly: hg19, hg38")
        }



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


        # Check if there there are traits
        if (nrow(traitsFromFile) > 0) {
            gwas <- subsetByTraits(gwas, tr = traitsFromFile$id)

        }else{
            cat("traits.txt is empty or absent. All
                   traits in the GWAS catalog will be analyzed")
        }

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

        for (i in seq_along(snpsGWAS)) {
            # Create an empty list of 2 elements to store the extracted
            # information
            snpsGWAS[[i]] <- vector("list", 2)
            names(snpsGWAS[[i]])[1] <- "targets"
            names(snpsGWAS[[i]])[2] <- "snps"

            # Fill the data frame with the target sequences
            snpsGWAS[[i]]$targets <- targets[[i]]

            # Make GR object for the upstream region of the circRNAs
            genRanges <- GenomicRanges::makeGRangesFromDataFrame(
                snpsGWAS[[i]]$targets,
                keep.extra.columns = TRUE,
                ignore.strand = FALSE,
                seqinfo = NULL,
                seqnames.field = c("chrom"),
                start.field = c("startGR"),
                end.field = c("endGR"),
                strand.field = "strand",
                starts.in.df.are.0based = FALSE
            )


            # Find the overlapping gwas snps
            #we can get some wornings if there are no sequence levels in common
            overlaps <-
                suppressWarnings(GenomicRanges::findOverlaps(gwas, genRanges, ignore.strand=TRUE))
            if (length(overlaps) == 0) {
                # no genomic ranges in common
                snpsGWAS[[i]]$snps <-
                    data.frame(matrix(nrow = 0, ncol = 11))
                colnames(snpsGWAS[[i]]$snps) <- c(
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


            } else{

                # Keep only targets where a hit is found
                snpsGWAS[[i]]$targets <-
                    snpsGWAS[[i]]$targets[S4Vectors::subjectHits(overlaps),]%>%
                    dplyr::filter(!duplicated(.))%>%
                    dplyr::arrange(.data$id)

                snpsGWAS[[i]]$snps <-
                    data.frame(genRanges[S4Vectors::subjectHits(overlaps)],
                        gwas[S4Vectors::queryHits(overlaps)])


                snpsGWAS[[i]]$snps <- snpsGWAS[[i]]$snps %>%
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
                    )%>%
                    dplyr::arrange(.data$id)


            }
        }

        return(snpsGWAS)
    }
