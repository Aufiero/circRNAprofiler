# Satisfy R CMD check

# data ebicat37 and ebicat38 from gwascat package are deprecated.
# if (getRversion() >= "3.1.0")
#     utils::globalVariables(c(
#         "gwascat", ".", "ebicat37","ebicat38")) 

if (getRversion() >= "3.1.0")
    utils::globalVariables(".")


#' @title Format annotation file
#'
#' @description The function formatGTF() formats the given annotation file.
#'
#' @param pathToGTF A string containing the path to the GTF file.
#' Use the same annotation file used during the RNA-seq mapping procedure.
#' If .gtf file is not present in the current working directory the full path
#' should be specified.
#'
#' @return A data frame.
#'
#' @examples
#' gtf <- formatGTF()
#'
#' @importFrom stringr str_extract
#' @importFrom rtracklayer import
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @import dplyr
#'
#' @export
formatGTF <- function(pathToGTF = NULL) {
    if (!is.null(pathToGTF)) {
        # Read GTF file
        gtf <- .readGTF(pathToGTF)
        # Get needed column in the GTF file
        needColumns <- .getNeededColumn()

        # For UCSC and Ncbi the exon_number column needs to be added
        if (!("exon_number" %in% colnames(gtf))) {
            pos <-  gtf  %>%
                dplyr::filter(strand == "+") %>%
                dplyr::group_by(transcript_id) %>%
                dplyr::arrange(start, .by_group = TRUE) %>%
                dplyr::mutate (exon_number = row_number())

            neg <-  gtf  %>%
                dplyr::filter(strand == "-") %>%
                dplyr::group_by(transcript_id) %>%
                dplyr::arrange(start, .by_group = TRUE) %>%
                dplyr::mutate (exon_number = rev(row_number()))

            formattedGTF <- rbind(pos, neg)
            formattedGTF <- formattedGTF %>%
                dplyr::select(all_of(needColumns), "exon_number") %>%
                dplyr::mutate(exon_number = as.numeric(exon_number)) %>%
                as.data.frame()

        } else {
            formattedGTF <- gtf %>%
                dplyr::select(all_of(needColumns), "exon_number") %>%
                dplyr::mutate(exon_number = as.numeric(exon_number)) %>%
                as.data.frame()
        }

        # For ensemble and ncbi the "chr" needs to be added to the chromosome number
        if (is.na(stringr::str_extract(formattedGTF$chrom[1], "chr"))) {
            formattedGTF$chrom <- paste0("chr", formattedGTF$chrom)
        }
        # Report mitocondrial chromosome as chrM (like in BSgenome) and not chrMT. 
        formattedGTFfinal <- formattedGTF %>%
            dplyr::mutate(
                chrom = ifelse(chrom == 'chrMT', 'chrM', chrom))
    } else{
        formattedGTF <- data.frame()
        cat("Specify path to GTF file.\n")
    }

    return(formattedGTF)
}


# Read GTF file
.readGTF<- function(pathToGTF = NULL){
    # suppress warning for closing unused connection
    gtf <- suppressWarnings(rtracklayer::import(pathToGTF)) %>%
        as.data.frame()

    # Keep only the needed rows
    # Rename seqnames to chrom
    gtf <- gtf %>%
        dplyr::filter(type == "exon") %>%
        dplyr::rename(chrom = seqnames) %>%
        dplyr::mutate(strand = as.character(strand),
            chrom = as.character(chrom))
    return(gtf)
}

# get needed column in the GTF file
.getNeededColumn<- function(){
    needColumns <-
        c(
            "chrom",
            "start",
            "end",
            "width",
            "strand",
            "type",
            "gene_name",
            "transcript_id"
        )
    return(needColumns)
}

# If the function you are looking for is not here check supportFunction.R
# Functions in supportFunction.R are used by multiple functions.
