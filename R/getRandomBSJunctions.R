#' @title Retrieve random back-spliced junctions
#'
#' @description The function getRandomBSJunctions() retrieves random
#' back-spliced junctions from the user genome annotation.
#'
#' @param gtf A dataframe containing genome annotation information This can be
#' generated with \code{\link{formatGTF}}.
#'
#' @param n Integer specifying the number of randomly selected transcripts
#' from which random back-spliced junctions are extracted. Default value = 100.
#'
#' @param f An integer specifying the fraction of single exon circRNAs that
#' have to be present in the output data frame. Default value is 10.
#' 
#' @param setSeed An integer which is used for selecting random back-spliced 
#' junctions. Default values is set to NULL. 
#'
#' @return A data frame.
#'
#' @examples
#' # Load short version of the gencode v19 annotation file
#' data("gtf")
#'
#' # Get 10 random back-spliced junctions
#' randomBSJunctions <- getRandomBSJunctions(gtf, n = 10, f = 10)
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @import dplyr
#' @export
getRandomBSJunctions <- function(gtf, n = 100, f = 10, setSeed=NULL) {
   
    if(!is.null(setSeed)){
        seed<-setSeed
    }else{
        seed = NULL
    }
      # Create an empty data frame
    randomBSJunctions <-.createRandomBSJunctionsDF(n)

   # Select random BSEs from gtf
    allBSEs <- .selectRandomBSEs(gtf, n, f,seed)

    # For negative strand
    allBSEsNeg <- allBSEs[allBSEs$strand == "-",]
    if (nrow(allBSEsNeg) > 0) {
        randomBSJunctions$gene[seq_along(allBSEsNeg$exon_number)] <-
            allBSEsNeg$gene_name
        randomBSJunctions$strand[seq_along(allBSEsNeg$exon_number)] <-
            allBSEsNeg$strand
        randomBSJunctions$chrom[seq_along(allBSEsNeg$exon_number)] <-
            allBSEsNeg$chrom
        randomBSJunctions$startUpBSE[seq_along(allBSEsNeg$exon_number)] <-
            allBSEsNeg$end
        randomBSJunctions$endDownBSE[seq_along(allBSEsNeg$exon_number)] <-
            allBSEsNeg$start1
    }

    # For positive strand
    allBSEsPos <- allBSEs[allBSEs$strand == "+",]
    if (nrow(allBSEsPos) > 0) {
        randomBSJunctions$gene[(nrow(allBSEsNeg) + 1):n] <-
            allBSEsPos$gene_name
        randomBSJunctions$strand[(nrow(allBSEsNeg) + 1):n] <-
            allBSEsPos$strand
        randomBSJunctions$chrom[(nrow(allBSEsNeg) + 1):n] <-
            allBSEsPos$chrom
        randomBSJunctions$startUpBSE[(nrow(allBSEsNeg) + 1):n] <-
            allBSEsPos$start
        randomBSJunctions$endDownBSE[(nrow(allBSEsNeg) + 1):n] <-
            allBSEsPos$end1
    }

    # Generate a unique identifier
    randomBSJunctions$id <- .getID(randomBSJunctions)
    return(randomBSJunctions)
}


# Create randomBSJunctions data frame
.createRandomBSJunctionsDF <- function(n=100){
    basicColumns <- .getBasicColNames()

    # Create an empty data frame
    randomBSJunctions <-
        data.frame(matrix(nrow = n, ncol = length(basicColumns)))
    colnames(randomBSJunctions) <- basicColumns
    return(randomBSJunctions)
}


# Select random BSEs from gtf file
.selectRandomBSEs<- function(gtf, n=100, f = 10,seed=NULL){
    
    if(!is.null(seed)){
      base::set.seed(seed)
    }
  # calculate the percentage of back-spliced junctions from single exons
  c <- round((n / 100) * f, 0)
  # Select two random back-spliced exons (up and down) from n randomly
  # selected transcript
  bsExons2 <- gtf %>%
    dplyr::filter(type == "exon") %>%
    dplyr::group_by(transcript_id) %>%
    dplyr::filter(n() >= 3) %>%
    dplyr::ungroup() %>%
    dplyr::filter(transcript_id %in%
                    sample(unique(transcript_id), n - c)) %>%
    dplyr::group_by(transcript_id) %>%
    dplyr::sample_n(2) %>%
    dplyr::arrange(exon_number)
  
    if(c !=0){
      # Select one random exon from n randomly selected transcript
      bsExons1 <- gtf %>%
        dplyr::filter(type == "exon") %>%
        dplyr::group_by(transcript_id) %>%
        dplyr::filter(n >= 3) %>%
        dplyr::ungroup() %>%
        dplyr::filter(transcript_id %in% base::sample(unique(transcript_id), c)) %>%
        dplyr::group_by(transcript_id) %>%
        dplyr::sample_n(1) %>%
        dplyr::arrange(exon_number)
      
      # If only one exon is picked then the row is duplicated. This step is only
      # made to simplify the code below without introducing an additional if
      # statement when a single exon is present in the isoform sampled in
      # the above step.
      bsExons1 <- dplyr::bind_rows(bsExons1, bsExons1)
      # Join the 2 data frame bsExons1 and bsExons2
      bsExons12 <-  dplyr::bind_rows(bsExons1, bsExons2)
    }else{
      bsExons12 <-  bsExons2
    }
 
    # Align duplicates on the same row
    indexDup <- which(duplicated(bsExons12$transcript_id))
    duplicates <- bsExons12[indexDup,]
    colnames(duplicates)<- paste0(colnames(duplicates),'1')
    
    cleanedDF <- bsExons12[-indexDup,]
    mt <- match(cleanedDF$transcript_id, duplicates$transcript_id1)

    allBSEs <- dplyr::bind_cols(cleanedDF, duplicates[mt,])

    return(allBSEs)
}


# If the function you are looking for is not here check supportFunction.R
# Functions in supportFunction.R are used by multiple functions.
