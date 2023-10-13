#' @title Screen target sequences for miR binding sites
#'
#' @description The function getMirSites() searches miRNA binding sites within
#' the circRNA sequences. The user can restrict the analisis only to a subset
#' of miRs. In this case, miR ids must go in miR.txt file. If the latter
#' is absent or empty, all miRs of the specified miRspeciesCode are considered
#' in the analysis.
#'
#' @param targets A list containing the target sequences to analyze.
#' This data frame can be generated with \code{\link{getCircSeqs}}.
#'
#' @param miRspeciesCode A string specifying the species code (3 letters) as
#' reported in miRBase db. E.g. to analyze the mouse microRNAs specify "mmu",
#' to analyze the human microRNAs specify "hsa". Type data(miRspeciesCode)
#' to see the available codes and the corresponding species reported
#' in miRBase 22 release. Default value is "hsa".
#'
#' @param miRBaseLatestRelease A logical specifying whether to download the
#' latest release of the mature sequences of the microRNAs from miRBase
#' (\url{http://www.mirbase.org/ftp.shtml}). If TRUE is specified then the
#' latest release is automatically downloaded. If FALSE is specified then a
#' file named mature.fa containing fasta format sequences of all mature miRNA
#' sequences previously downloaded by the user from mirBase must be present
#' in the working directory. Default value is TRUE.
#'
#' @param totalMatches An integer specifying the total number of matches that
#' have to be found between the seed region of the miR and the seed site of the
#' target sequence. If the total number of matches found is less than the
#' cut-off, the seed site is discarded. The maximun number of possible matches
#' is 7. Default value is 7.
#'
#' @param maxNonCanonicalMatches An integer specifying the max number of
#' non-canonical matches (G:U) allowed between the seed region of the miR and
#' the seed site of the target sequence. If the max non-canonical matches found
#' is greater than the cut-off, the seed site is discarded. Default value is 1.
#'
#' @param pathToMiRs A string containing the path to the miRs.txt file.
#' The file miRs.txt contains the microRNA ids from miRBase
#' specified by the user. It must have one column with header id. The first row
#' must contain the miR name starting with the ">", e.g >hsa-miR-1-3p. The
#' sequences of the miRs will be automatically retrieved from the mirBase latest
#' release or from the given mature.fa file, that should be present in the
#' working directory. By default pathToMiRs is set to NULL and the file it is
#' searched in the working directory. If miRs.txt is located in a different
#' directory then the path needs to be specified. If this file is absent or
#' empty, all miRs of the specified species are considered in the
#' miRNA analysis.
#'
#' @return A list.
#'
#' @examples
#' # Load a data frame containing detected back-spliced junctions
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
#' # Retrieve target sequences.
#' targets <- getCircSeqs(
#'     annotatedBSJs,
#'     gtf,
#'     genome)
#'
#' # Screen target sequence for miR binding sites.
#' pathToMiRs <- system.file("extdata", "miRs.txt", package="circRNAprofiler")
#'
#' #miRsites <- getMiRsites(
#' #   targets,
#' #   miRspeciesCode = "hsa",
#' #   miRBaseLatestRelease = TRUE,
#' #   totalMatches = 6,
#' #   maxNonCanonicalMatches = 1,
#' #   pathToMiRs )
#' }
#' 
#'
#' @importFrom IRanges reverse
#' @importFrom readr read_tsv
#' @importFrom utils read.table
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom seqinr swap
#' @importFrom stringr str_count
#' @importFrom utils data
#' @import dplyr
#' @export
getMiRsites <- function(targets, miRspeciesCode = "hsa",
    miRBaseLatestRelease = TRUE, totalMatches = 7,
    maxNonCanonicalMatches = 1, pathToMiRs = NULL) {


    # Retrieve miR sequences from miRBase or from mature.fa file
    microRNAs <- .getMiRseqs(miRBaseLatestRelease, miRspeciesCode, pathToMiRs)

    if (length(targets) == 1 & names(targets)[[1]] == "circ") {
        # Create list to store miR results
        miRsites <- .createMiRsitesList(targets, microRNAs)
        # Scan the target sequences, by using a sliding windows of 1
        for (i in seq_along(miRsites$targets$id)) {
            cat(paste0('\nAnalysing: ', miRsites$targets$id[i]))
            analysisStart <- 30
            circSeq <- miRsites$targets$seq[i]
            circLen <- miRsites$targets$length[i]
            for (j in seq_along(miRsites$microRNAs$id)) {
                # miR sequence to analyze
                mirSeq <-  miRsites$microRNAs$seqRev[j]
                # The analysis start at position 30 to retrieve upstream features
                indexTargetSeq <- analysisStart
                # Create createMiRsitesTempDF data frame
                tempDF <- .createMiRsitesTempDF()
                k <- 1
                # CircRNAs seq is X nucleotides longer than the predicted length.
                while (indexTargetSeq <= (circLen + (analysisStart -  1))) {
                    # FIND MATCHES WITH THE SEED REGION
                    seedMatches <- .getSeedMatches(circSeq,
                            indexTargetSeq,
                            mirSeq,
                            maxNonCanonicalMatches)
                    # Retrieve all info if seedMatches are found
                    if (nrow(seedMatches) != 0) {
                        if (seedMatches$maxTotalMatches >= totalMatches) {
                            tempDF[k, ] <- .createMiRsitesTempDF()
                            tempDF$totMatchesInSeed[k] <- seedMatches$maxTotalMatches

                            # check if there are continuous WC matches
                            tempDF$cwcMatchesInSeed[k] <-.checkMer(seedMatches)
                            # Keep t1 nucleotide info
                            tempDF$t1[k] <-  .getT1nucleotide(circSeq, indexTargetSeq, mirSeq)
                            # Keep location info
                            tempDF$seedLocation[k] <- indexTargetSeq
                            # FIND MATCHES WITH THE CENTRAL REGION
                            centralMatches <-  .getCentralMatches(circSeq, indexTargetSeq, mirSeq)
                            tempDF$totMatchesInCentral[k] <-  centralMatches$maxTotalMatches
                            tempDF$cwcMatchesInCentral[k] <- centralMatches$maxContinuousMatches
                            # FIND MATCHES WITH THE COMPENSATORY REGION
                            compensatoryMatches <- .getCompensatoryMatches(circSeq, indexTargetSeq, mirSeq)
                            tempDF$totMatchesInCompensatory[k] <- compensatoryMatches$maxTotalMatches
                            tempDF$cwcMatchesInCompensatory[k] <-  compensatoryMatches$maxContinuousMatches
                            # FOR LOCAL AU CONTENT around the seed region 2-8
                            tempDF$localAUcontent[k] <-  .getAUcontent(circSeq, indexTargetSeq, mirSeq)
                            # Count seed site
                            tempDF$count[k] <- 1
                            k <- k + 1
                        }
                    }
                    indexTargetSeq <- indexTargetSeq + 1 # increase of 1
                }
                # Keep the following info if at least 1 seed match is found
                if (sum(tempDF$count) >= 1) {
                    mirId <- miRsites$microRNAs$id[j]
                    rowId <- i
                    miRsites <- .fillMiRsites(miRsites, rowId, mirId, tempDF)
                }
            }
        }
    } else{
        stop("target sequences not valid, only circRNA sequences are allowed.")
    }
    return(miRsites)
}


# Retrive miRNA sequences
.getMiRseqs <- function(miRBaseLatestRelease = TRUE,
    miRspeciesCode = "hsa",
    pathToMiRs = NULL) {
    # 271 options are possible the argument miRspeciesCode
    .checkMiRspeciesCode(miRspeciesCode)

    # Read experiment information
    options(readr.num_columns = 0)
    # Retrieve miR sequences from miRBase or read from a mature.fa file
    miRBase <- .getMiRsFromMiRBase(miRBaseLatestRelease)

    # Get miRNAs ids to analyze
    microRNAids <- .getMiRids(miRBase, pathToMiRs, miRspeciesCode)
    # Fill the microRNAids data frame in the miRsites list
    microRNAs  <- data.frame(matrix(nrow = length(microRNAids), ncol = 4))
    colnames(microRNAs) <-  c("id", "length", "seq", "seqRev")
    for (i in seq_along(microRNAs$id)) {
        indexMir <- which(miRBase == microRNAids[i])
        microRNAs$id[i] <- microRNAids[i]
        microRNAs$seq[i] <- miRBase[indexMir + 1]

        # The direction of the sequences of the microRNAs downloaded
        # from miRBase are from 5' to 3'. To be able to find matching
        # regions within the target sequences also reported in direction
        # 5' to 3', the reverse of the mir sequences must be used.
        microRNAs$seqRev[i] <- IRanges::reverse(miRBase[indexMir + 1])
        microRNAs$length[i] <-  nchar(microRNAs$seq[i])
    }
    return(microRNAs)
}


# Check miR species code
.checkMiRspeciesCode <- function(miRspeciesCode){
    miRspeciesCodes <- NULL
    data("miRspeciesCodes", package= "circRNAprofiler",envir = environment())
    if (!miRspeciesCode %in% miRspeciesCodes$code) {
        stop(
            "miR species code not correct, type data(miRspeciesCode) to see
            the available codes and the corresponding species."
        )
    }
}

#Retrieve mir sequences from miRBase or read from a  mature.fa file
.getMiRsFromMiRBase <- function(miRBaseLatestRelease = TRUE){

    if (miRBaseLatestRelease) {
        url<-"ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz"

        miRBase <-
            tryCatch(
            readr::read_tsv(url, col_names = FALSE), error = function(e) NULL)
        if(is.null(miRBase)){
            stop(
                "Unable to establish a connection with mirbase.org.
                Try later."
            )
        }

    } else if (file.exists("mature.fa")) {
        miRBase <-
            utils::read.table(
                "mature.fa",
                header = FALSE,
                stringsAsFactors = FALSE,
                sep = "\t"
            )
    } else {
        stop(
            "mature.fa not present in the project folder.
            Set miRBaseLatestRelease = TRUE to automatically download
            the latest release of the mature sequences of the
            microRNAs from miRBase"
        )
    }

   
    colnames(miRBase)[1] <- "seq"

    firstChar <- base::substr(miRBase$seq[1], 1, 1)
    if (firstChar != ">") {
        stop(
            "The microRNAs in mature.fa must be reported in fasta
            format. The first row must start with >miRBaseID then the
            second row must contain the miR sequence."
        )
    }

    miRBaseCleaned <- gsub(" .*", "", miRBase$seq)
    return(miRBaseCleaned)
}



# Retrive miRNA to analyze
.getMiRids <- function(miRBase, pathToMiRs = NULL, miRspeciesCode){
    # Read miRs.txt
    miRsFromFile <- .readMiRs(pathToMiRs)
    if (nrow(miRsFromFile) == 0) {
        cat("miRs.txt is empty or absent. All miRNAs of the
            specified species will be analyzed")
    }
    if (nrow(miRsFromFile) > 0) {
        # The first line should start with ">"
        firstChar <- base::substr(miRsFromFile$id[1], 1, 1)
        if (firstChar != ">") {
            stop("The microRNA ids in miRs.txt must start with " > ".
                E.g. >miRid")
        }
        microRNAids <-  miRBase[miRBase %in% miRsFromFile$id]

        if (length(miRsFromFile$id) - length(microRNAids) > 0) {
            notFound <- miRsFromFile$id[!(miRsFromFile$id %in% microRNAids)]
            cat(paste("miRs not found:", paste(notFound, collapse = ", ")))
        }
        } else{
            microRNAids <- miRBase[grep(miRspeciesCode, miRBase)]
        }
    return(microRNAids)
}



# Create list to store miR results
.createMiRsitesList <- function(targets, microRNAs) {
    if (length(targets) == 1 & names(targets)[[1]] == "circ") {
        miRsites <- vector("list", 12)
        miRsites <- .nameMiRsitesListDF(miRsites)
        # Fill the target dataframe in the miRsites list
        miRsites$targets <- targets[[1]]
        # Fill the microRNAs data frame in the miRsites list
        miRsites$microRNAs  <-  microRNAs
        # Dataframe to store the number of miR binding sites
        miRsites$counts <- createMiRsitesInternalDF(miRsites, microRNAs)
        colnames(miRsites$counts) <- c("id", microRNAs$id)
        miRsites$counts$id <- miRsites$targets$id
        # Dataframe to store the total matches
        miRsites$totMatchesInSeed <- createMiRsitesInternalDF(miRsites, microRNAs)
        colnames(miRsites$totMatchesInSeed) <- c("id", microRNAs$id)
        miRsites$totMatchesInSeed$id <- miRsites$targets$id
        # Dataframe to store the continuous WC matches
        miRsites$cwcMatchesInSeed <- createMiRsitesInternalDF(miRsites, microRNAs)
        colnames(miRsites$cwcMatchesInSeed) <- c("id", microRNAs$id)
        miRsites$cwcMatchesInSeed$id <- miRsites$targets$id
        # Dataframe to store the location of the seed match
        miRsites$seedLocation <- createMiRsitesInternalDF(miRsites, microRNAs)
        colnames(miRsites$seedLocation) <- c("id", microRNAs$id)
        miRsites$seedLocation$id <- miRsites$targets$id
        # Dataframe to store the info about the t1 nucleotide
        miRsites$t1 <- createMiRsitesInternalDF(miRsites, microRNAs)
        colnames(miRsites$t1) <- c("id", microRNAs$id)
        miRsites$t1$id <- miRsites$targets$id
        # Dataframe to store the total matches with the central region
        miRsites$totMatchesInCentral <- createMiRsitesInternalDF(miRsites, microRNAs)
        colnames(miRsites$totMatchesInCentral) <- c("id", microRNAs$id)
        miRsites$totMatchesInCentral$id <- miRsites$targets$id
        # Dataframe to store continuous WC matches with the central region
        miRsites$cwcMatchesInCentral <- createMiRsitesInternalDF(miRsites, microRNAs)
        colnames(miRsites$cwcMatchesInCentral) <- c("id", microRNAs$id)
        miRsites$cwcMatchesInCentral$id <- miRsites$targets$id
        # Dataframe to store total matches found with the compensatory region
        miRsites$totMatchesInCompensatory <- createMiRsitesInternalDF(miRsites, microRNAs)
        colnames(miRsites$totMatchesInCompensatory) <- c("id", microRNAs$id)
        miRsites$totMatchesInCompensatory$id <- miRsites$targets$id
        # Dataframe to store continuous WC matches with the compensatory region
        miRsites$cwcMatchesInCompensatory <- createMiRsitesInternalDF(miRsites, microRNAs)
        colnames(miRsites$cwcMatchesInCompensatory) <- c("id", microRNAs$id)
        miRsites$cwcMatchesInCompensatory$id <- miRsites$targets$id
        # Dataframe to store the percentage of A/U nucleotides
        miRsites$localAUcontent <-  createMiRsitesInternalDF(miRsites, microRNAs)
        colnames(miRsites$localAUcontent) <- c("id", microRNAs$id)
        miRsites$localAUcontent$id <- miRsites$targets$id
    } else{
        stop("target sequences not valid, only circRNA sequences are allowed.")
    }
    return(miRsites)
}


# get names miRsiteslist
.nameMiRsitesListDF <- function(miRsites){
    names(miRsites)[1] <- "targets"
    names(miRsites)[2] <- "microRNAs"
    names(miRsites)[3] <- "counts"
    names(miRsites)[4] <- "totMatchesInSeed"
    names(miRsites)[5] <- "cwcMatchesInSeed"
    names(miRsites)[6] <- "seedLocation"
    names(miRsites)[7] <- "t1"
    names(miRsites)[8] <- "totMatchesInCentral"
    names(miRsites)[9] <- "cwcMatchesInCentral"
    names(miRsites)[10] <- "totMatchesInCompensatory"
    names(miRsites)[11] <- "cwcMatchesInCompensatory"
    names(miRsites)[12] <- "localAUcontent"
    return(miRsites)
}

# Create miRsites internal data frame
createMiRsitesInternalDF <- function(miRsites, microRNAs) {
    intDF <- data.frame(matrix(
        nrow = nrow(miRsites$targets),
        ncol = length(microRNAs$id) + 1
    ))

    return(intDF)
}


# Create createMiRsitesTempDF data frame
.createMiRsitesTempDF <- function() {
    tempDF <-
        data.frame(matrix(nrow = 1,
            ncol = 10))
    colnames(tempDF) <-
        c(
            "count",
            "totMatchesInSeed",
            "cwcMatchesInSeed",
            "seedLocation",
            "t1",
            "totMatchesInCentral",
            "cwcMatchesInCentral",
            "totMatchesInCompensatory",
            "cwcMatchesInCompensatory",
            "localAUcontent"
        )

    tempDF$count <- 0
    return(tempDF)
}

# The function getMatches() retrives the number of total matches,
# continuous canonical Watson and Crick matches and a non-canonical matches
# (G:U) between 2 given sequences of the same length. The analysis is
# performed by comparing the 2 sequences without gap and by introducing a
# gap in each position of one of the 2 sequences at the time.
.getMatches <-  function(seq1, seq2, isSeed = TRUE) {
    # Exclude the t1 nucleotide from the analysis if the seq2 is the seed
    if (isSeed) {
        seq1 <-  base::substr(seq1, 1, nchar(seq2))
    }
    # Introduce a gap WITHIN the sequence.
    seq1[2] <- paste0("-", seq1)
    seq1Char <- unlist(base::strsplit(seq1[2], ""))

    seq2[2] <- paste0("-", seq2)
    seq2Char <- unlist(base::strsplit(seq2[2], ""))

    i <- 1
    # seq2Char and seq1Char have the same length
    # We do not consider the gap in the last position since it correspond to a
    # shift of 1 nucleotide. It will be analyzed when we increase the
    # indexTargetSeq.
    while (i < (length(seq2Char) - 1)) {
        # Swap the  gap position "-" within the sequence of the microRNA
        swap(seq2Char[i], seq2Char[i + 1])
        seq2[i + 2] <-
            paste(seq2Char, collapse = "")

        # Swap the  gap position "-" within the sequence of the circRNA
        swap(seq1Char[i], seq1Char[i + 1])
        seq1[i + 2] <- paste(seq1Char, collapse = "")

        i <- i + 1
    }

    # Compare the 2 sequences and store the results in matches data frame
    matches <- .fillMatches(seq1, seq2)
    return(matches)
}

# Create matches data frame
.createMatchesDF <- function(seq1, seq2){
    # cwcm stands for continuous watson-crik matches tm stands for total matches
    # ncm stands for non canonical matches
    matches <-  data.frame(matrix(nrow = 1 + ((length(
        seq1
    ) - 2) + (length(
        seq2
    ) - 2)), ncol = 3))
    colnames(matches) <- c("cwcm", "tm", "ncm")
    return(matches)
}

# Compare the 2 sequences and find the max number of total matches and the
# max number of continuous matches between the seed region and the seed site
# considering also the presence of a buldge within the seed region of the miR
# sequence or the target sequences. Store the information of each iteration
# matches data frame.
.fillMatches<- function(seq1, seq2){
    matches <- .createMatchesDF(seq1, seq2)
    # n is a non canonical pair, w is a  watson- crick match, m is a mistmatch
    comparedSeq <- .compareSequences(seq1[1], seq2[1], isGUMatch = FALSE)
    continuousMatches <- max(nchar(base::strsplit(comparedSeq, "m")[[1]]))
    matches$cwcm[1] <- continuousMatches

    comparedSeq <- .compareSequences(seq1[1], seq2[1], isGUMatch = TRUE)
    totalMatches <- stringr::str_count(comparedSeq, "w") +
        stringr::str_count(comparedSeq, "n")
    nonCanonicalpairs <- stringr::str_count(comparedSeq, "n")
    matches$tm[1] <- totalMatches
    matches$ncm[1] <- nonCanonicalpairs

    # Find the max number of total matches and the max number of continuous
    # matches considering also the presence of a buldge
    i <- 2
    k <- 2
    while (i < length(seq1)) {
        comparedSeq <- .compareSequences(seq1[2], seq2[i + 1], isGUMatch = FALSE)
        continuousMatches <- max(nchar(base::strsplit(comparedSeq, "m")[[1]]))
        matches$cwcm[k] <- continuousMatches

        comparedSeq <- .compareSequences(seq1[2], seq2[i + 1], isGUMatch = TRUE)
        totalMatches <- stringr::str_count(comparedSeq, "w") + stringr::str_count(comparedSeq, "n")
        matches$tm[k] <- totalMatches
        nonCanonicalpairs <- stringr::str_count(comparedSeq, "n")
        matches$ncm[k] <- nonCanonicalpairs
        i <- i + 1
        k <- k + 1
    }
    # Find the max number of total matches and continuous matches considering
    # also  the presence of a buldge within the seed site of the target sequence.
    i <- 2
    while (i < length(seq2)) {
        comparedSeq <- .compareSequences(seq2[2], seq1[i + 1], isGUMatch = FALSE)
        continuousMatches <- max(nchar(base::strsplit(comparedSeq, "m")[[1]]))
        matches$cwcm[k] <- continuousMatches

        comparedSeq <- .compareSequences(seq1[2], seq2[i + 1], isGUMatch = TRUE)
        totalMatches <- stringr::str_count(comparedSeq, "w") + stringr::str_count(comparedSeq, "n")
        matches$tm[k] <- totalMatches
        nonCanonicalpairs <- stringr::str_count(comparedSeq, "n")
        matches$ncm[k] <- nonCanonicalpairs
        i <- i + 1
        k <- k + 1
    }
    return(matches)
}


# The function compareSequences() combines two equal-length
# strings that are assumed to be aligned into a single character string.
# The Watson-Crick matches are replaced with "w", the mismatches with "m",
# G:U matches are replaced either with "n" or with "w" based on the value of
# isGUMatch argument.
.compareSequences <- function(seq1, seq2, isGUMatch = TRUE) {
    # n is a non canonical pair
    # w is a watson- crick match
    # m is a mistmatch
    GU <- as.character()
    if (isGUMatch) {
        GU <- "n"
    } else{
        GU <- "m"
    }

    charsSeq1 <- unlist(base::strsplit(seq1, ""))
    charsSeq2 <- unlist(base::strsplit(seq2, ""))
    newSeq <- as.character()

    for (i in seq_along(charsSeq1)) {
        if ((charsSeq1[i] == "A" &  charsSeq2[i] == "U") |
                (charsSeq1[i] == "U" &  charsSeq2[i] == "A") |
                (charsSeq1[i] == "C" &  charsSeq2[i] == "G") |
                (charsSeq1[i] == "G" &  charsSeq2[i] == "C")) {
            newSeq[i] <- "w"
        } else if ((charsSeq1[i] == "G" &  charsSeq2[i] == "U")  |
                (charsSeq1[i] == "U" & charsSeq2[i] == "G")) {
            newSeq[i] <- GU
        } else {
            newSeq[i] <- "m"
        }
    }

    comparedSeq <- paste(newSeq, collapse = "")

    return(comparedSeq)
}

# Find match between miR seed region and the miR seed region binding site
.getSeedMatches <-
    function(circSeq,
        indexTargetSeq,
        mirSeq,
        maxNonCanonicalMatches) {
        seedSeq <-
            base::substr(mirSeq, nchar(mirSeq) - 7, nchar(mirSeq) - 1)

        circSeqForSeed <-
            base::substr(circSeq,
                indexTargetSeq,
                indexTargetSeq + nchar(seedSeq))

        seedMatches <-
            .getMatches(circSeqForSeed,
                seedSeq,
                isSeed = TRUE)

        # Filter based on the number of non-canonical matches
        # allowed and keep the max total matches found
        seedMatchesFiltered <- seedMatches %>%
            dplyr::filter(ncm <= maxNonCanonicalMatches) %>%
            dplyr::arrange(desc(tm)) %>%
            dplyr::slice(1) %>%
            dplyr::rename(maxTotalMatches = tm,
                maxContinuousMatches = cwcm)


        return(seedMatchesFiltered)
    }


# check if there are continuous WC matches
# Write mer if there are at least 2 continuous
# wc matches
.checkMer <- function(seedMatches) {
    cwcMatchesInSeed <- as.character()
    if (seedMatches$maxContinuousMatches > 1) {
        cwcMatchesInSeed <- paste(seedMatches$maxContinuousMatches,
            "mer",
            sep = "")
    } else{
        cwcMatchesInSeed <-  seedMatches$maxContinuousMatches
    }
    return(cwcMatchesInSeed)
}

# Get t1 nucleotide
.getT1nucleotide <- function(circSeq, indexTargetSeq, mirSeq) {
    seedSeq <-
        base::substr(mirSeq, nchar(mirSeq) - 7, nchar(mirSeq) - 1)

    circSeqForSeed <-
        base::substr(circSeq,
            indexTargetSeq,
            indexTargetSeq + nchar(seedSeq))

    t1 <-
        base::substr(circSeqForSeed,
            nchar(circSeqForSeed),
            nchar(circSeqForSeed))
    return(t1)
}




# Find match between miR central region and the sequence upstream the
# miR seed region binding site
.getCentralMatches <- function(circSeq, indexTargetSeq, mirSeq) {
    circSeqForCentral <-
        base::substr(circSeq,
            indexTargetSeq - 4,
            indexTargetSeq - 1)

    centralSeq <-
        base::substr(mirSeq,
            nchar(mirSeq) - 11,
            nchar(mirSeq) - 8)

    centralMatches <-
        .getMatches(circSeqForCentral,
            centralSeq,
            isSeed = FALSE)

    # Filter based on the number of non-canonical
    # matches allowed and keep the max total matches
    # found
    centralMatchesFiltered <- centralMatches %>%
        dplyr::arrange(desc(tm)) %>%
        dplyr::slice(1) %>%
        dplyr::rename(maxTotalMatches = tm,
            maxContinuousMatches = cwcm)

    return(centralMatchesFiltered)
}



# Find match between miR compensatory region and the sequence upstream the
# miR seed region binding site
.getCompensatoryMatches <-
    function(circSeq, indexTargetSeq, mirSeq) {
        circSeqForCompensatory <-
            base::substr(circSeq,
                indexTargetSeq - 8,
                indexTargetSeq - 5)

        compensatorySeq <-
            base::substr(mirSeq,
                nchar(mirSeq) - 15,
                nchar(mirSeq) - 12)

        compensatoryMatches <-
            .getMatches(circSeqForCompensatory,
                compensatorySeq,
                isSeed = FALSE)

        # Filter based on the number fo non-canonical
        # matches allowed and keep the max total matches
        # found
        compensatoryMatchesFiltered <-
            compensatoryMatches %>%
            dplyr::arrange(desc(tm)) %>%
            dplyr::slice(1) %>%
            dplyr::rename(maxTotalMatches = tm,
                maxContinuousMatches = cwcm)

        return(compensatoryMatchesFiltered)

    }

# Find the local AU content around the miR seed region binding site
.getAUcontent <- function(circSeq, indexTargetSeq, mirSeq) {
    seedSeq <-
        base::substr(mirSeq, nchar(mirSeq) - 7, nchar(mirSeq) - 1)
    upstreamRegion <-
        base::substr(circSeq,
            indexTargetSeq - 10,
            indexTargetSeq - 1)
    downstreamRegion <-
        base::substr(circSeq,
            indexTargetSeq + nchar(seedSeq[1]),
            indexTargetSeq + nchar(seedSeq[1]) + 9)

    localAUcontent <-
        Biostrings::letterFrequency(Biostrings::RNAStringSet(paste0(upstreamRegion,
            downstreamRegion)),
            letters = "AU") / 20

    return(localAUcontent)
}


# Fill dataframe in miRsites list
.fillMiRsites <- function(miRsites, rowId, mirId, tempDF) {
    # Number of miR seed binding sites
    miRsites$counts[rowId, mirId] <- sum(tempDF$count)
    # Location of the miR seed binding sites
    miRsites$seedLocation[rowId, mirId] <-
        paste(tempDF$seedLocation, collapse = ",")
    # Check t1 nucleotide
    miRsites$t1[rowId, mirId] <- paste(tempDF$t1, collapse = ",")
    # Number of total matches with miR seed region
    miRsites$totMatchesInSeed[rowId, mirId] <-
        paste(tempDF$totMatchesInSeed, collapse = ",")
    # Number of continuous WC matches miR seed region
    miRsites$cwcMatchesInSeed[rowId, mirId] <-
        paste(tempDF$cwcMatchesInSeed, collapse = ",")
    # Number of total matches with miR central region
    miRsites$totMatchesInCentral[rowId, mirId] <-
        paste(tempDF$totMatchesInCentral, collapse = ",")
    # Number of continuous WC matches with miR central region
    miRsites$cwcMatchesInCentral[rowId, mirId] <-
        paste(tempDF$cwcMatchesInCentral, collapse = ",")
    # Number of total matches with miR compensatory region
    miRsites$totMatchesInCompensatory[rowId, mirId] <-
        paste(tempDF$totMatchesInCompensatory,
            collapse = ",")
    # Number of of continuous WC matches with miR compensatory region
    miRsites$cwcMatchesInCompensatory[rowId, mirId] <-
        paste(tempDF$cwcMatchesInCompensatory,
            collapse = ",")
    # Content of AU nucleotide around the miR seed binding sites
    miRsites$localAUcontent[rowId, mirId] <-
        paste(tempDF$localAUcontent, collapse = ",")
    return(miRsites)
}


#' @title Rearrange miR results
#'
#' @description The function rearrangeMiRres() rearranges the results of the
#' getMiRsites() function. Each element of the list contains the miR results
#' relative to one circRNA. For each circRNA only miRNAs for which at least 1
#' miRNA binding site is found are reported.
#'
#' @param  miRsites A list containing the miR sites found in the RNA
#' target sequence. it can be generated with \code{\link{getMiRsites}}.
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
#' genome <- BSgenome::getBSgenome("BSgenome.Hsapiens.UCSC.hg19")
#'
#' # Retrieve target sequences.
#'  targets <- getCircSeqs(
#'      annotatedBSJs,
#'      gtf,
#'      genome)
#'
#' # Screen target sequence for miR binding sites.
#' pathToMiRs <- system.file("extdata", "miRs.txt",package="circRNAprofiler")
#'
#' # miRsites <- getMiRsites(
#' #    targets,
#' #    miRspeciesCode = "hsa",
#' #    miRBaseLatestRelease = TRUE,
#' #    totalMatches = 6,
#' #    maxNonCanonicalMatches = 1,
#' #    pathToMiRs)
#'
#' # Rearrange miR results
#' # rearragedMiRres <- rearrangeMiRres(miRsites)
#'
#'
#' @import dplyr
#' @importFrom magrittr %>%
#' @export
rearrangeMiRres <- function(miRsites) {
    # Create an empty list of 9 data frames
    rearragedMiRres <- vector("list", nrow(miRsites$targets))
    for (i in seq_along(rearragedMiRres)) {
        rearragedMiRres[[i]] <- vector("list", 2)
        # First dataframe is filled with the target information
        rearragedMiRres[[i]][[1]] <- data.frame(matrix(
                nrow = 1,
                ncol = length(.getTargetsColNames())
            ))
        colnames(rearragedMiRres[[i]][[1]]) <- colnames(miRsites$targets)
        rearragedMiRres[[i]][[1]] <- miRsites$targets[i,]
        # Second dataframe is filled with the results
        rearragedMiRres[[i]][[2]] <- data.frame(matrix(
                nrow = nrow(miRsites$microRNAs),
                ncol = 3 + length(names(miRsites)[3:9])
            ))
        colnames(rearragedMiRres[[i]][[2]]) <-
            c("miRid", "miRseq", "miRseqRev", names(miRsites)[3:9])

        rearragedMiRres[[i]][[2]]$miRid <- miRsites$microRNAs$id
        rearragedMiRres[[i]][[2]]$miRseq <-miRsites$microRNAs$seq
        rearragedMiRres[[i]][[2]]$miRseqRev <- miRsites$microRNAs$seqRev

        idM <- miRsites$microRNAs$id
        counts <- as.character(miRsites$counts[i, idM])
        counts[which(is.na(miRsites$counts[i, idM]))] <- 0
        rearragedMiRres[[i]][[2]]$counts <- as.numeric(counts)

        rearragedMiRres[[i]][[2]]$totMatchesInSeed <-
            as.character(miRsites$totMatchesInSeed[i, idM])

        rearragedMiRres[[i]][[2]]$cwcMatchesInSeed <-
            as.character(miRsites$cwcMatchesInSeed[i, idM])

        rearragedMiRres[[i]][[2]]$seedLocation <-
            as.character(miRsites$seedLocation[i, idM])

        rearragedMiRres[[i]][[2]]$t1 <- as.character(miRsites$t1[i, idM])

        rearragedMiRres[[i]][[2]]$totMatchesInCentral <-
            as.character(miRsites$totMatchesInCentral[i, idM])
        rearragedMiRres[[i]][[2]]$cwcMatchesInCentral <-
            as.character(miRsites$cwcMatchesInCentral[i, idM])

        rearragedMiRres[[i]][[2]]$totMatchesInCompensatory <-
            as.character(miRsites$totMatchesInCompensatory[i, idM])
        rearragedMiRres[[i]][[2]]$cwcMatchesInCompensatory <-
            as.character(miRsites$cwcMatchesInCompensatory[i, idM])

        rearragedMiRres[[i]][[2]]$localAUcontent <-
            as.character(miRsites$localAUcontent[i, idM])

        #Remove miR with zero counts
        rearragedMiRres[[i]][[2]] <- rearragedMiRres[[i]][[2]] %>%
            dplyr::filter(counts  != 0)
    }
    return(rearragedMiRres)
}


# If the function you are looking for is not here check supportFunction.R
# Functions in supportFunction.R are used by multiple functions.
