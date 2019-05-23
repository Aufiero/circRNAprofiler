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
#' (\url{http://www.mirbase.org/ftp.shtml}). If TRUE is specifed then the
#' latest release is automatically downloaded. If FALSE is specified then a
#' file named mature.fa containing fasta format sequences of all mature miRNA
#' sequences previously downloaded by the user from mirBase must be present
#' in the working directory. Default value is TRUE.
#'
#' @param totalMatches An integer specifying the total number of matches that
#' have to be found between the seed region of the miR and the seed site of the
#' target sequence. If the total number of matches found is less than the
#' cut-off, the seed site is discarded. Default value is 7.
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
#' # Annotate the first 3 back-spliced junctions
#' annotatedBSJs <- annotateBSJs(mergedBSJunctions[1:3, ], gtf,
#'     isRandom = FALSE)
#'
#' # Get genome
#' genome <- BSgenome::getBSgenome("BSgenome.Hsapiens.UCSC.hg19")
#'
#' # Retrieve target sequences.
#' targets <- getCircSeqs(
#'     annotatedBSJs,
#'     gtf,
#'     genome)
#'
#' # Screen target sequence for miR binding sites.
#' pathToMiRs <- system.file("extdata", "miRs.txt",package="circRNAprofiler")
#'
#' #miRsites <- getMiRsites(
#' #    targets,
#' #    miRspeciesCode = "hsa",
#' #    miRBaseLatestRelease = TRUE,
#' #    totalMatches = 6,
#' #    maxNonCanonicalMatches = 1,
#' #    pathToMiRs)
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom seqinr swap
#' @importFrom stringr str_count
#' @import dplyr
#' @export
getMiRsites <- function(targets,
    miRspeciesCode = "hsa",
    miRBaseLatestRelease = TRUE,
    totalMatches = 7,
    maxNonCanonicalMatches = 1,
    pathToMiRs = NULL) {
    # 271 options are possible the argument miRspeciesCode
    miRspeciesCodes <- circRNAprofiler::miRspeciesCodes
    if (!miRspeciesCode %in% miRspeciesCodes$code) {
        stop(
            "miR species code not correct, type data(miRspeciesCode) to see
            the available codes and the corresponding species."
        )
    }
    # Retrieve miR sequences from miRBase or read from a given file
    # that must be named mature.fa
    microRNAs <-
        getMiRseqs(miRBaseLatestRelease,
            miRspeciesCode,
            pathToMiRs)

    if (length(targets) == 1 & names(targets)[[1]] == "circ") {
        # Create list to store miR results
        miRsites <- createMiRsitesList(targets, microRNAs)

        # Find the matching regions.
        # To identify possible matches between the miR seqs and the
        # circRNA seqs, each circRNA is scanned, using a sliding window of 1,
        # for the presence of putative matches with the reverse
        # of the miR sequences.

        # Scan the target sequences, by using a sliding windows of 1
        for (i in seq_along(miRsites$targets$seq)) {
            # Analysis starts at:
            analysisStart <- 30

            for (j in seq_along(miRsites$microRNAs$seqRev)) {
                # miR sequence to analyze
                mirSeq <-  miRsites$microRNAs$seqRev[j]
                seedSeq <-
                    base::substr(mirSeq, nchar(mirSeq) - 7, nchar(mirSeq) - 1)

                # The seed region of the miR sequences is searched starting
                # from position 30 of the target sequence. This is done because
                # a stretch of nucleotides needs to be left to search for AU
                # content upstream the seed site and for possible matches with
                # the central and compensatoy regions of the miR sequence.
                # Initialize the variables

                indexTargetSeq <- analysisStart
                counts <- 0
                k <- 1
                totMatchesInSeed <- as.character()
                cwcMatchesInSeed <- as.character()
                seedLocation <- as.character()
                t1 <- as.character()
                totMatchesInCentral <- as.character()
                cwcMatchesInCentral <- as.character()
                totMatchesInCompensatory <- as.character()
                cwcMatchesInCompensatory <- as.character()
                localAUcontent <- as.character()

                # The length of the circRNAs is greater because there are
                # additional nucleotides (50) taken from the first
                # back-spliced exon. This is done so that also the
                # back-spliced junction is analyzed. By using a sliding
                # window of 1 nucleotide the circRNA sequence is scanned
                # and searched for the presence of the miR seed sequence.

                while (indexTargetSeq <=
                        (miRsites$targets$length[i] + (analysisStart -  1))) {
                    # Allow  asymmetric mismatches/ buldge in the sequences

                    # FIND MATCHES WITH THE SEED REGION
                    circSeqForSeed <-
                        base::substr(miRsites$targets$seq[i],
                            indexTargetSeq,
                            indexTargetSeq + nchar(seedSeq))

                    seedMatches <-
                        getMatches(circSeqForSeed,
                            seedSeq,
                            isSeed = TRUE)

                    # Filter based on the number fo non-canonical matches
                    # allowed and keep the max total matches found
                    seedMatches <- seedMatches %>%
                        dplyr::filter(.data$ncm <= maxNonCanonicalMatches) %>%
                        dplyr::arrange(desc(.data$tm)) %>%
                        dplyr::slice(1) %>%
                        dplyr::rename(
                            maxTotalMatches = .data$tm,
                            maxContinuousMatches = .data$cwcm)

                    # Check seed type
                    # If there are at least totalMatches then retrive all
                    # the other info.
                    if (nrow(seedMatches) != 0) {
                        if (seedMatches$maxTotalMatches >= totalMatches) {
                            totMatchesInSeed[k] <- seedMatches$maxTotalMatches

                            # Write mer if there are at least 2 continuous
                            # wc matches
                            if (seedMatches$maxContinuousMatches > 1) {
                                cwcMatchesInSeed[k] <-
                                    paste(seedMatches$maxContinuousMatches,
                                        "mer", sep = "")
                            } else{
                                cwcMatchesInSeed[k] <-
                                    seedMatches$maxContinuousMatches
                            }

                            t1[k] <-
                                base::substr(circSeqForSeed,
                                    nchar(circSeqForSeed),
                                    nchar(circSeqForSeed))

                            seedLocation[k] <- indexTargetSeq

                            # FIND MATCHES WITH THE CENTRAL REGION
                            circSeqForcentral <-
                                base::substr(miRsites$targets$seq[i],
                                    indexTargetSeq - 4,
                                    indexTargetSeq - 1)

                            centralSeq <-
                                base::substr(mirSeq,
                                    nchar(mirSeq) - 11,
                                    nchar(mirSeq) - 8)

                            centralMatches <-
                                getMatches(circSeqForcentral,
                                    centralSeq,
                                    isSeed = FALSE)

                            # Filter based on the number of non-canonical
                            # matches allowed and keep the max total matches
                            # found
                            centralMatches <- centralMatches %>%
                                dplyr::arrange(desc(.data$tm)) %>%
                                dplyr::slice(1) %>%
                                dplyr::rename(
                                    maxTotalMatches = .data$tm,
                                    maxContinuousMatches = .data$cwcm
                                )

                            if (nrow(centralMatches) != 0) {
                                totMatchesInCentral[k] <-
                                    centralMatches$maxTotalMatches
                                cwcMatchesInCentral[k] <-
                                    centralMatches$maxContinuousMatches
                            } else{
                                totMatchesInCentral[k] <- 0
                                cwcMatchesInCentral[k] <- 0
                            }

                            # FIND MATCHES WITH THE COMPENSATORY REGION
                            circSeqForcompensatory <-
                                base::substr(
                                    miRsites$targets$seq[i],
                                    indexTargetSeq - 8,
                                    indexTargetSeq - 5
                                )

                            compensatorySeq <-
                                base::substr(mirSeq,
                                    nchar(mirSeq) - 15,
                                    nchar(mirSeq) - 12)

                            compensatoryMatches <-
                                getMatches(circSeqForcompensatory,
                                    compensatorySeq,
                                    isSeed = FALSE)

                            # Filter based on the number fo non-canonical
                            # matches allowed and keep the max total matches
                            # found
                            compensatoryMatches <-
                                compensatoryMatches %>%
                                dplyr::arrange(desc(.data$tm)) %>%
                                dplyr::slice(1) %>%
                                dplyr::rename(
                                    maxTotalMatches = .data$tm,
                                    maxContinuousMatches = .data$cwcm)

                            if (nrow(compensatoryMatches) != 0) {
                                totMatchesInCompensatory[k] <-
                                    compensatoryMatches$maxTotalMatches
                                cwcMatchesInCompensatory[k] <-
                                    compensatoryMatches$maxContinuousMatches
                            } else{
                                totMatchesInCompensatory[k] <- 0
                                cwcMatchesInCompensatory[k] <- 0
                            }

                            # FOR LOCAL AU CONTENT around the seed region 2-8
                            downstreamRegion <-
                                base::substr(
                                    miRsites$targets$seq[i],
                                    indexTargetSeq - 10,
                                    indexTargetSeq - 1)
                            upstreamRegion <-
                                base::substr(
                                    miRsites$targets$seq[i],
                                    indexTargetSeq + nchar(seedSeq[1]),
                                    indexTargetSeq + nchar(seedSeq[1]) + 9)

                            localAUcontent[k] <-
                                Biostrings::letterFrequency(Biostrings::RNAStringSet(
                                    paste0(upstreamRegion, downstreamRegion)),
                                    letters = "AU") / 20

                            k <- k + 1

                            # Count seed site
                            counts <- counts + 1
                        }
                    }

                    # Sliding window - increase 1 nucleotide at a time
                    indexTargetSeq <- indexTargetSeq + 1
                }

                # The target sequence was analyzed for the presence of all
                # miR sequences and now the information needs to be stored
                # into the corresponding dataframes

                # Keep the following info if at least 1 seed match is found
                if (counts >= 1) {
                    # Number of seed sites found for each mir
                    miRsites$counts[i, miRsites$microRNAs$id[j]] <-
                        counts
                    # Location of the seed sites found for each mir
                    miRsites$seedLocation[i, miRsites$microRNAs$id[j]] <-
                        paste(seedLocation, collapse = ",")
                    # Check t1 nucleotide
                    miRsites$t1[i, miRsites$microRNAs$id[j]] <-
                        paste(t1, collapse = ",")
                    # Number of total matches found with the seed region of
                    # the miR
                    miRsites$totMatchesInSeed[i, miRsites$microRNAs$id[j]] <-
                        paste(totMatchesInSeed, collapse = ",")
                    # Number of continuous WC matches found with the seed
                    # region of the miR
                    miRsites$cwcMatchesInSeed[i, miRsites$microRNAs$id[j]] <-
                        paste(cwcMatchesInSeed, collapse = ",")
                    # Number of total matches found with central region of
                    # the miR
                    miRsites$totMatchesInCentral[i, miRsites$microRNAs$id[j]] <-
                        paste(totMatchesInCentral, collapse = ",")
                    # Number of continuous WC matches found with the central
                    # region of the miR
                    miRsites$cwcMatchesInCentral[i, miRsites$microRNAs$id[j]] <-
                        paste(cwcMatchesInCentral, collapse = ",")
                    # Number of total matches found with compensatory region
                    # of the miR
                    miRsites$totMatchesInCompensatory[i, miRsites$microRNAs$id[j]] <-
                        paste(totMatchesInCompensatory, collapse = ",")
                    # Number of of continuous WC matches found with
                    # compensatory region of the miR
                    miRsites$cwcMatchesInCompensatory[i, miRsites$microRNAs$id[j]] <-
                        paste(cwcMatchesInCompensatory, collapse = ",")
                    # Content of AU nucletide in the 10 nucleotide upstream and
                    # downstream the seed region
                    miRsites$localAUcontent[i, miRsites$microRNAs$id[j]] <-
                        paste(localAUcontent, collapse = ",")

                }

            }
        }

    } else{
        stop("target sequences not valid, only circRNA sequences
            are allowed.")
    }

    return(miRsites)
    }


#' @title Retrive miRNA sequences
#'
#' @description The function getMiRseqs() retrieves miRNA sequences from
#' miRbase data base
#'
#' @param miRspeciesCode A string specifying the species code (3 letters) as
#' reported in miRBase. E.g. to analyze the mouse microRNAs specify "mmu",
#' to analyze the human microRNAs specify "hsa". Type data(miRspeciesCode)
#' to see the available codes and the corresponding species reported
#' in miRBase 22 release. Default value is "hsa".
#'
#' @param miRBaseLatestRelease A logical specifying whether to download the
#' latest release of the mature sequences of the microRNAs from miRBase
#' (\url{http://www.mirbase.org/ftp.shtml}). If TRUE is specifed then the
#' latest release is automatically downloaded. If FALSE is specified then a
#' file named mature.fa containing fasta format sequences of all mature miRNA
#' sequences previously downloaded by the user from mirBase must be present
#' in the working directory. Default value is TRUE.
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
#' empty, all miRs of the species are retrieved.
#'
#' @return A data frame.
#'
#' @keywords internal
#'
#' @examples
#' # Inner function
#' microRNAs <- getMiRseqs(miRBaseLatestRelease = TRUE,
#'     miRspeciesCode = "hsa")
#'
#' @importFrom IRanges reverse
#' @importFrom readr read_tsv
#' @importFrom magrittr %>%
#' @importFrom utils read.table
#' @export
getMiRseqs <- function(miRBaseLatestRelease = TRUE,
    miRspeciesCode = "hsa",
    pathToMiRs = NULL) {
    # 271 options are possible the argument miRspeciesCode
    miRspeciesCodes <- circRNAprofiler::miRspeciesCodes
    if (!miRspeciesCode %in% miRspeciesCodes$code) {
        stop(
            "miR species code not correct, type data(miRspeciesCode) to see
            the available codes and the corresponding species."
        )
    }

    # Read experiment information
    options(readr.num_columns = 0)


    # 2) Retrieve mir sequences from miRBase or read from a given file
    # that must be named mature.fa
    if (miRBaseLatestRelease) {
        miRBase <-
            readr::read_tsv("ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz",
                col_names = FALSE)
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

    # The first line shuld start with ">"
    firstChar <- base::substr(miRBase$seq[1], 1, 1)
    if (firstChar != ">") {
        stop(
            "The microRNAs in mature.fa must be reported in fasta
            format. The first row must start with >miRBaseID then the
            second row must contain the miR sequence."
        )
    }

    # Clean mir names from unwanted characters
    miRBaseCleaned <- gsub(" .*", "", miRBase$seq)

    if (is.null(pathToMiRs)) {
        pathToMiRs <- "miRs.txt"
    }


    if (file.exists(pathToMiRs)) {
        # Read the user given miR sequences
        miRsFromFile <-
            utils::read.table(pathToMiRs,
                header = TRUE,
                stringsAsFactors = FALSE)

        # colnames(miRsFromFile)[1] <- "id"
    } else{
        miRsFromFile <- data.frame()

    }


    if (nrow(miRsFromFile) == 0) {
        cat("miRs.txt is empty or absent. All miRNAs of the
            specified species will be analyzed")
    }

    # microRNAs used in the analysis
    if (nrow(miRsFromFile) > 0) {
        # The first line should start with ">"
        firstChar <- base::substr(miRsFromFile$id[1], 1, 1)
        if (firstChar != ">") {
            stop("The microRNA ids in miRs.txt must start with " > ".
                E.g. >miRid")
        }

        microRNAids <-
            miRBaseCleaned[miRBaseCleaned %in% miRsFromFile$id]

        if (length(miRsFromFile$id) - length(microRNAids) > 0) {
            notFound <-
                miRsFromFile$id[!(miRsFromFile$id %in% microRNAids)]
            cat(paste("miRs not found:",
                paste(notFound, collapse = ", ")))

        }
        } else{
            microRNAids <-
                miRBaseCleaned[grep(miRspeciesCode, miRBaseCleaned)]

        }

    # Fill the microRNAids data frame in the miRsites list
    microRNAs  <-
        data.frame(matrix(nrow = length(microRNAids), ncol = 4))
    colnames(microRNAs) <-
        c("id", "length", "seq", "seqRev")


    for (i in seq_along(microRNAs$id)) {
        indexMir <-
            which(miRBaseCleaned == microRNAids[i])

        microRNAs$id[i] <- microRNAids[i]
        microRNAs$seq[i] <-
            miRBaseCleaned[indexMir + 1]

        # The direction of the sequences of the microRNAs downloaded
        # from miRBase are from 5' to 3'. To be able to find matching
        # regions within the target sequences also reported in direction
        # 5' to 3', the reverse of the mir sequences must be used.

        microRNAs$seqRev[i] <-
            IRanges::reverse(miRBaseCleaned[indexMir + 1])


        microRNAs$length[i] <-
            nchar(microRNAs$seq[i])
    }
    return(microRNAs)
    }



# Create list to store miR results
createMiRsitesList <- function(targets, microRNAs) {
    if (length(targets) == 1 & names(targets)[[1]] == "circ") {
        miRsites <- vector("list", 12)
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


        # Fill the target dataframe in the miRsites list
        miRsites$targets <- targets[[1]]

        # Fill the microRNAs data frame in the miRsites list
        miRsites$microRNAs  <-  microRNAs

        # Create the empty dataframe to store the number of miR sites found
        # in the target sequences
        miRsites$counts <-
            data.frame(matrix(
                nrow = nrow(miRsites$targets),
                ncol = length(microRNAs$id) + 1
            ))
        colnames(miRsites$counts) <- c("id", microRNAs$id)
        miRsites$counts$id <- miRsites$targets$id

        # Create the empty dataframe to store the info about the total
        # matches found between the seed site of the target and the seed
        # region of the miR
        miRsites$totMatchesInSeed <-
            data.frame(matrix(
                nrow = nrow(miRsites$targets),
                ncol = length(microRNAs$id) + 1
            ))
        colnames(miRsites$totMatchesInSeed) <-
            c("id", microRNAs$id)
        miRsites$totMatchesInSeed$id <- miRsites$targets$id

        # Create the empty dataframe to store the info about the continuous
        # WC matches found between the seed site of the target and the
        # seed region of the miR
        miRsites$cwcMatchesInSeed <-
            data.frame(matrix(
                nrow = nrow(miRsites$targets),
                ncol = length(microRNAs$id) + 1
            ))
        colnames(miRsites$cwcMatchesInSeed) <-
            c("id", microRNAs$id)
        miRsites$cwcMatchesInSeed$id <- miRsites$targets$id

        # Create the empty dataframe to store the location of the seed match
        miRsites$seedLocation <-
            data.frame(matrix(
                nrow = nrow(miRsites$targets),
                ncol = length(microRNAs$id) + 1
            ))
        colnames(miRsites$seedLocation) <- c("id", microRNAs$id)
        miRsites$seedLocation$id <- miRsites$targets$id

        # Create the empty dataframe to store the info about the t1
        # nucleotide
        miRsites$t1 <-
            data.frame(matrix(
                nrow = nrow(miRsites$targets),
                ncol = length(microRNAs$id) + 1
            ))
        colnames(miRsites$t1) <- c("id", microRNAs$id)
        miRsites$t1$id <- miRsites$targets$id

        # Create the empty dataframe to store the total matches found
        # with the central region of the miR
        miRsites$totMatchesInCentral <-
            data.frame(matrix(
                nrow = nrow(miRsites$targets),
                ncol = length(microRNAs$id) + 1
            ))
        colnames(miRsites$totMatchesInCentral) <-
            c("id", microRNAs$id)
        miRsites$totMatchesInCentral$id <- miRsites$targets$id

        # Create the empty dataframe to store continuous WC matches
        # found with the central region of the miR
        miRsites$cwcMatchesInCentral <-
            data.frame(matrix(
                nrow = nrow(miRsites$targets),
                ncol = length(microRNAs$id) + 1
            ))
        colnames(miRsites$cwcMatchesInCentral) <-
            c("id", microRNAs$id)
        miRsites$cwcMatchesInCentral$id <- miRsites$targets$id

        # Create the empty dataframe to store total matches found with the
        # compensatory region
        miRsites$totMatchesInCompensatory <-
            data.frame(matrix(
                nrow = nrow(miRsites$targets),
                ncol = length(microRNAs$id) + 1
            ))
        colnames(miRsites$totMatchesInCompensatory) <-
            c("id", microRNAs$id)
        miRsites$totMatchesInCompensatory$id <-
            miRsites$targets$id

        # Create the empty dataframe to store continuous WC matches
        # found with the compensatory region
        miRsites$cwcMatchesInCompensatory <-
            data.frame(matrix(
                nrow = nrow(miRsites$targets),
                ncol = length(microRNAs$id) + 1
            ))
        colnames(miRsites$cwcMatchesInCompensatory) <-
            c("id", microRNAs$id)
        miRsites$cwcMatchesInCompensatory$id <-
            miRsites$targets$id

        # Create the empty dataframe to store the info about the
        # percentage of A/U nucleotides found in the 10 nucleotides
        # upstream and downstream the seed match
        miRsites$localAUcontent <-
            data.frame(matrix(
                nrow = nrow(miRsites$targets),
                ncol = length(microRNAs$id) + 1
            ))
        colnames(miRsites$localAUcontent) <- c("id", microRNAs$id)
        miRsites$localAUcontent$id <- miRsites$targets$id
    } else{
        stop("target sequences not valid, only circRNA sequences
            are allowed.")
    }

    return(miRsites)

    }



# The function getMatches() retrives the number of total matches,
# continuous canonical Watson and Crick matches and a non-canonical matches
# (G:U) between 2 given sequences of the same length. The analysis is
# performed by comparing the 2 sequences without gap and by introducing a
# gap in each position of one of the 2 sequences at the time.
getMatches <-  function(seq1, seq2, isSeed = TRUE) {
    # Exclude the t1 nucleotide from the analysis if the seq2 is the
    # seed
    if (isSeed) {
        seq1 <-
            base::substr(seq1, 1, nchar(seq2))
    }

    # Introduce a gap WITHIN the sequence.
    # This is done to analyze also whether the presence of a single buldge
    # WITHIN one of the 2 sequence increase the number of matches

    seq1[2] <- paste0("-", seq1)
    seq1Char <- unlist(base::strsplit(seq1[2], ""))

    seq2[2] <- paste0("-", seq2)
    seq2Char <- unlist(base::strsplit(seq2[2], ""))

    i <- 1
    # seq2Char and seq1Char have the same length
    # We do not consider the gap in the last position since it correspond to a
    # scift of 1 nucleotide. It will be analyzed when we increase the
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

    # Initialize the dataframe
    # cwcm stands for continuous watson-crik matches
    # tm stands for total matches
    # ncm stands for non canonical matches
    matches <-  data.frame(matrix(nrow = 1 + ((length(
        seq1
    ) - 2) + (length(
        seq2
    ) - 2)), ncol = 3))
    colnames(matches) <- c("cwcm", "tm", "ncm")


    # Compare the 2 sequences
    # n is a non canonical pair
    # w is a watson- crick match
    # m is a mistmatch
    comparedSeq <-
        compareSequences(seq1[1], seq2[1], isGUMatch = FALSE)
    continuousMatches <-
        max(nchar(base::strsplit(comparedSeq, "m")[[1]]))
    matches$cwcm[1] <- continuousMatches

    comparedSeq <-
        compareSequences(seq1[1], seq2[1], isGUMatch = TRUE)
    totalMatches <-
        stringr::str_count(comparedSeq, "w") +
        stringr::str_count(comparedSeq, "n")
    nonCanonicalpairs <- stringr::str_count(comparedSeq, "n")
    matches$tm[1] <- totalMatches
    matches$ncm[1] <- nonCanonicalpairs

    # Find the max number of total matches and the max number of continuous
    # matches beteeen the seed region and the seed site considering also the
    # presence of a buldge within the seed region of the mir sequence
    i <- 2
    k <- 2
    while (i < length(seq1)) {
        comparedSeq <-
            compareSequences(seq1[2], seq2[i + 1], isGUMatch = FALSE)
        continuousMatches <-
            max(nchar(base::strsplit(comparedSeq, "m")[[1]]))
        matches$cwcm[k] <- continuousMatches

        comparedSeq <-
            compareSequences(seq1[2], seq2[i + 1], isGUMatch = TRUE)
        totalMatches <-
            stringr::str_count(comparedSeq, "w") + stringr::str_count(comparedSeq, "n")
        matches$tm[k] <- totalMatches
        nonCanonicalpairs <- stringr::str_count(comparedSeq, "n")
        matches$ncm[k] <- nonCanonicalpairs

        i <- i + 1
        k <- k + 1
    }

    # Find the max number of total matches and the max number of continuous
    # matches beteeen the seed region and the seed site considering also
    # the presence of a buldge within the seed site of the target sequence.
    i <- 2
    while (i < length(seq2)) {
        comparedSeq <-
            compareSequences(seq2[2], seq1[i + 1], isGUMatch = FALSE)
        continuousMatches <-
            max(nchar(base::strsplit(comparedSeq, "m")[[1]]))
        matches$cwcm[k] <- continuousMatches

        comparedSeq <-
            compareSequences(seq1[2], seq2[i + 1], isGUMatch = TRUE)
        totalMatches <-
            stringr::str_count(comparedSeq, "w") + stringr::str_count(comparedSeq, "n")
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
compareSequences <- function(seq1, seq2, isGUMatch = TRUE) {
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



#' @title Rearrange miR results
#'
#' @description The function rearrangeMiRres() rearranges the results of the
#' getMiRsites() function. Each element of the list contains the miR resutls
#' relative to one circRNA.
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
#' # Annotate the first 3 back-spliced junctions
#' annotatedBSJs <- annotateBSJs(mergedBSJunctions[1:3, ], gtf,
#'     isRandom = FALSE)
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
#' #miRsites <- getMiRsites(
#' #    targets,
#' #    miRspeciesCode = "hsa",
#' #    miRBaseLatestRelease = TRUE,
#' #    totalMatches = 6,
#' #    maxNonCanonicalMatches = 1,
#' #    pathToMiRs)
#'
#' # Rearrange miR results
#' #rearragedMiRres <- rearrangeMiRres(miRsites)
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
        rearragedMiRres[[i]][[1]] <-
            data.frame(matrix(
                nrow = 1,
                ncol = length(getTargetsColNames())
            ))
        colnames(rearragedMiRres[[i]][[1]]) <-
            colnames(miRsites$targets)
        rearragedMiRres[[i]][[1]] <- miRsites$targets[i, ]

        # Second dataframe is filled with the results
        rearragedMiRres[[i]][[2]] <-
            data.frame(matrix(
                nrow = nrow(miRsites$microRNAs),
                ncol = 3 + length(names(miRsites)[3:9])
            ))
        colnames(rearragedMiRres[[i]][[2]]) <-
            c("miRid", "miRseq", "miRseqRev", names(miRsites)[3:9])


        rearragedMiRres[[i]][[2]]$miRid <- miRsites$microRNAs$id
        rearragedMiRres[[i]][[2]]$miRseq <-
            miRsites$microRNAs$seq
        rearragedMiRres[[i]][[2]]$miRseqRev <-
            miRsites$microRNAs$seqRev

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

        rearragedMiRres[[i]][[2]]$t1 <-
            as.character(miRsites$t1[i, idM])

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

        #Remove miR with no counts
        rearragedMiRres[[i]][[2]] <-
            rearragedMiRres[[i]][[2]] %>%
            dplyr::filter(counts  != 0)
    }

    return(rearragedMiRres)
}
