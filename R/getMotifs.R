#' @title Screen target sequences for recurrent motifs
#'
#' @description The function getMotifs() scans the target sequences for the
#' presence of recurrent motifs of a specific length defined in input.
#' By setting rbp equals to TRUE, the identified motifs are matched with motifs
#' of known RNA Binding Proteins (RBPs) deposited in the ATtRACT database and
#' with motifs specified by the user. The user motifs must go in the file
#' motifs.txt. If this file is absent or  empty, only motifs from the ATtRACT
#' database are considered in the analysis. By setting rbp equals to FALSE,
#' only motifs that do not match with any motifs deposited in the ATtRACT
#' database or user motifs are reported in the final output. Location of the
#' selected motifs is also reported. This corresponds to the start position of
#' the motif within the sequence (1-index based).
#'
#' @param targets A list containing the target sequences to analyze.
#' It can be generated with \code{\link{getCircSeqs}},
#' \code{\link{getSeqsAcrossBSJs}} or \code{\link{getSeqsFromGRs}}.
#'
#' @param width An integer specifying the length of all possible motifs to
#' extract from the target sequences. Default value is 6.
#'
#' @param species A string specifying the species of the RBP motifs to use.
#' Type data(attractSpecies) to see the possible options.
#' Default value is "Hsapiens".
#'
#' @param rbp A logical specifying whether to report only motifs matching
#' with known RBP motifs from ATtRACT database or user motifs specified in
#' motifs.txt. If FALSE is specified only motifs that do not match with any of
#' these motifs are reported. Default values is TRUE.
#'
#' @param reverse A logical specifying whether to reverse the motifs collected
#' from ATtRACT database and from motifs.txt. If TRUE is specified all the
#' motifs are reversed and analyzed together with the direct motifs as they are
#' reported in the ATtRACT db and motifs.txt. Default value is FALSE.
#'
#' @param pathToMotifs A string containing the path to the motifs.txt
#' file. The file motifs.txt contains motifs/regular expressions specified
#' by the user. It must have 3 columns with headers:
#' - id (1st column): name of the motif. - e.g. RBM20 or motif1.
#' - motif (2nd column): motif/pattern to search.
#' - length (3rd column): length of the motif.
#' By default pathToMotifs is set to NULL and the file it is searched in the
#' working directory. If motifs.txt is located in a different directory then
#' the path needs to be specified. If this file is absent or empty only the
#' motifs of RNA Binding Proteins in the ATtRACT database are considered in
#' the motifs analysis.
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
#' # Example with the first back-spliced junctions.
#' # Multiple back-spliced junctions can also be analyzed at the same time.
#'
#' # Annotate detected back-spliced junctions
#' annotatedBSJs <- annotateBSJs(mergedBSJunctions[1:10, ], gtf, isRandom = FALSE)
#'
#' # Retrieve target sequences.
#' targets <- getSeqsFromGRs(
#'     annotatedBSJs,
#'     lIntron = 200,
#'     lExon = 10,
#'     type = "ie",
#'     species = "Hsapiens",
#'     genome = "hg19")
#'
#' # Get motifs
#' motifs <- getMotifs(
#'     targets,
#'     width = 6,
#'     species = "Hsapiens",
#'     rbp = TRUE,
#'     reverse = FALSE)
#'
#'
#' @importFrom stringr str_count
#' @importFrom stringi stri_locate_all
#' @importFrom Biostrings RNAStringSet
#' @export
getMotifs <-
    function(targets,
        width = 6,
        species  = "Hsapiens",
        rbp = TRUE,
        reverse = FALSE,
        pathToMotifs = NULL) {
        # Compute all motifs of a given length (width) and keep the ones
        # present at least one time in the target sequences

        computedMotifs <- computeMotifs(targets, width)

        # Filter motifs matching with RBPs
        filteredMotifs <-
            filterMotifs(computedMotifs, species, rbp, reverse, pathToMotifs)

        # Crete list to store motif reults
        motifs <- createMotifsList(targets, filteredMotifs)

        for (i in seq_along(motifs)) {
            # Put the string NA so that the we do not get error when calling
            # RNAStringSet
            motifs[[i]]$targets$seq[is.na(motifs[[i]]$targets$seq)] <-
                "NA"

            # Adjust sequences - (!!! See also computeMotifs())
            if (names(motifs)[1] == "circ") {
                # We adjust the length of circ seqs  based on the length of
                # the motifs and remove the first 30 nucleotides since
                # they are analyzed at the end with together with
                # the BSJs (see getCircSeq() function)
                ss <- base::substring(motifs[[i]]$targets$seq, 30,
                    ((motifs[[i]]$targets$length + 30) + (width - 1)))

                # Put the string NA so that the we do not get error when
                # calling RNAStringSet
                ss[is.na(ss)] <- "NA"
                rnaSS <- Biostrings::RNAStringSet(ss)

            } else if (names(motifs)[1] == "bsj") {
                # We adjust the length of BSJ seqs  based on the length of
                # the motifs in a way that the motifs must have at least
                # 1 nucleotide crossing the BSJ.
                r1 <-
                    (base::nchar(motifs[[i]]$targets$seq) / 2) - (width - 2)
                r2 <-
                    (base::nchar(motifs[[i]]$targets$seq) / 2) + (width - 1)
                rnaSS <-
                    Biostrings::RNAStringSet(base::substring(motifs[[i]]$targets$seq, r1, r2))

            } else{
                rnaSS <- Biostrings::RNAStringSet(motifs[[i]]$targets$seq)
            }

            # Count occurence and location of the filtered motifs

            for (j in seq_along(motifs[[i]]$motif$motif)) {
                stringMotif <- motifs[[i]]$motif$motif[j]
                # Find location. Consider also overlapping
                # patterns (?=pattern)
                locations <-
                    stringi::stri_locate_all(rnaSS,
                        regex = paste("(?=", stringMotif, ")", sep = ""))

                if (names(motifs)[1] == "circ") {
                    motifs[[i]]$locations[stringMotif] <-
                        as.character(lapply(
                            locations,
                            FUN = function(locations)
                                paste(locations[, 1] + 29, collapse = ",")
                        ))
                } else{
                    motifs[[i]]$locations[stringMotif] <-
                        as.character(lapply(
                            locations,
                            FUN = function(locations)
                                paste(locations[, 1], collapse = ",")
                        ))
                }

                # Count occurences
                motifs[[i]]$counts[stringMotif] <-
                    stringr::str_count(rnaSS,
                        paste("(?=", stringMotif, ")", sep = ""))

            }

            # Put 0 count where the string "NA" is found
            #index <- which(motifs[[i]]$targets$seq == "NA")
            #motifs[[i]]$counts[index,-1] <- 0

        }


        return(motifs)

    }


#' @title Compute and search motifs of a given length in the target sequences.
#'
#' @description The function computeMotifs() first computes all possible motifs
#' of a given length and subsequently searches those motifs in the target
#' sequences. Only motifs present at least 1 time are kept.
#'
#' @param targets A list containing the target sequences to analyze.
#' It can be generated with: \code{\link{getCircSeqs}} or
#' \code{\link{getSeqsAcrossBSJs}} or \code{\link{getSeqsFromGRs}}.
#'
#' @param width An integer specifying the length of all possible motifs to
#' extract from the target sequences. Default value is 6.
#'
#' @return A data frame.
#'
#' @examples
#' # Load data frame containing detected back-spliced junctions
#' data("mergedBSJunctions")
#'
#' # Load short version of the gencode v19 annotation file
#' data("gtf")
#'
#' # Example with the first back-spliced junctions.
#' # Multiple back-spliced junctions can also be analyzed at the same time.
#'
#' # Annotate detected back-spliced junctions
#' annotatedBSJs <- annotateBSJs(mergedBSJunctions[1, ], gtf)
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
#' # Inner function
#' # Compute motifs
#' computedMotifs <- computeMotifs(
#'     targets,
#'     width = 6)
#'
#'
#' @import magrittr
#' @importFrom Biostrings RNAStringSet
#' @importFrom Biostrings oligonucleotideFrequency
#' @export
computeMotifs <-
    function(targets, width = 6) {
        if (width < 4 | width > 20) {
            stop("width must be an integer between 4 and 20")
        }

        # Compute all motifs of length width
        computedMotifs <- as.character()
        for (i in seq_along(targets)) {
            # Put the string NA so that the we do not get error when calling
            # RNAStringSet
            targets[[i]]$seq[is.na(targets[[i]]$seq)] <-
                "NA"

            # Adjust sequences
            if (names(targets)[1] == "circ") {
                # We adjust the length of circ seqs  based on the length of
                # the motifs and remove the first 30 nucleotides since
                # they are analyzed at the end with together with
                # the BSJs (see getCircSeq() function)
                ss <- base::substring(targets[[i]]$seq, 30,
                    ((targets[[i]]$length + 30) + (width - 1)))

                # Put the string NA so that the we do not get error when
                # calling RNAStringSet
                ss[is.na(ss)] <- "NA"
                rnaSS <- Biostrings::RNAStringSet(ss)

            } else if (names(targets)[1] == "bsj") {
                # We adjust the length of BSJ seqs  based on the length of
                # the motifs in a way that the motifs must have at least
                # 1 nucleotide crossing the BSJ.
                r1 <-
                    (base::nchar(targets[[i]]$seq) / 2) - (width - 2)
                r2 <-
                    (base::nchar(targets[[i]]$seq) / 2) + (width - 1)
                rnaSS <-
                    Biostrings::RNAStringSet(base::substring(targets[[i]]$seq, r1, r2))

            } else{
                rnaSS <- Biostrings::RNAStringSet(targets[[i]]$seq)
            }

            computedMotifs <- c(
                computedMotifs,
                Biostrings::oligonucleotideFrequency(rnaSS, width, step = 1, simplify.as =
                        "collapsed") %>%
                    .[. > 0] %>%
                    names()
            ) %>%
                unique()
        }


        return(computedMotifs)

    }


#' @title Retrieve ATtRACT motifs
#'
#' @description The function getUserAttractMotifs() reads the user motifs
#' in motifs.txt and retrieves the motifs deposited in the ATtRACT database
#' (\url{http://attract.cnic.es}).
#'
#' @param species A string specifying the species of the RBP motifs to use.
#' Type data(attractSpecies) to see the possible options.
#' Default value is "Hsapiens".
#'
#' @param reverse A logical specifying whether to reverse the motifs collected
#' from ATtRACT database and from motifs.txt. If TRUE is specified all the
#' motifs are reversed and analyzed together with the direct motifs as they are
#' reported in the ATtRACT db and motifs.txt. Default value is FALSE.
#'
#' @param pathToMotifs A string containing the path to the motifs.txt
#' file. The file motifs.txt contains motifs/regular expressions specified
#' by the user. It must have 3 columns with headers:
#' - id (1st column): name of the motif. - e.g. RBM20 or motif1.
#' - motif (2nd column): motif/pattern to search.
#' - length (3rd column): length of the motif.
#' By default pathToMotifs is set to NULL and the file it is searched in the
#' working directory. If motifs.txt is located in a different directory then
#' the path needs to be specified. If this file is absent or empty only the
#' motifs of RNA Binding Proteins in the ATtRACT database are considered in
#' the motifs analysis.
#'
#' @return A data frame.
#'
#' @examples
#' # Inner function
#' userAttractMotifs <- getUserAttractMotifs(species = "Hsapiens",
#'     reverse = FALSE)
#'
#' @importFrom readr read_tsv
#' @importFrom utils download.file
#' @importFrom rlang .data
#' @importFrom IRanges reverse
#' @export
getUserAttractMotifs <-
    function(species  = "Hsapiens",
        reverse = FALSE,
        pathToMotifs = NULL) {
        options(readr.num_columns = 0)


        # Create a temporary directory
        td = tempdir()
        # Create the placeholder file
        tf = tempfile(tmpdir=td, fileext="ATtRACT.zip")
        # download into the placeholder file
        utils::download.file("https://attract.cnic.es/attract/static/ATtRACT.zip",
            tf)
        db <- suppressWarnings(read_tsv(unz(tf, "ATtRACT_db.txt")))


        # # Check if ATtRACT_db.txt is in the working directory.
        # # If not it will be downloaded from ATtRACT database
        # if (!file.exists("ATtRACT.zip")) {
        #     utils::download.file("https://attract.cnic.es/attract/static/ATtRACT.zip",
        #         "./ATtRACT.zip")
        # }
        #
        # db <-
        #     suppressWarnings(readr::read_tsv(unz("./ATtRACT.zip", "ATtRACT_db.txt")))

        # Reformat how the species name is reported in ATtRACT database so
        # that it can be compared with the species given in input.
        # E.g Mus_musculus in ATtRACT db becomes Mmusculus. the last one is
        # how it is reported in BSgenome.
        el1 <-
            substr(unlist(lapply(
                base::strsplit(db$Organism, "_"), "[", 1
            )), 1, 1)
        el2 <-
            tolower(unlist(lapply(
                base::strsplit(db$Organism, "_"), "[", 2
            )))
        db$Organism <- paste0(el1, el2)

        # if the organims given in input is not present in ATtRACT db
        # then take the human motifs.

        if (species %in% db$Organism) {
            attractRBPmotifs <- db %>%
                dplyr::filter(.data$Organism == species) %>%
                dplyr::select(.data$Gene_name,
                    .data$Motif) %>%
                dplyr::rename(id = .data$Gene_name,
                    motif = .data$Motif)
        } else{
            cat(
                paste(
                    "Organism not found in ATtRACT db, the human RBP
                    motifs in the ATtRACT database will be used.",
                    sep = " "
                )
            )
            attractRBPmotifs <- db %>%
                dplyr::filter(.data$Organism == "Hsapiens") %>%
                dplyr::select(.data$Gene_name,
                    .data$Motif) %>%
                dplyr::rename(id = .data$Gene_name,
                    motif = .data$Motif)
        }


        if (is.null(pathToMotifs)) {
            pathToMotifs <- "motifs.txt"
        }


        if (file.exists(pathToMotifs)) {
            # Read file containing user given patterns
            motifsFromFile <-
                utils::read.table(
                    pathToMotifs,
                    stringsAsFactors = FALSE,
                    header = TRUE,
                    sep = "\t"
                )
            # colnames(motifsFromFile) <-
            #     c("id", "motif", "length")


        } else{
            motifsFromFile <- data.frame(matrix(nrow = 0, ncol = 3))
            colnames(motifsFromFile) <- c("id", "motif", "length")
        }

        if (nrow(motifsFromFile) == 0) {
            cat("motifs.txt is empty or absent. Only
                ATtRACT motifs will be analyzed")
        }

        # If the file is empty then only the ATtRACT motifs are analyzed
        if (nrow(motifsFromFile) > 0 & reverse) {
            # reverse motifs given in input
            reversedMotifsFromFile <- motifsFromFile

            for (m in seq_along(reversedMotifsFromFile$motif)) {
                reversedMotifsFromFile$motif[m] <-
                    reversedMotifsFromFile$motif[m] %>%
                    gsub("\\[", "Z", .) %>%
                    gsub("]", "X", .) %>%
                    IRanges::reverse() %>%
                    gsub("X", "\\[", .) %>%
                    gsub("Z", "]", .)
            }

            newMotifsFromFile <-
                dplyr::bind_rows(motifsFromFile, reversedMotifsFromFile)

        } else{
            newMotifsFromFile <- motifsFromFile
        }


        if (reverse) {
            # we revere the motifs so that they can be analyzed also
            # in the other orientation
            reverseAttractRBPmotifs <- attractRBPmotifs
            reverseAttractRBPmotifs$motif <-
                IRanges::reverse(reverseAttractRBPmotifs$motif)

            newAttractRBPmotifs <-
                dplyr::bind_rows(attractRBPmotifs, reverseAttractRBPmotifs)
            newAttractRBPmotifs <-
                newAttractRBPmotifs[!duplicated(newAttractRBPmotifs),]
        } else{
            newAttractRBPmotifs <- attractRBPmotifs
        }

        userAttractMotifs <-
            dplyr::bind_rows(newMotifsFromFile[, c(1, 2)], newAttractRBPmotifs)


        return(userAttractMotifs)
    }



#' @title Filter motifs
#'
#' @description The function filterMotifs() finds motifs that match with motifs
#' of known RNA Binding Proteins (RBPs) deposited in the ATtRACT database and
#' with motifs specified by the user reported in motifs.txt. Subsequently,
#' they are filtered based on the value of rbp argument.
#'
#' @param computedMotifs A character vector containing motifs extracted from
#' the target sequences. It can be generated with \code{\link{computeMotifs}}.
#'
#' @param species A string specifying the species of the RBP motifs to use.
#' Type data(attractSpecies) to see the possible options.
#' Default values is "Hsapiens".
#'
#' @param rbp A logical specifying whether to report only motifs matching
#' with known RBP motifs from ATtRACT database or user motifs specified in
#' motifs.txt. If FALSE is specified only motifs that do not match with any of
#' these motifs are reported. Default value is TRUE.
#'
#' @param reverse A logical specifying whether to reverse the motifs collected
#' from ATtRACT database and from motifs.txt. If TRUE is specified all the
#' motifs are reversed and analyzed together with the direct motifs as they are
#' reported in the ATtRACT db and motifs.txt. Default value is FALSE.
#'
#' @param pathToMotifs A string containing the path to the motifs.txt
#' file. The file motifs.txt contains motifs/regular expressions specified
#' by the user. It must have 3 columns with headers:
#' - id (1st column): name of the motif. - e.g. RBM20 or motif1.
#' - motif (2nd column): motif/pattern to search.
#' - length (3rd column): length of the motif.
#' By default pathToMotifs is set to NULL and the file it is searched in the
#' working directory. If motifs.txt is located in a different directory then
#' the path needs to be specified. If this file is absent or empty only the
#' motifs of RNA Binding Proteins in the ATtRACT database are considered in
#' the motifs analysis.
#'
#' @return A data frame
#'
#' @examples
#' # Load data frame containing detected back-spliced junctions
#' data("mergedBSJunctions")
#'
#' # Load short version of the gencode v19 annotation file
#' data("gtf")
#'
#' # Example with the first back-spliced junctions.
#' # Multiple back-spliced junctions can also be analyzed at the same time.
#'
#' # Annotate detected back-spliced junctions
#' annotatedBSJs <- annotateBSJs(mergedBSJunctions[1, ], gtf)
#'
#' # Retrieve targets
#' targets <- getSeqsFromGRs(
#'     annotatedBSJs,
#'     lIntron = 200,
#'     lExon = 10,
#'     type = "ie",
#'     species = "Hsapiens",
#'     genome = "hg19" )
#'
#' # Inner function
#' # Compute motifs
#' computedMotifs <- computeMotifs(
#'     targets,
#'     width = 6)
#'
#' # Inner function
#' # Filter motifs
#' filteredMotifs <- filterMotifs(
#'    computedMotifs,
#'    species = "Hsapiens",
#'    rbp = TRUE,
#'    reverse = FALSE)
#'
#' @import dplyr
#' @import magrittr
#' @importFrom stringr str_extract
#' @export
filterMotifs <-
    function(computedMotifs,
        species  = "Hsapiens",
        rbp = TRUE,
        reverse = FALSE,
        pathToMotifs = NULL) {
        # Identify motifs matching with RBPs
        filteredMotifs <-
            data.frame(matrix(nrow = length(computedMotifs),
                ncol = 2))
        colnames(filteredMotifs) <-  c("motif", "id")
        filteredMotifs$motif <- computedMotifs

        widthCompMotifs <- nchar(computedMotifs[1])

        # Get user and ATtRACT RBP motifs
        userATtRACTmotifs <-
            getUserAttractMotifs(species,  reverse, pathToMotifs)

        . <- NULL # For R CMD check
        for (j in seq_along(filteredMotifs$motif)) {
            # Check whether the motifs matches with or it is contanined within
            # any RBP motifs

            # Grep do not work with pattern. In motifs.txt the user can reports
            # patterns.
            # By using grep there is a hit if there is a perfect match between
            # 2 strings or the first string it is contained as substring within
            # the second
            grepedM <-
                userATtRACTmotifs[base::grep(filteredMotifs$motif[j], userATtRACTmotifs$motif), ] %>%
                dplyr::mutate(motif = as.character(.data$motif),
                    id = as.character(.data$id))

            # str_extract works with pattern. In motifs.txt the user can reports
            # patterns.
            extractedM <-
                base::cbind(
                    stringr::str_extract(filteredMotifs$motif[j], userATtRACTmotifs$motif),
                    userATtRACTmotifs$id
                ) %>%
                magrittr::set_colnames(c("motif", "id")) %>%
                as.data.frame() %>%
                dplyr::select(.data$id, .data$motif) %>%
                dplyr::filter(!is.na(.data$motif)) %>%
                dplyr::mutate(motif = as.character(.data$motif),
                    id = as.character(.data$id)) %>%
                dplyr::filter(base::nchar(.data$motif) >= widthCompMotifs)

            joinedM <- dplyr::bind_rows(grepedM, extractedM) %>%
                dplyr::filter(!duplicated(.))

            if (nrow(joinedM) > 0) {
                filteredMotifs$id[j] <- paste(unique(joinedM$id), collapse = ",")
            }
        }

        # Filter
        if (rbp) {
            # Keep motifs matching with known RBP motifs
            filteredMotifs <- filteredMotifs %>%
                dplyr::filter(!is.na(.data$id))
        } else{
            # Keep unknown motifs
            filteredMotifs <- filteredMotifs %>%
                dplyr::filter(is.na(.data$id)) %>%
                dplyr::mutate(id = paste("m", seq(1, length(.data$id)), sep = ""))
        }

        return(filteredMotifs)
    }



#' @title Create list to store motif results
#'
#' @param targets A list containing the target sequences to analyze.
#' It can be generated with \code{\link{getCircSeqs}}.
#'
#' @param filteredMotifs A data frame containing motifs to analyze.
#' It can be generated with \code{\link{filterMotifs}}.
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
#' # Example with the first back-spliced junctions.
#' # Multiple back-spliced junctions can also be analyzed at the same time.
#'
#' # Annotate detected back-spliced junctions
#' annotatedBSJs <- annotateBSJs(mergedBSJunctions[1, ], gtf)
#'
#' # Retrieve targets
#' targets <-
#' getSeqsFromGRs(
#'     annotatedBSJs,
#'     lIntron = 200,
#'     lExon = 10,
#'     type = "ie",
#'     species = "Hsapiens",
#'     genome = "hg19")
#'
#' # Inner function
#' # Compute motifs
#' computedMotifs <- computeMotifs(
#'     targets,
#'     width = 6)
#'
#' # Inner function
#' # Filter motifs
#' filteredMotifs <- filterMotifs(
#'     computedMotifs,
#'     species = "Hsapiens",
#'     rbp = TRUE,
#'     reverse = FALSE)
#'
#' # Inner function
#' # Create list
#' motifs <- createMotifsList(
#'     targets,
#'     filteredMotifs)
#'
#'
#' @export
createMotifsList <-
    function(targets, filteredMotifs) {
        if (length(targets) == 2 & names(targets)[[1]] == "upGR") {
            # Create an empty list of 2 elements
            motifs <- vector("list", 2)
            names(motifs)[1] <- "upGR"
            names(motifs)[2] <- "downGR"

        } else if (length(targets) == 1 &
                names(targets)[[1]] == "bsj") {
            # Create a enmpty list of 1 elements
            motifs <- vector("list", 1)
            names(motifs)[1] <- "bsj"

        } else if (length(targets) == 1 &
                names(targets)[[1]] == "circ") {
            # Create a enmpty list of 1 elements
            motifs <- vector("list", 1)
            names(motifs)[1] <- "circ"
        }



        for (i in seq_along(motifs)) {
            # Put the string NA so that the we do not get error when calling
            # RNAStringSet
            targets[[i]]$seq[is.na(targets[[i]]$seq)] <- "NA"

            # Create an empty list of 4 elements to store the extracted
            # information
            motifs[[i]] <- vector("list", 4)
            names(motifs[[i]])[1] <- "targets"
            names(motifs[[i]])[2] <- "counts"
            names(motifs[[i]])[3] <- "locations"
            names(motifs[[i]])[4] <- "motifs"

            # Fill the dataframe with the target sequences
            motifs[[i]]$targets <- targets[[i]]

            # Fill the data frame with the found motifs
            motifs[[i]]$motifs <-  filteredMotifs

            if (nrow(filteredMotifs) > 0) {
                # Create the empty dataframe to store the occurences of motifs
                # found in the target sequences
                motifs[[i]]$counts <-
                    data.frame(matrix(
                        nrow = nrow(targets[[i]]),
                        ncol = length(filteredMotifs$motif) + 1
                    ))
                colnames(motifs[[i]]$counts) <-
                    c("id", c(filteredMotifs$motif))
                motifs[[i]]$counts$id <- targets[[i]]$id

                # Create the empty dataframe to store the location of motifs
                # found in the target sequences
                motifs[[i]]$locations <-
                    data.frame(matrix(
                        nrow = nrow(targets[[i]]),
                        ncol = length(filteredMotifs$motif) + 1
                    ))
                colnames(motifs[[i]]$locations) <-
                    c("id", c(filteredMotifs$motif))
                motifs[[i]]$locations$id <- targets[[i]]$id

            }
        }

        return(motifs)

    }



#' @title Group motifs shared by multiple RBPs
#'
#' @description A same RBP can recognize multiple motifs, the function
#' mergeMotifs() groups all the motifs found for each RBP and report the
#' total counts.
#'
#' @param motifs A data frame generated with \code{\link{getMotifs}}.
#'
#' @return A data frame.
#'
#' @examples
#' # Load data frame containing detected back-spliced junctions
#' data("mergedBSJunctions")
#'
#' # Load short version of the gencode v19 annotation file
#' data("gtf")
#'
#' # Example with the first back-spliced junctions.
#' # Multiple back-spliced junctions can also be analyzed at the same time.
#'
#' # Annotate detected back-spliced junctions
#' annotatedBSJs <- annotateBSJs(mergedBSJunctions[1, ], gtf)
#'
#' # Retrieve targets
#' targets <-
#' getSeqsFromGRs(
#'     annotatedBSJs,
#'     lIntron = 200,
#'     lExon = 10,
#'     type = "ie",
#'     species = "Hsapiens",
#'     genome = "hg19")
#'
#' # Get motifs
#' motifs <-
#' getMotifs(
#'     targets,
#'     width = 6,
#'     species = "Hsapiens",
#'     rbp = TRUE,
#'     reverse = FALSE)
#'
#' # Group motifs
#' mergedMotifs <- mergeMotifs(motifs)
#'
#' @importFrom rlang .data
#' @importFrom reshape2 melt
#' @import dplyr
#' @import magrittr
#' @export
mergeMotifs <- function(motifs) {
    if (nrow(motifs[[1]]$motifs) == 0) {
        # Make empty data frame
        mergedMotifsAll <- data.frame(matrix(nrow = 0, ncol = 3))
        colnames(mergedMotifsAll) <-  c("id", "count", "motif")
    } else{
        mergedMotifs <- list()
        for (i in seq_along(motifs)) {
            # Reshape the data frame
            mergedMotifs[[i]] <- motifs[[i]]$counts %>%
                reshape2::melt(
                    id.vars = c("id"),
                    variable.name = "motif",
                    value.name = "count"
                ) %>%
                # dplyr::mutate_all(funs(replace(., is.na(.), 0)))%>%
                dplyr::group_by(.data$motif) %>%
                dplyr::summarise(count = sum(.data$count, na.rm = TRUE)) %>%
                dplyr::mutate(motif = as.character(.data$motif))

        }

        # Check whether a same motif is shared by multiple RBPs.
        # If so retrive and duplicate those motifs by reporting one
        # RBP name.

        # In case the up and the down sequences have been analyzed the
        # motifs[[1]]$motif and motifs[[2]]$motif are the same,
        # so we use the motifs reported in motifs[[1]]$motif.
        toSplit <-
            motifs[[1]]$motifs[base::grep(",", motifs[[1]]$motif$id), ]

        if (nrow(toSplit) >= 1) {
            # Remove motifs shared by multiple RBPs
            cleanedMotifs <-
                motifs[[1]]$motif[-base::grep(",", motifs[[1]]$motif$id), ]

            # Ducplicate motifs
            rbpsWithSharedMotifs <-
                as.data.frame(matrix(nrow = 0, ncol = 2))
            colnames(rbpsWithSharedMotifs) <- c("id", "motif")
            for (j in seq_along(toSplit$motif)) {
                id <- base::strsplit(toSplit$id[j], ",")[[1]]
                for (b in seq_along(id)) {
                    rbpsWithSharedMotifs <- dplyr::bind_rows(
                        rbpsWithSharedMotifs,
                        as.data.frame(cbind(id[b], toSplit$motif[j])) %>%
                            magrittr::set_colnames(c("id", "motif")) %>%
                            dplyr::mutate(
                                id = as.character(.data$id),
                                motif = as.character(.data$motif)
                            )
                    )
                }
            }

            splittedRBPs <-
                dplyr::bind_rows(rbpsWithSharedMotifs, cleanedMotifs)
        } else{
            splittedRBPs <- motifs[[1]]$motif
        }


        . <- NULL # For R CMD check
        if (length(mergedMotifs) == 2) {
            mergedMotifsAll <-
                dplyr::bind_rows(mergedMotifs[[1]], mergedMotifs[[2]]) %>%
                dplyr::group_by(.data$motif) %>%
                dplyr::summarise(count = sum(.data$count)) %>%
                base::merge(
                    .,
                    splittedRBPs,
                    by = "motif",
                    all = TRUE,
                    sort = FALSE
                ) %>%
                dplyr::ungroup() %>%
                dplyr::group_by(.data$id) %>%
                dplyr::summarise(
                    count = sum(.data$count),
                    motif = paste(.data$motif, collapse = ",")
                ) %>%
                as.data.frame()


        } else {
            mergedMotifsAll <- mergedMotifs[[1]] %>%
                base::merge(
                    .,
                    splittedRBPs,
                    by = "motif",
                    all = TRUE,
                    sort = FALSE
                ) %>%
                dplyr::ungroup() %>%
                dplyr::group_by(.data$id) %>%
                dplyr::summarise(
                    count = sum(.data$count),
                    motif = paste(.data$motif, collapse = ",")
                ) %>%
                as.data.frame()
        }
    }


    return(mergedMotifsAll)
}
