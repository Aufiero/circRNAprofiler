#' @title Screen target sequences for recurrent motifs
#'
#' @description The function getMotifs() scans the target sequences for the
#' presence of recurrent motifs of a specific length defined in input.
#' By setting rbp equals to TRUE, the identified motifs are matched with motifs
#' of known RNA Binding Proteins (RBPs) deposited in the ATtRACT
#' (\url{http://attract.cnic.es}) or MEME database (\url{http://meme-suite.org/})
#' and with motifs specified by the user.
#' The user motifs must go in the file motifs.txt. If this file is absent or
#' empty, only motifs from the ATtRACT or MEME database are considered in the analysis.
#' By setting rbp equals to FALSE, only motifs that do not match with any motifs
#' deposited in the databases or user motifs are reported in the final
#' output. Location of the selected motifs is also reported. This corresponds
#' to the start position of the motif within the sequence (1-index based).
#'
#' @param targets A list containing the target sequences to analyze.
#' It can be generated with \code{\link{getCircSeqs}},
#' \code{\link{getSeqsAcrossBSJs}} or \code{\link{getSeqsFromGRs}}.
#'
#' @param width An integer specifying the length of all possible motifs to
#' extract from the target sequences. Default value is 6.
#'
#' @param database A string specifying the RBP database to use. Possible options
#' are ATtRACT or MEME. Default database is "ATtRACT".
#'
#' @param species A string specifying the species of the ATtRACT RBP motifs to
#' use. Type data(attractSpecies) to see the possible options.
#' Default value is "Hsapiens".
#'
#' @param memeIndexFilePath An integer specifying the index of the file path
#' of the meme file to use.Type data(memeDB) to see the possible options.
#' Default value is 18 corresponding to the following file:
#' motif_databases/RNA/Ray2013_rbp_Homo_sapiens.meme
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
#' \describe{
#' \item{id:}{ (1st column) - name of the motif. - e.g. RBM20 or motif1).}
#' \item{motif:}{(2nd column) -motif/pattern to search.}
#' \item{length:}{(3rd column) - length of the motif.}
#' }
#'
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
#' # Example with the first back-spliced junction
#' # Multiple back-spliced junctions can also be analyzed at the same time
#'
#' # Annotate the first back-spliced junction
#' annotatedBSJs <- annotateBSJs(mergedBSJunctions[1, ], gtf)
#'
#' # Get genome
#' if (requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)){
#' 
#' genome <- BSgenome::getBSgenome("BSgenome.Hsapiens.UCSC.hg19")
#'
#' # Retrieve target sequences
#' targets <- getSeqsFromGRs(
#'     annotatedBSJs,
#'     genome,
#'     lIntron = 200,
#'     lExon = 10,
#'     type = "ie"
#'     )
#'
#' # Get motifs
#' motifs <- getMotifs(
#'     targets,
#'     width = 6,
#'     database = 'ATtRACT',
#'     species = "Hsapiens",
#'     rbp = TRUE,
#'     reverse = FALSE)
#' 
#' }
#' 
#'
#' @importFrom readr read_tsv
#' @importFrom utils download.file
#' @importFrom utils untar
#' @importFrom rlang .data
#' @importFrom IRanges reverse
#' @importFrom stringr str_count
#' @importFrom stringr str_extract
#' @importFrom stringi stri_locate_all
#' @importFrom Biostrings RNAStringSet
#' @importFrom Biostrings oligonucleotideFrequency
#' @importFrom universalmotif read_meme
#' @import dplyr
#' @import magrittr
#' @export
getMotifs <-
    function(targets,
             width = 6,
             database = 'ATtRACT',
             species  = "Hsapiens",
             memeIndexFilePath = 18,
             rbp = TRUE,
             reverse = FALSE,
             pathToMotifs = NULL) {
        if (width < 4 | width > 20) {
            stop("width must be an integer between 4 and 20")
        }
        # Compute all motifs of a given length (width)
        computedMotifs <- .computeMotifs(targets, width)
        # Filter motifs matching with RBPs
        filteredMotifs <-
            .filterMotifs(computedMotifs,
                          database,
                          species,
                          memeIndexFilePath,
                          rbp,
                          reverse,
                          pathToMotifs)
        # Create list to store motif reults
        motifs <- .createMotifsList(targets, filteredMotifs)
        for (i in seq_along(motifs)) {
            targetsToAnalyze <- motifs[[i]]$targets
            # Get sequences in RNAStringSet format
            rnaSS <- .getRNAss(targetsToAnalyze, width)

            # Count occurence and location of the filtered motifs
            for (j in seq_along(motifs[[i]]$motif$motif)) {
                stringMotif <- motifs[[i]]$motif$motif[j]
                # Find location. Consider also overlapping patterns (?=pattern)
                locations <- stringi::stri_locate_all(as.character(rnaSS),
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
                    stringr::str_count(as.character(rnaSS),
                                       paste("(?=", stringMotif, ")", sep = ""))
            }
        }
        return(motifs)
    }


# The function computeMotifs() first computes all possible motifs
.computeMotifs <-
    function(targets, width = 6) {
        # Compute all motifs of length width
        computedMotifs <- as.character()
        for (i in seq_along(targets)) {
            targetsToAnalyze <- targets[[i]]
            # Get sequences in RNAStringSet format
            rnaSS <- .getRNAss(targetsToAnalyze, width)

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





# retrieve target sequence to analyze in RNAStringSet format
.getRNAss <- function(targetsToAnalyze, width = 6) {
    # Put NA to avoid error when calling RNAStringSet
    targetsToAnalyze$seq[is.na(targetsToAnalyze$seq)] <- "NA"

    # Adjust sequences
    if (targetsToAnalyze$type[1] == "circ") {
        # Adjust the length of circ seqs
        ss <- base::substring(targetsToAnalyze$seq, 30,
                              ((targetsToAnalyze$length + 30) + (width - 1)))

        # Put NA to avoid error when calling RNAStringSet
        ss[is.na(ss)] <- "NA"
        rnaSS <- Biostrings::RNAStringSet(ss)

    } else if (targetsToAnalyze$type[1] == "bsj") {
        # We adjust the length of BSJ seqs. Motifs must have at least
        # 1 nucleotide crossing the BSJ.
        r1 <-
            (base::nchar(targetsToAnalyze$seq) / 2) - (width - 2)
        r2 <-
            (base::nchar(targetsToAnalyze$seq) / 2) + (width - 1)
        rnaSS <-
            Biostrings::RNAStringSet(base::substring(targetsToAnalyze$seq, r1, r2))

    } else{
        rnaSS <- Biostrings::RNAStringSet(targetsToAnalyze$seq)
    }
    return(rnaSS)
}


# The function getUserDBmotifs() reads the user motifs
# in motifs.txt and retrieves the motifs deposited in the ATtRACT database
# (\url{http://attract.cnic.es}).
.getUserDBmotifs <-
    function(database = 'ATtRACT',
             species  = "Hsapiens",
             memeIndexFilePath = 18,
             reverse = FALSE,
             pathToMotifs = NULL) {
        options(readr.num_columns = 0)


        if (database == 'ATtRACT') {
            # Get motifs from attract data base (it contains motifs for 159 human RBPs)
            rbpMotifsFromDB <- .getRBPmotifsAttract(species)
        } else if (database == 'MEME') {
            # Get MEME motifs (it contains motifs for 80 human RBPs)
            rbpMotifsFromDB <- .getRBPmotifsMEME(memeIndexFilePath)
        } else{
            stop("database not correct, only ATtRACT or MEME are allowed")
        }

        # Read motifs.txt
        motifsFromFile <- .readMotifs(pathToMotifs)

        if (nrow(motifsFromFile) == 0) {
            cat("motifs.txt is empty or absent. Only",
                database,
                "motifs will be analyzed if available")
        }

        # we reverse the motifs so that they can be analyzed also
        # in the other orientation
        if (reverse) {
            # If the file is empty then only the ATtRACT or MEME motifs are analyzed
            if (nrow(motifsFromFile) > 0) {
                # reverse motifs given in input
                motifsFromFileNew <-
                    getReverseMotifs(motifsFromFile)
            }

            if (database == 'ATtRACT' & nrow(rbpMotifsFromDB) > 0) {
                rbpMotifsFromDBnew <-
                    .getReverseAttractRBPmotifs(rbpMotifsFromDB)
            } else if (nrow(rbpMotifsFromDB) > 0) {
                rbpMotifsFromDBnew <- getReverseMotifs(rbpMotifsFromDB)
            }


        } else{
            motifsFromFileNew <- motifsFromFile
            rbpMotifsFromDBnew <- rbpMotifsFromDB
        }

        userDBmotifs <-
            rbind(motifsFromFileNew[, c(1, 2)], rbpMotifsFromDBnew)

        return(userDBmotifs)
    }


# The function filterMotifs() finds motifs that match with motifs
# of known RNA Binding Proteins (RBPs) deposited in the ATtRACT database and
# with motifs specified by the user reported in motifs.txt. Subsequently,
# they are filtered based on the value of rbp argument.
.filterMotifs <-
    function(computedMotifs,
             database = 'ATtRACT',
             species  = "Hsapiens",
             memeIndexFilePath = 18,
             rbp = TRUE,
             reverse = FALSE,
             pathToMotifs = NULL) {
        # Identify motifs matching with RBPs
        filteredMotifs <- .createFilteredMotifsDF(computedMotifs)

        # Get user and ATtRACT RBP motifs
        userDBmotifs <-
            .getUserDBmotifs(database,
                             species,
                             memeIndexFilePath,
                             reverse,
                             pathToMotifs)

        # Check whether the motifs matches with or it is contanined within
        # any RBP motifs
        if(nrow(userDBmotifs)>0){
            filteredMotifs <-
                .matchWithKnowRBPs(filteredMotifs, userDBmotifs, computedMotifs)
        }

        # Filter
        if (rbp) {
            # Keep motifs matching with known RBP motifs
            filteredMotifs <- filteredMotifs %>%
                dplyr::filter(!is.na(id))
        } else{
            # Keep unknown motifs
            filteredMotifs <- filteredMotifs %>%
                dplyr::filter(is.na(id)) %>%
                dplyr::mutate(id = paste("m", seq(1, length(id)), sep = ""))
        }

        return(filteredMotifs)
    }

# Create filteredMotifs data frame
.createFilteredMotifsDF <- function(computedMotifs) {
    filteredMotifs <-
        data.frame(matrix(nrow = length(computedMotifs),
                          ncol = 2))
    colnames(filteredMotifs) <-  c("motif", "id")
    filteredMotifs$motif <- computedMotifs

    return(filteredMotifs)

}

# Check whether the motifs matches with or it is contanined within
# any RBP motifs
.matchWithKnowRBPs <-
    function(filteredMotifs,
             userDBmotifs,
             computedMotifs) {
        widthCompMotifs <- nchar(computedMotifs[1])
        for (j in seq_along(filteredMotifs$motif)) {
            # Grep do not work with pattern. In motifs.txt the user can reports
            # patterns.
            # By using grep there is a hit if there is a perfect match between
            # 2 strings or the first string it is contained as substring within
            # the second
            grepedM <-
                userDBmotifs[base::grep(filteredMotifs$motif[j], userDBmotifs$motif), ] %>%
                dplyr::mutate(motif = as.character(motif),
                              id = as.character(id))

            # str_extract works with pattern.
            extractedM <-
                base::cbind(
                    stringr::str_extract(filteredMotifs$motif[j], userDBmotifs$motif),
                    userDBmotifs$id
                ) %>%
                magrittr::set_colnames(c("motif", "id")) %>%
                as.data.frame() %>%
                dplyr::select(id, motif) %>%
                dplyr::filter(!is.na(motif)) %>%
                dplyr::mutate(motif = as.character(motif),
                              id = as.character(id)) %>%
                dplyr::filter(base::nchar(motif) >= widthCompMotifs)

            joinedM <- dplyr::bind_rows(grepedM, extractedM) %>%
                dplyr::filter(!duplicated(.))

            if (nrow(joinedM) > 0) {
                filteredMotifs$id[j] <- paste(unique(joinedM$id), collapse = ",")
            }
        }
        return(filteredMotifs)
    }


# Create list to store motif results
.createMotifsList <-
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
#' # Get genome
#' genome <- BSgenome::getBSgenome("BSgenome.Hsapiens.UCSC.hg19")
#'
#' # Retrieve target sequences
#' targets <- getSeqsFromGRs(
#'     annotatedBSJs,
#'     genome,
#'     lIntron = 200,
#'     lExon = 10,
#'     type = "ie"
#'     )
#'
#' # Get motifs
#' motifs <-
#' getMotifs(
#'     targets,
#'     width = 6,
#'     species = "Hsapiens",
#'     rbp = TRUE,
#'    reverse = FALSE)
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
            counts <- motifs[[i]]$counts
            mergedMotifs[[i]] <- .reshapeCounts(counts)
        }

        # For the up and the down sequences motifs[[1]]$motifs and
        # motifs[[2]]$motifs are the same, so we use the motifs reported
        # in motifs[[1]]$motif.
        motifsToSplit <- motifs[[1]]$motifs

        # Check whether a same motif is shared by multiple RBPs.
        # If so retrive and duplicate those motifs by reporting one
        # RBP name.
        splittedRBPs <- .splitRBPs(motifsToSplit)

        if (length(mergedMotifs) == 2) {
            mergedMotifsAll <-
                dplyr::bind_rows(mergedMotifs[[1]], mergedMotifs[[2]]) %>%
                dplyr::group_by(motif) %>%
                dplyr::summarise(count = sum(count)) %>%
                base::merge(
                    .,
                    splittedRBPs,
                    by = "motif",
                    all = TRUE,
                    sort = FALSE
                ) %>%
                dplyr::ungroup() %>%
                dplyr::group_by(id) %>%
                dplyr::summarise(
                    count = sum(count),
                    motif = paste(motif, collapse = ",")
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
                dplyr::group_by(id) %>%
                dplyr::summarise(
                    count = sum(count),
                    motif = paste(motif, collapse = ",")
                ) %>%
                as.data.frame()
        }
    }
    return(mergedMotifsAll)
}


# Reshape the data frame
.reshapeCounts <- function(counts) {
    reshapedCounts <- counts %>%
        reshape2::melt(
            id.vars = c("id"),
            variable.name = "motif",
            value.name = "count"
        ) %>%
        # dplyr::mutate_all(funs(replace(., is.na(.), 0)))%>%
        dplyr::group_by(motif) %>%
        dplyr::summarise(count = sum(count, na.rm = TRUE)) %>%
        dplyr::mutate(motif = as.character(motif))
    return(reshapedCounts)
}

# Check whether a same motif is shared by multiple RBPs.
# If so retrive and duplicate those motifs by reporting one
# RBP name.
# For the up and the down sequences motifs[[1]]$motifs and
# motifs[[2]]$motifs are the same, so we use the motifs reported
# in motifs[[1]]$motif.
.splitRBPs <- function(motifsToSplit) {
    toSplit <-
        motifsToSplit[base::grep(",", motifsToSplit$id), ]

    if (nrow(toSplit) >= 1) {
        # Remove motifs shared by multiple RBPs
        cleanedMotifs <-
            motifsToSplit[-base::grep(",", motifsToSplit$id), ]

        # Ducplicate motifs
        rbpsWithSharedMotifs <-
            as.data.frame(matrix(nrow = 0, ncol = 2))
        colnames(rbpsWithSharedMotifs) <- c("id", "motif")
        for (j in seq_along(toSplit$motif)) {
            id <- base::strsplit(toSplit$id[j], ",")[[1]]
            for (b in seq_along(id)) {
                rbpsWithSharedMotifs <- rbind(
                    rbpsWithSharedMotifs,
                    as.data.frame(cbind(id[b], toSplit$motif[j])) %>%
                        magrittr::set_colnames(c("id", "motif")) %>%
                        dplyr::mutate(
                            id = as.character(id),
                            motif = as.character(motif)
                        )
                )
            }
        }

        splittedRBPs <-
            dplyr::bind_rows(rbpsWithSharedMotifs, cleanedMotifs)
    } else{
        splittedRBPs <- motifsToSplit

    }

    return(splittedRBPs)
}





# Get RBP motifs from attract data base.
# Download the https://attract.cnic.es/attract/static/ATtRACT.zip.
# Motifs from the following file were used: ATtRACT_db.txt
.getRBPmotifsAttract <- function(species) {
    # Create a temporary directory
    td = tempdir()
    # Create the placeholder file
    tf = tempfile(tmpdir = td, fileext = "ATtRACT.zip")
    # download into the placeholder file

    url <- "https://attract.cnic.es/attract/static/ATtRACT.zip"
    tc <-
        tryCatch(
            utils::download.file(url, tf),
            warning = function(w)
                NULL
        )

    if (!is.null(tc)) {
        db <- suppressWarnings(read_tsv(unz(tf, "ATtRACT_db.txt"), show_col_types=FALSE))

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
                dplyr::filter(Organism == species) %>%
                dplyr::select(Gene_name,
                              Motif) %>%
                dplyr::rename(id = Gene_name,
                              motif = Motif)
        } else{
            cat(
                paste(
                    "Organism not found in ATtRACT db, the human RBP
                    motifs in the ATtRACT database will be analyzed.",
                    sep = " "
                )
            )
            attractRBPmotifs <- db %>%
                dplyr::filter(Organism == "Hsapiens") %>%
                dplyr::select(Gene_name,
                              Motif) %>%
                dplyr::rename(id = Gene_name,
                              motif = Motif)
        }
    } else{
        attractRBPmotifs <- data.frame(matrix(nrow = 0, ncol = 2))
        colnames(attractRBPmotifs) <-  c("id", "motif")
        cat('URL can not be reached: ',
            url,
            '.\nATtRACT motif can not be analyzed.')
    }

    return(attractRBPmotifs)
}



# Rerieve MEME consensus sequences and convert it to a regular expression.
# Download the http://meme-suite.org/meme-software/Databases/motifs/motif_databases.12.19.tgz.
# Motifs from the following file are used: RNA\Ray2013_rbp_All_Species.meme
.getRBPmotifsMEME <-
    function(memeIndexFilePath = 5,
             isDNA = FALSE) {
        # Create a temporary directory
        td = tempdir()
        # Create the placeholder file
        tf = tempfile(tmpdir = td, fileext = 'db')
        # download into the placeholder file
        url <-
            "http://meme-suite.org/meme-software/Databases/motifs/motif_databases.12.19.tgz"
        tc <-
            tryCatch(
                utils::download.file(url, tf),
                warning = function(w)
                    NULL
            )

        if (!is.null(tc)) {
            # Get the name of the first file in the zip archive
            memeDB <- utils::untar(tf, list = TRUE)
            memeFile <- grep('rbp', memeDB, value = TRUE) %>%
                data.frame() %>%
                magrittr::set_colnames('path') %>%
                dplyr::filter(., !grepl("dna_encoded",path)) %>%
                dplyr::mutate(index = seq_along(.data$path)) %>%
                dplyr::filter(.data$index == memeIndexFilePath)


            if (nrow(memeFile)) {
                # unzip the file to the temporary directory
                utils::untar(tf, files = memeFile$path, exdir = td)
                # fpath is the full path to the extracted file
                fpath = file.path(td, memeFile$path)

                # Read meme motifs from file
                memeMotifs <- .readMemeMotifs(fpath)
            } else{
                cat(
                    'Index not found, MEME motifs can not be analyzed.\nType(memeDB) to see all possible options.'
                )
            }


        } else{
            memeMotifs <- data.frame(matrix(nrow = 0, ncol = 3))
            colnames(memeMotifs) <-  c("id", "motif", "length")
            cat('URL can not be reached: ',
                url,
                '.\nMEme motifs can not be analyzed.')
        }

        return(memeMotifs)
    }

# Read meme motifs from file
.readMemeMotifs <- function(fpath) {
    meme <-
        universalmotif::read_meme(
            fpath,
            skip = 0,
            readsites = FALSE,
            readsites.meta = FALSE
        )
    # Create empty data frame
    memeMotifs <- data.frame(matrix(nrow = length(meme), ncol = 3))
    colnames(memeMotifs) <-  c("id", "motif",'length')

    # Convert the sequence to UPPER CASE

    for (i in seq_along(meme)) {
        memeMotifs$id[i] <- meme[[i]]@altname
        memeMotifs$length[i] <- nchar(meme[[i]]@consensus)
        memeMotifs$motif[i] <-
            getRegexPattern(meme[[i]]@consensus, isDNA = FALSE)
    }


    memeMotifs$motif <-  gsub('T', 'U', memeMotifs$motif)

    return(memeMotifs)
}




# get reverse motifs from file
getReverseMotifs <- function(motifsWithRegExp) {
    # reverse motifs given in input
    reversedMotifs <- motifsWithRegExp

    for (m in seq_along(reversedMotifs$motif)) {
        reversedMotifs$motif[m] <-
            reversedMotifs$motif[m] %>%
            gsub("\\[", "Z", .) %>%
            gsub("]", "X", .) %>%
            IRanges::reverse() %>%
            gsub("X", "\\[", .) %>%
            gsub("Z", "]", .)
    }

    motifsNew <-
        dplyr::bind_rows(motifsWithRegExp, reversedMotifs)

    return(motifsNew)
}


# Reverse motifs from attract data base
.getReverseAttractRBPmotifs <- function(rbpMotifsFromDB) {
    reverseAttractRBPmotifs <- rbpMotifsFromDB
    reverseAttractRBPmotifs$motif <-
        IRanges::reverse(reverseAttractRBPmotifs$motif)

    rbpMotifsFromDBnew <-
        dplyr::bind_rows(rbpMotifsFromDB, reverseAttractRBPmotifs)
    rbpMotifsFromDBnew <-
        rbpMotifsFromDBnew[!duplicated(rbpMotifsFromDBnew),]

    return(rbpMotifsFromDBnew)
}


# If the function you are looking for is not here check supportFunction.R
# Functions in supportFunction.R are used by multiple functions.
