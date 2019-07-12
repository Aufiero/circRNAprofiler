
# Functions in this file are used by multiple functions.


# The function getBasicColNames() returns the basic column names.
.getBasicColNames <- function() {
    basicColumns <- c("id",
        "gene",
        "strand",
        "chrom",
        "startUpBSE", # back-spliced junction
        "endDownBSE") # back-spliced junction
    return(basicColumns)

}


# The function checkBSJsDF() verifies that the functions to import
# the detected circRNAs generate the correct data structure and content.
.checkBSJsDF <- function(backSplicedJunctions, addColNames = NULL) {

    # Remove rows where strand is not available.
    # (Note: This step is done if the backSplicedJunctions data frame
    # contains liftedover coordinates. In this case if the strand is NA
    # means that the coordinates were not found in the new species:
    # conversion failed)
    backSplicedJunctions <- backSplicedJunctions[!is.na(backSplicedJunctions$strand),]

    colNames <- c(.getBasicColNames(), addColNames)

    if (!all(colNames  %in% colnames(backSplicedJunctions))) {
        missingNamesId <-
            which(!colnames(backSplicedJunctions) %in% colNames)
        stop("missing or wrong column names: ",
            paste(colnames(backSplicedJunctions)[missingNamesId],
                collapse = " \t"))
    }

    for (i in seq_along(backSplicedJunctions$strand)) {
        # For gene located on the positive strand the coordinates of an
        # acceptor site of an exon must be smaller than the coordinates of
        # the donor site
        if (backSplicedJunctions$strand[i] == "+" &
                backSplicedJunctions$startUpBSE[i] > backSplicedJunctions$endDownBSE[i]) {
            stop(
                "found startUpBSE > endDownBSE - for gene on the positive
                strand the coordinates of an startUpBSE (acceptor site)
                of an exon must be smaller than the coordinates
                of the endDownBSE(donor site)"
            )
        }

        # For gene located on the negative strand the coordinates of an
        # acceptor site of an exon must be greater than the coordinates
        # of the donor site
        if (backSplicedJunctions$strand[i] == "-" &
                backSplicedJunctions$startUpBSE[i] < backSplicedJunctions$endDownBSE[i]) {
            stop(
                "found startUpBSE < endDownBSE - for gene on the negative
                strand the coordinates of an startUpBSE (acceptor site)
                of an exon must be greater than the coordinates
                of the endDownBSE (donor site)"
            )
        }
        }
    return(backSplicedJunctions)
        }



# The function getTargetsColNames() returns the column names.
.getTargetsColNames <- function() {
    colNames <- c(
        "id",
        "gene",
        "transcript",
        "strand",
        "chrom",
        "startGR",
        "endGR",
        "length",
        "seq",
        "type"
    )
    return(colNames)

}

# create target data frame
.getTargetsDF <- function(annotatedBSJs){
    targets <- data.frame(matrix(
        nrow = nrow(annotatedBSJs),
        ncol = length(.getTargetsColNames())
    ))
    colnames(targets) <- .getTargetsColNames()
    return(targets)

}


# get exon sequences from BS genome
.getExonSeqs <- function(transcript, genome) {

    exonSeqs <- as.character()
    for (i in seq_along(transcript$exon_number)) {
        exonSeqs[i] <- gsub("T",
            "U",
            as.character(
                BSgenome::getSeq(
                    genome,
                    names = transcript$chrom[i],
                    transcript$start[i],
                    transcript$end[i]
                )
            ))
    }
    return(exonSeqs)
}


# Generate a unique identifier by combining the values of the following
# columns: gene, strand, chrom, endDownBSE and startUpBSE.
# The values are separated by a semicolumns (:)
.getID <- function(circTable) {
    # Generate a unique identifier by combining the values of the following
    # columns: gene, strand, chrom, endDownBSE and startUpBSE.
    # The values are separated by a semicolumns (:)
    id <- paste(
        circTable$gene,
        circTable$strand,
        circTable$chrom,
        circTable$startUpBSE,
        circTable$endDownBSE,
        sep = ":"
    )
    return(id)
}


# Fix the start and the end coordinates for the positive and negative
# strand.
# For a cicRNA arising from  gene located on the negative strand the
# coordinates of the upstream exon in the pre-mRNA have to be greater
# than the coordinates of downstream exon. The "for" statment swithces
# the coordinates.
.fixCoords <- function(circTable) {
    for (i in seq_along(circTable$id)) {
         if (!is.na(circTable$strand[i])  &
                !is.na(circTable$chrom[i]) &
                !is.na(circTable$startUpBSE[i]) &
                !is.na(circTable$endDownBSE[i]) &
                circTable$strand[i] == "-" &
                circTable$startUpBSE[i] <
                circTable$endDownBSE[i]) {
            x <- circTable$startUpBSE[i]
            circTable$startUpBSE[i] <-
                circTable$endDownBSE[i]
            circTable$endDownBSE[i] <- x

        } else if (!is.na(circTable$strand[i])  &
                !is.na(circTable$chrom[i]) &
                !is.na(circTable$startUpBSE[i]) &
                !is.na(circTable$endDownBSE[i]) &
                circTable$strand[i] == "+" &
                circTable$startUpBSE[i] >
                circTable$endDownBSE[i]) {
            x <- circTable$startUpBSE[i]
            circTable$startUpBSE[i] <-
                circTable$endDownBSE[i]
            circTable$endDownBSE[i] <- x
        }
    }
    return(circTable)
}

# Read experiments.txt
.readExperiment <- function(pathToExperiment = NULL) {
    if (is.null(pathToExperiment)) {
        pathToExperiment <- "experiment.txt"
    }
    if (file.exists(pathToExperiment)) {
        # Read from path given in input
        experiment <-
            utils::read.table(
                pathToExperiment,
                stringsAsFactors = FALSE,
                header = TRUE,
                sep = "\t"
            )
    } else{
        experiment <- data.frame()
    }
    return(experiment)
}

# Read traits.txt
.readTraits <- function(pathToTraits = NULL) {
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
    } else{
        traitsFromFile <- data.frame()
    }

    return(traitsFromFile)
}

# get motifs.txt
.readMotifs <- function(pathToMotifs = NULL) {
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

    } else{
        motifsFromFile <- data.frame(matrix(nrow = 0, ncol = 3))
        colnames(motifsFromFile) <- c("id", "motif", "length")
    }

    return(motifsFromFile)
}

# Read pathToTranscript.txt
.readTranscripts <-
    function(pathToTranscripts = NULL, isRandom = FALSE) {
        if (is.null(pathToTranscripts)) {
            pathToTranscripts <- "transcripts.txt"
        }

        # Read from working directory
        if (file.exists(pathToTranscripts) & !isRandom) {
            transcriptsFromFile <-
                utils::read.table(
                    pathToTranscripts,
                    stringsAsFactors = FALSE,
                    header = TRUE,
                    sep = "\t"
                )

        } else {
            transcriptsFromFile <- data.frame()
        }
        return(transcriptsFromFile)
    }

# Read miRs.txt
.readMiRs <- function(pathToMiRs = NULL){
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
    return(miRsFromFile)
}

# Clean targets dataframe by removing the rows with NA values to avoid
# getting errors with the makeGRangesFromDataFrame function that does
# not handle NA values
.cleanTargets <- function(targets) {
    index1 <- which(
        !is.na(targets[[1]]$transcript) &
            !is.na(targets[[1]]$startGR) &
            !is.na(targets[[1]]$endGR)
    )
    index2 <- which(
        !is.na(targets[[2]]$transcript) &
            !is.na(targets[[2]]$startGR) &
            !is.na(targets[[2]]$endGR)
    )

    targets[[1]] <- targets[[1]][intersect(index1, index2), ]
    targets[[2]] <- targets[[2]][intersect(index1, index2), ]

    return(targets)
}
