#' @title Check project folder
#'
#' @description The function checkProjectFolder() verifies that the
#' project folder is set up correctly. It checks that the mandatory files
#' (.gtf file, the folders with the circRNAs_X.txt files and experiemnt.txt)
#' are present in the working directory.The function
#' \code{\link{initCircRNAprofiler}} can be used to initialize the project folder.
#'
#' @param pathToExperiment A string containing the path to the experiment.txt
#' file. The file experiment.txt contains the experiment design information.
#' It must have at least 3 columns with headers:
#' - label (1st column): unique names of the samples (short but informative).
#' - fileName (2nd column): name of the input files - e.g. circRNAs_X.txt, where
#' x can be can be 001, 002 etc.
#' - group (3rd column): biological conditions - e.g. A or B; healthy or
#' diseased if you have only 2 conditions.
#' By default pathToExperiment is set to NULL and the file it is searched in the
#' working directory. If experiment.txt is located in a different directory
#' then the path needs to be specified.
#'
#' @param pathToGTF A string containing the path to the the GTF file.
#' Use the same annotation file used during the RNA-seq mapping procedure.
#' By default pathToGTF is set to NULL and the file it is searched in the
#' working directory. If .gtf is located in a different directory then the
#' path needs to be specified.
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
#' motifs of RNA Binding Proteins in the ATtRACT or MEME database are considered
#' in the motifs analysis.
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
#' empty, all miRs of the species specified in input are considered in the
#' miRNA analysis.
#'
#' @param pathToTranscripts A string containing the path to the transcripts.txt
#' file. The file transcripts.txt contains the transcript ids of the
#' circRNA host gene to analyze. It must have one column with header id.
#' By default pathToTranscripts is set to NULL and the file it is searched in
#' the working directory. If transcripts.txt is located in a different
#' directory then the path needs to be specified. If this file is empty or
#' absent the longest transcript of the circRNA host gene containing the
#' back-spliced junctions are considered in the annotation analysis.
#'
#' @param pathToTraits A string containing the path to the traits.txt
#' file. contains diseases/traits specified by the user. It must
#' have one column with header id. By default pathToTraits is set to NULL and
#' the file it is searched in the working directory. If traits.txt is located
#' in a different directory then the path needs to be specified. If this file is
#' absent or empty SNPs associated with all diseases/traits in
#' the GWAS catalog are considered in the SNPs analysis.
#'
#' @return An integer. If equals to 0 the project folder is correctly
#' set up.
#'
#' @examples
#' checkProjectFolder()
#'
#' @importFrom utils read.table
#' @export
checkProjectFolder <-
    function(pathToExperiment = NULL,
        pathToGTF = NULL,
        pathToMotifs = NULL,
        pathToMiRs = NULL,
        pathToTranscripts = NULL,
        pathToTraits = NULL) {
        # Check optional files
        .checkMotifs(pathToMotifs)
        .checkTraits(pathToTraits)
        .checkMiRs(pathToMiRs)
        .checkTranscripts(pathToTranscripts)

        # Check mandatory files
        # Check GTF
        check1 <- .checkGTF(pathToGTF)
        # check experiment.txt and prediction results
        check2 <- .checkExperiment(pathToExperiment)

        checks <- check1 + check2
        return(checks)
    }



# Check motifs.txt
.checkMotifs <- function(pathToMotifs = NULL) {
    # Check  optional files
    # check motifs.txt
    motifsFromFile <- .readMotifs(pathToMotifs)

    if (nrow(motifsFromFile) > 0) {
        cnm <- c("id", "motif", "length")

        if (!all(cnm %in% colnames(motifsFromFile))) {
            missingNamesId <- which(!cnm %in%
                    colnames(motifsFromFile))
            cat(
                "(!) missing or wrong column names in motifs.txt: ",
                paste(cnm[missingNamesId], collapse = " \t"),
                "\n"
            )

        } else if (ncol(motifsFromFile) != 3) {
            cat("(!) motifs.txt must have 3 column with header id, motif and length\n")
        }

    } else{
        cat(
            "Missing or empty motifs.txt file.
            Optional file. If absent or empty only
            ATtRACT motifs will be analyzed\n"
        )
    }
    }


# check traits.txt
.checkTraits <- function(pathToTraits = NULL) {
    # Read traits.txt
    traitsFromFile <- .readTraits(pathToTraits)
    # Check if there there are traits
    if (nrow(traitsFromFile) > 0) {
        # Check if column id
        if (!"id" %in% colnames(traitsFromFile)) {
            cat("(!) missing or wrong column name in traits.txt: id\n ")
        } else if (ncol(traitsFromFile) != 1) {
            cat("(!) traits.txt must have 1 column with header id\n ")
        }


    } else {
        cat(
            "Missing or empty traits.txt file.
            Optional file. If absent or empty all
            traits in the GWAS catalog will be analyzed\n"
        )
    }
    }

# check miRs.txt
.checkMiRs <- function(pathToMiRs = NULL) {
    # Read miRs.txt
    miRsFromFile <- .readMiRs(pathToMiRs)

    if (nrow(miRsFromFile) > 0) {
        # Check if column id
        if (!"id" %in% colnames(miRsFromFile)) {
            cat("(!) missing or wrong column name in traits.txt: id\n ")
        } else if (ncol(miRsFromFile) != 1) {
            cat("(!) miRs.txt must have 1 column with header id\n ")
        }

    } else{
        cat(
            "Missing or empty miRs.txt file.
            Optional file. If absent or empty all miRNAs of the
            specified species will be analyzed\n"
        )
    }
    }

# check transcripts.txt
.checkTranscripts <- function(pathToTranscripts = NULL) {
    # check transcripts.txt
    transcriptsFromFile <- .readTranscripts(pathToTranscripts)

    if (nrow(transcriptsFromFile) > 0) {
        # Check if column id
        if (!"id" %in% colnames(transcriptsFromFile)) {
            cat("(!) missing or wrong column name in traits.txt: id\n ")
        } else if (ncol(transcriptsFromFile) != 1) {
            cat("(!) transcripts.txt must have 1 column with header id\n ")
        }

    } else{
        cat(
            "Missing or empty transcripts.txt.
            Optional file. If absent or empty the longest
            transcripts for all circRNAs will be analyzed\n"
        )
    }
    }


# Check GTF file
.checkGTF <- function(pathToGTF = NULL) {
    fileNames <- list.files()
    check <- 0
    # check GTF file
    if (is.null(pathToGTF)) {
        pathToGTF <- grep("gtf", fileNames, value = TRUE)[1]
    }

    if (is.na(pathToGTF)) {
        cat("(!): missing gtf file\n")
        check <- check + 1

    }

    return(check)
}

# check experiment.txt and prediction results
.checkExperiment <- function(pathToExperiment = NULL) {
    fileNames <- list.files()
    check <- 0
    # Read experiment.txt
    experiment <- .readExperiment(pathToExperiment)
    if (nrow(experiment) > 0) {
        cne <- c("label", "fileName", "condition")
        if (!all(cne %in%  colnames(experiment))) {
            missingNamesId <- which(!cne %in% colnames(experiment))
            cat(
                "(!): missing or wrong column names in experiment.txt: ",
                paste(cne[missingNamesId], collapse = " \t", "\n")
            )
            check <- check + 1
        }else if (ncol(experiment) != 3) {
            cat("(!) experiment must have 3 column with header label, fileName and condition\n")
        }

        # check folders with circRNA predictions
        predictionToolsAll <- getDetectionTools()
        if (sum(predictionToolsAll$name  %in% fileNames) >= 1) {
            pt <-
                predictionToolsAll$name[which(predictionToolsAll$name %in% fileNames)]

            for (i in seq_along(pt)) {
                if (!all(experiment$fileName %in% list.files(pt[i]))) {
                    missingFilesId <- which(!experiment$fileName %in% list.files(pt[i]))
                    cat(
                        "(!): .txt file reported in experiment.txt is not
                        present in folder named",
                        pt[i],
                        "\n"
                    )
                    cat(
                        "Missing files:",
                        paste(experiment$fileName[missingFilesId],
                            collapse = " \t"),
                        "\n"
                    )
                    check <- check + 1
                }
            }

        } else {
            cat("(!): missing folders containing circRNA predictions\n")
            cat(
                "Folders containing .txt files with circRNA predictions
                must be present in the wd\n"
            )
            check <- check + 1
        }
    } else {
        cat("(!): experiment.txt is absent or empty\n")
        check <- check + 1
    }
    return(check)
    }


# If the function you are looking for is not here check supportFunction.R
# Functions in supportFunction.R are used by multiple functions.
