#' @title Initialize the project folder
#'
#' @description The function initCircRNAprofiler() initializes the project
#' forlder.
#'
#' @param projectFolderName A character string specifying the name of the
#' project folder.
#'
#' @param detectionTools A character vector specifying the tools
#' used to predict circRNAs. The following options are allowed: mapsplice,
#' nclscan, knife, circexplorer2, circmarker and uroborus. If the tool is
#' not mapsplice, nclscan, knife, circexplorer2, uroborus or circmarker
#' then use the option other. The user can choose 1 or multiple tools.
#' Subfolders named as the specified tools will be generated under
#' the working directory.
#'
#' @return A NULL object
#'
#' @examples
#' initCircRNAprofiler(projectFolderName = "circProject",
#'     detectionTools = "mapsplice")
#'
#'
#' @importFrom magrittr %>%
#' @importFrom utils write.table
#' @export
initCircRNAprofiler <-
    function(projectFolderName,
        detectionTools) {
        tools <- getDetectionTools()
        if (!all(detectionTools %in% tools$name)) {
            stop(paste(
                detectionTools[detectionTools %in% tools$name == FALSE],
                "not supported. Supported prediction tools:",
                paste(tools$name, collapse = ", ")
            ))
        }
        if (file.exists(projectFolderName)) {
            stop(
                paste(
                    projectFolderName,
                    "already exists, rename existing project folder or
                    choose another name",
                    sep = " "
                )
            )
        } else {
            dir.create(projectFolderName)
        }

        # Create a folder for each prediction tool
        for (i in seq_along(detectionTools)) {
            dir.create(paste(projectFolderName, detectionTools[i], sep = "/"))
        }

        # Create a template for the experiment.txt
        getExperimentTamplate(projectFolderName)

        # Create a tamplate for the motifs.txt
        getMotifsTamplate(projectFolderName)

        # Create a tamplate for the traits.txt
        getTraitsTamplate(projectFolderName)


        # Create a tamplate for the miRs.txt
        getMiRsTamplate(projectFolderName)

        # Create a tamplate for the transcripts.txt
        getTranscriptsTamplate(projectFolderName)
    }

# Experiment template
getExperimentTamplate <- function(projectFolderName) {
    experiment <-
        data.frame(matrix(nrow = 0, ncol = 3))
    colnames(experiment) <-
        c("label", "fileName", "condition")
    file.create
    utils::write.table(
        experiment,
        paste0(projectFolderName, "/experiment.txt"),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )
}


# Create a tamplate for the motifs.txt
getMotifsTamplate <- function(projectFolderName) {
    motifs <-
        data.frame(matrix(nrow = 0, ncol = 3))
    colnames(motifs) <- c("id", "motif", "length")
    utils::write.table(
        motifs,
        paste0(projectFolderName, "/motifs.txt"),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

}

# Create a tamplate for the traits.txt
getTraitsTamplate <- function(projectFolderName) {
    traits <-
        data.frame(matrix(nrow = 0, ncol = 1))
    colnames(traits) <- c("id")
    utils::write.table(
        traits,
        paste0(projectFolderName, "/traits.txt"),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )
}

# Create a tamplate for the miRs.txt
getMiRsTamplate <- function(projectFolderName) {
    miRs <-
        data.frame(matrix(nrow = 0, ncol = 1))
    colnames(miRs) <- c("id")
    utils::write.table(
        miRs,
        paste0(projectFolderName, "/miRs.txt"),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )
}

# Create a tamplate for the transcripts.txt
getTranscriptsTamplate <- function(projectFolderName) {
    transcripts <-
        data.frame(matrix(nrow = 0, ncol = 1))
    colnames(transcripts) <- c("id")
    utils::write.table(
        transcripts,
        paste0(projectFolderName, "/transcripts.txt"),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

}


# If the function you are looking for is not here check supportFunction.R
# Functions in supportFunction.R are used by multiple functions.
