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
#' \dontrun{
#' initCircRNAprofiler(projectFolderName = "circProject",
#'     detectionTools = "mapsplice")
#' }
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
            dir.create(file.path(projectFolderName, detectionTools[i]))
        }

        # Create a template for the experiment.txt
        .getExperimentTemplate(projectFolderName)

        # Create a template for the motifs.txt
        .getMotifsTemplate(projectFolderName)

        # Create a template for the traits.txt
        .getTraitsTemplate(projectFolderName)


        # Create a template for the miRs.txt
        .getMiRsTemplate(projectFolderName)

        # Create a template for the transcripts.txt
        .getTranscriptsTemplate(projectFolderName)
    }

# Experiment template
.getExperimentTemplate <- function(projectFolderName) {
    experiment <-
        data.frame(matrix(nrow = 0, ncol = 3))
    colnames(experiment) <-
        c("label", "fileName", "condition")
    file.create
    utils::write.table(
        experiment,
        file.path(projectFolderName, "experiment.txt"),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )
}


# Create a template for the motifs.txt
.getMotifsTemplate <- function(projectFolderName) {
    motifs <-
        data.frame(matrix(nrow = 0, ncol = 3))
    colnames(motifs) <- c("id", "motif", "length")
    utils::write.table(
        motifs,
        file.path(projectFolderName, "motifs.txt"),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

}

# Create a template for the traits.txt
.getTraitsTemplate <- function(projectFolderName) {
    traits <-
        data.frame(matrix(nrow = 0, ncol = 1))
    colnames(traits) <- c("id")
    utils::write.table(
        traits,
        file.path(projectFolderName, "traits.txt"),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )
}

# Create a template for the miRs.txt
.getMiRsTemplate <- function(projectFolderName) {
    miRs <-
        data.frame(matrix(nrow = 0, ncol = 1))
    colnames(miRs) <- c("id")
    utils::write.table(
        miRs,
        file.path(projectFolderName, "miRs.txt"),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )
}

# Create a template for the transcripts.txt
.getTranscriptsTemplate <- function(projectFolderName) {
    transcripts <-
        data.frame(matrix(nrow = 0, ncol = 1))
    colnames(transcripts) <- c("id")
    utils::write.table(
        transcripts,
        file.path(projectFolderName, "transcripts.txt"),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

}


# If the function you are looking for is not here check supportFunction.R
# Functions in supportFunction.R are used by multiple functions.
