#' @title Initialize the project folder
#'
#' @description The function initCircRNAprofiler() initializes the project
#' forlder.
#'
#' @param projectFolderName A character string specifying the name of the
#' project folder.
#'
#' @param detectionTools A character vector specifying the tools
#' used for circRNA detection. The following options are allowed: mapsplice,
#' nclscan, knife, circexplorer2, circmarker, uroborus and other. The user can
#' choose 1 or multiple tools. Subfolders named with the tools specified will
#' be generated under the project folder.
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

        # Create a template for the experiment file that the user needs to
        # fill in
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

        # Create a tamplate for the motifs file that the user needs to fill in
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

        # Create a tamplate for the traits file that the user needs to fill in
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


        # Create a tamplate for the miRs file that the user needs to fill in
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

        # Create a tamplate for the transcripts file that the user needs to
        # fill in
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
