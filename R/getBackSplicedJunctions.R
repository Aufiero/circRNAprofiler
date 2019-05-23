#' @title Import detected circRNAs
#'
#' @description The function getBackSplicedJunctions() reads the
#' circRNAs_X.txt with the detected circRNAs, adapts the content and generates
#' a unique data frame with all circRNAs identified by each circRNA detection
#' tool and the occurrences found in each sample (named as reported in the
#' column label in experiment.txt).
#'
#' @param gtf A data frame containing the annotation information. It can be
#' generated with \code{\link{formatGTF}}.
#'
#' @param pathToExperiment A string containing the path to the experiment.txt
#' file. The file experiment.txt contains the experiment design information.
#' It must have at least 3 columns with headers:
#' - label (1st column): unique names of the samples (short but informative).
#' - fileName (2nd column): name of the input files - e.g. circRNAs_X.txt, where
#' x can be can be 001, 002 etc.
#' - group (3rd column): biological conditions - e.g. A or B; healthy or
#' diseased if you have only 2 conditions.
#'
#' By default pathToExperiment i set to NULL and the file it is searched in
#' the working directory. If experiment.txt is located in a different directory
#' then the path needs to be specified.
#'
#' @return A data frame.
#'
#' @examples
#' check <- checkProjectFolder()
#'
#' if(check == 0){
#' # Create gtf object
#' gtf <- formatGTF(pathToGTF)
#'
#' # Read and adapt detected circRNAs
#' backSplicedJunctions<- getBackSplicedJunctions(gtf)}
#'
#' @seealso
#' \code{\link{backSplicedJunctions}} for a description of the data frame
#' containing back-spliced junctions coordinates.
#'
#' @import dplyr
#' @importFrom magrittr %>%
#' @importFrom utils read.table
#' @importFrom rlang .data
#' @export
getBackSplicedJunctions <-  function(gtf, pathToExperiment = NULL) {
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

        fileNames <- list.files()

        # Retrieve the code for each circRNA prediction tool
        detectionTools <- getDetectionTools()

        # Keep only the circRNA prediction tools that have been used for
        # circRNA detection
        detectionToolsUsed <- detectionTools %>%
            dplyr::filter(.data$name %in% fileNames)

        if (nrow(detectionToolsUsed) > 0) {
            # Get basic colum names
            basicColumns <- getBasicColNames()
            # Create the data frame that will be filled with the circRNA prediction
            # perfomed by the prediction tools used.
            backSplicedJunctions <-
                data.frame(matrix(nrow = 0, ncol = length(c(
                    basicColumns, "tool"
                ))))
            colnames(backSplicedJunctions) <-
                c(basicColumns, "tool")


            for (j in seq_along(detectionToolsUsed$name)) {
                # Create en empty data frame to be filled the circRNA prediction of
                # each tool at a time.
                backSplicedJunctionsTool <-
                    data.frame(matrix(
                        nrow = 0,
                        ncol = length(basicColumns)
                    ))
                colnames(backSplicedJunctionsTool) <- basicColumns

                for (i in seq_along(experiment$fileName)) {
                    # Read the file contaning the prediction one at the time
                    pathToFile <-
                        paste(detectionToolsUsed$name[j],
                            experiment$fileName[i],
                            sep = "/")

                    # A specific import function is called
                    adaptedPatientBSJunctions <-
                        switch(
                            detectionToolsUsed$name[j],
                            mapsplice = importMapSplice(pathToFile, gtf),
                            nclscan = importNCLscan(pathToFile),
                            knife = importKnife(pathToFile),
                            circexplorer2 = importCircExplorer2(pathToFile),
                            other = importOther(pathToFile),
                            circmarker = importCircMarker (pathToFile, gtf),
                            uroborus <- importUroborus(pathToFile)
                        )

                    # Check validity of adaptedPatientBSJunctions.
                    adaptedPatientBSJunctions <-
                        checkBSJsDF(adaptedPatientBSJunctions, addColNames = "coverage")

                    patientBSJunctions <- adaptedPatientBSJunctions

                    indexCoverage <-
                        which(colnames(patientBSJunctions) == "coverage")

                    colnames(patientBSJunctions)[indexCoverage] <-
                        experiment$label[i]

                    # Add a new row if a new circRNA is found and the corresponding
                    # columns with all circRNA coverages
                    backSplicedJunctionsTool <- base::merge(
                        backSplicedJunctionsTool,
                        patientBSJunctions,
                        by = basicColumns,
                        all = TRUE,
                        sort = FALSE
                    )


                }

                tool <-
                    rep(detectionToolsUsed$code[j],
                        nrow(backSplicedJunctionsTool))

                # Repalce NA values with 0 (zero)
                backSplicedJunctionsTool <-
                    backSplicedJunctionsTool %>%
                    dplyr::mutate_at(experiment$label, ~ replace(., is.na(.), 0)) %>%
                    dplyr::mutate(tool = tool)

                #bind the data frame containing circRNA prediction perfomed by the
                # circRNA detection tools
                backSplicedJunctions <-
                    dplyr::bind_rows(backSplicedJunctions,
                        backSplicedJunctionsTool)

            }

            # Round
            backSplicedJunctions[, experiment$label] <-
                round(backSplicedJunctions[, experiment$label], 0)
        } else{
            backSplicedJunctions <- data.frame()
            cat(
                "Missing folders with circRNA_X.txt files. Check working directory
                or type getDetectionTools() to see the name of the circRNA
                detection tools."
            )
        }


        } else{
            backSplicedJunctions <- data.frame()
            cat("experiment.txt not found. The analysis can not start.")
        }


    return(backSplicedJunctions)
    }



#' @title Create data frame with circRNA detection codes
#'
#' @description The function getDetectionTools() creates a data frame
#' containing the codes corresponding to each circRNA detection tool for which
#' a specific import function has been developed.
#'
#' @return A data frame
#'
#' @examples
#' getDetectionTools()
#'
#' @export
getDetectionTools <- function() {
    # Create a data frame
    detectionTools <- data.frame(matrix(nrow = 7, ncol = 2))
    colnames(detectionTools) <- c("name", "code")

    detectionTools$name[1] <- "mapsplice"
    detectionTools$code[1] <- "ms"

    detectionTools$name[2] <- "nclscan"
    detectionTools$code[2] <- "ns"

    detectionTools$name[3] <- "circexplorer2"
    detectionTools$code[3] <- "ce"

    detectionTools$name[4] <- "knife"
    detectionTools$code[4] <- "kn"

    detectionTools$name[5] <- "other"
    detectionTools$code[5] <- "ot"

    detectionTools$name[6] <- "circmarker"
    detectionTools$code[6] <- "cm"

    detectionTools$name[7] <- "uroborus"
    detectionTools$code[7] <- "ur"

    return(detectionTools)

}


#' @title Group circRNAs identified by multiple prediction tools
#'
#' @description The function mergeBSJunctions() shrinks the data frame by
#' grouping back-spliced junctions commonly identified by multiple
#' detection tools. The read counts of the samples reported in the final
#' data frame will be the ones of the tool that detected the highest total mean
#' across all samples. All the tools that detected the back-spliced junctions
#' are then listed in the column "tool" of the final data frame.
#' See \code{\link{getDetectionTools}}  for more detail about the code
#' corresponding to each circRNA detection tool.
#'
#' @param backSplicedJunctions A data frame containing back-spliced junction
#' coordinates and counts generated with \code{\link{getBackSplicedJunctions}}.
#'
#' @param gtf A data frame containing genome annotation information,
#' generated with \code{\link{formatGTF}}.
#'
#' @param pathToExperiment A string containing the path to the experiment.txt
#' file. The file experiment.txt contains the experiment design information.
#' It must have at least 3 columns with headers:
#' - label (1st column): unique names of the samples (short but informative).
#' - fileName (2nd column): name of the input files - e.g. circRNAs_X.txt, where
#' x can be can be 001, 002 etc.
#' - group (3rd column): biological conditions - e.g. A or B; healthy or
#' diseased if you have only 2 conditions.
#'
#' By default pathToExperiment i set to NULL and the file it is searched in
#' the working directory. If experiment.txt is located in a different directory
#' then the path needs to be specified.
#'
#' @param antisense A logical specifying whether to export the identified antisense
#' circRNAs in a file named antisenseCircRNAs.txt. Default value is FALSE.
#' A circRNA is defined antisense if the strand reported in the prediction
#' results is different from the strand reported in the genome annotation file.
#' The antisense circRNAs are removed from the returned data frame.
#'
#' @return A data frame.
#'
#' @examples
#' # Load detected back-soliced junctions
#' data("backSplicedJunctions")
#'
#' # Load short version of the gencode v19 annotation file
#' data("gtf")
#'
#' pathToExperiment <- system.file("extdata", "experiment.txt",
#'     package ="circRNAprofiler")
#'
#' # Merge commonly identified circRNAs
#' mergedBSJunctions <- mergeBSJunctions(backSplicedJunctions, gtf,
#'     pathToExperiment)
#'
#' @importFrom magrittr %>%
#' @importFrom utils read.table
#' @importFrom utils write.table
#' @importFrom rlang .data
#' @import dplyr
#' @export
mergeBSJunctions <-
    function(backSplicedJunctions,
        gtf,
        pathToExperiment = NULL,
        antisense = FALSE) {
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
            detectionTools <- getDetectionTools()
            # Find and merge commonly identified back-spliced junctions
            mergedBSJunctions <- backSplicedJunctions %>%
                dplyr::mutate(mean = rowMeans(.[, experiment$label])) %>%
                dplyr::group_by(.data$strand,
                    .data$chrom,
                    .data$startUpBSE,
                    .data$endDownBSE) %>%
                dplyr::arrange(desc(mean)) %>%
                dplyr::mutate(mergedTools = paste(sort(.data$tool), collapse = ",")) %>%
                dplyr::filter(row_number() == 1) %>%
                dplyr::ungroup() %>%
                dplyr::select(-c(.data$tool, .data$mean)) %>%
                dplyr::rename(tool = .data$mergedTools) %>%
                dplyr::select(
                    .data$id,
                    .data$gene,
                    .data$strand,
                    .data$chrom,
                    .data$startUpBSE,
                    .data$endDownBSE,
                    .data$tool,
                    everything()
                ) %>%
                as.data.frame()
            # dplyr::mutate(mergedTools = paste(sort(detectionTools$code[match(.data$tool, detectionTools$name)]), collapse = ","))

            # For some circRNAs the strand reported in prediction results is
            # sometimes different from the strand reported in the gtf file.
            shrinkedGTF <- gtf %>%
                dplyr::select(.data$gene_name, .data$strand) %>%
                dplyr::group_by(.data$gene_name) %>%
                dplyr::filter(row_number() == 1)

            mt <-
                match(mergedBSJunctions$gene, shrinkedGTF$gene_name)
            antisenseCircRNAs <-
                dplyr::bind_cols(mergedBSJunctions, shrinkedGTF[mt, ]) %>%
                dplyr::filter(.data$strand != .data$strand1) %>%
                dplyr::select(-c(.data$gene_name, .data$strand1))

            if(antisense){
                utils::write.table(
                    antisenseCircRNAs,
                    "antisenseCircRNAs.txt",
                    quote = FALSE,
                    row.names = FALSE,
                    col.names = TRUE,
                    sep = "\t"
                )
            }


            # Remove from the dataframe the antisense circRNAs
            mergedBSJunctionsClenead <- mergedBSJunctions %>%
                dplyr::filter(!(mergedBSJunctions$id %in% antisenseCircRNAs$id))

        } else{
            mergedBSJunctionsClenead <- backSplicedJunctions
            cat("experiment.txt not found, data frame can not be merged")
        }

        return(mergedBSJunctionsClenead)
    }
