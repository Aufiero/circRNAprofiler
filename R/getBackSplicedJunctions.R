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
#' \describe{
#' \item{label:}{(1st column) - unique names of the samples (short but informative).}
#' \item{fileName:}{(2nd column) - name of the input files - e.g. circRNAs_X.txt, where
#' x can be can be 001, 002 etc.}
#' \item{group:}{ (3rd column) - biological conditions - e.g. A or B; healthy or diseased
#' if you have only 2 conditions.}
#' }
#'
#' By default pathToExperiment is set to NULL and the file it is searched in
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
    # Read experiment.txt
    experiment <- .readExperiment(pathToExperiment)
    if (nrow(experiment)) {
        fileNames <- list.files()
        # Retrieve the code for each circRNA prediction tool
        detectionTools <- getDetectionTools()
        # Keep tools that have been used for circRNA detection
        detectionToolsUsed <- detectionTools %>%
            dplyr::filter(.data$name %in% fileNames)
        if (nrow(detectionToolsUsed) > 0) {
            # Create backSplicedJunctions data frame
            backSplicedJunctions <-
                .createBackSplicedJunctionsDF(addColNames = "tool")
            for (j in seq_along(detectionToolsUsed$name)) {
                # data frame to store circRNA prediction
                backSplicedJunctionsTool <- .createBackSplicedJunctionsDF()

                for (i in seq_along(experiment$fileName)) {
                    # Read the file contaning the prediction one at the time
                    pathToFile <- file.path(detectionToolsUsed$name[j],
                        experiment$fileName[i])

                    nameTool <- detectionToolsUsed$name[j]
                    # A specific import function is called
                    adaptedPatientBSJunctions <-
                        .getAdaptedPatientBSJunctions(nameTool, pathToFile, gtf)

                    # Check validity of adaptedPatientBSJunctions.
                    adaptedPatientBSJunctions <-
                        .checkBSJsDF(adaptedPatientBSJunctions, addColNames = "coverage")

                    patientBSJunctions <- adaptedPatientBSJunctions
                    indexCoverage <- which(colnames(patientBSJunctions) == "coverage")
                    colnames(patientBSJunctions)[indexCoverage] <- experiment$label[i]

                    # Merge circRNAs
                   basicColumns <- .getBasicColNames()
                    backSplicedJunctionsTool <- base::merge(
                        backSplicedJunctionsTool,
                        patientBSJunctions,
                        by = basicColumns,
                        all = TRUE,
                        sort = FALSE
                    )
                }
                tool <- rep(detectionToolsUsed$code[j],
                        nrow(backSplicedJunctionsTool))
                # Repalce NA values with 0 (zero)
                backSplicedJunctionsTool <- backSplicedJunctionsTool %>%
                    dplyr::mutate_at(experiment$label, ~ replace(., is.na(.), 0)) %>%
                    dplyr::mutate(tool = tool)

                # rbind the data frame containing circRNA prediction
                backSplicedJunctions <- dplyr::bind_rows(backSplicedJunctions,
                        backSplicedJunctionsTool)
            }
            # Round
            backSplicedJunctions[, experiment$label] <-
                round(backSplicedJunctions[, experiment$label], 0)
        } else{
            backSplicedJunctions <- data.frame()
            cat( "Missing folders with circRNA_X.txt files. Check working directory
                or type getDetectionTools() to see the name of the circRNA
                detection tools."
            )
        }
    } else{
        backSplicedJunctions <- data.frame()
        cat("experiment.txt not found or empty. The analysis can not start.")
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
#' @param exportAntisense A logical specifying whether to export the identified
#' antisense circRNAs in a file named antisenseCircRNAs.txt. Default value is
#' FALSE. A circRNA is defined antisense if the strand reported in the prediction
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
        exportAntisense = FALSE) {
        # Read experiment.txt
        experiment <- .readExperiment(pathToExperiment)
        if (nrow(experiment) > 0) {
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

            # Identified antisense circRNAs
            antisenseCircRNAs <-
                .getAntisenseCircRNAs(mergedBSJunctions, gtf, exportAntisense)
            # Remove from the dataframe the antisense circRNAs
            mergedBSJunctionsClenead <- mergedBSJunctions %>%
                dplyr::filter(!(mergedBSJunctions$id %in% antisenseCircRNAs$id))

        } else{
            mergedBSJunctionsClenead <- backSplicedJunctions
            cat("experiment.txt not found or empty, data frame can not be merged")
        }

        return(mergedBSJunctionsClenead)
    }


# Create backSplicedJunctions data frame
.createBackSplicedJunctionsDF <- function(addColNames = NULL) {
    # Get basic colum names
    basicColumns <- .getBasicColNames()
    # Create the data frame that will be filled with the circRNA prediction
    # perfomed by the prediction tools used.
    backSplicedJunctions <-
        data.frame(matrix(nrow = 0, ncol = length(c(
            basicColumns, addColNames
        ))))
    colnames(backSplicedJunctions) <-
        c(basicColumns, addColNames)

    return(backSplicedJunctions)
}





# Get the adaptedPatientBSJunctions data frame with circRNA predictions
.getAdaptedPatientBSJunctions <-
    function(nameTool, pathToFile, gtf) {
        # A specific import function is called
        adaptedPatientBSJunctions <-
            switch(
                nameTool,
                mapsplice = importMapSplice(pathToFile),
                nclscan = importNCLscan(pathToFile),
                knife = importKnife(pathToFile),
                circexplorer2 = importCircExplorer2(pathToFile),
                other = importOther(pathToFile),
                circmarker = importCircMarker(pathToFile, gtf),
                uroborus <- importUroborus(pathToFile)
            )

        return(adaptedPatientBSJunctions)
    }

# For some circRNAs the strand reported in prediction results is
# sometimes different from the strand reported in the gtf file.
# With this function we identified the antisense circRNAs
.getAntisenseCircRNAs <-
    function(mergedBSJunctions, gtf, exportAntisense = FALSE) {
        shrinkedGTF <- gtf %>%
            dplyr::select(.data$gene_name, .data$strand) %>%
            dplyr::group_by(.data$gene_name) %>%
            dplyr::filter(row_number() == 1)

        mt <-
            match(mergedBSJunctions$gene, shrinkedGTF$gene_name)
        antisenseCircRNAs <-
            dplyr::bind_cols(mergedBSJunctions, shrinkedGTF[mt,]) %>%
            dplyr::filter(.data$strand != .data$strand1) %>%
            dplyr::select(-c(.data$gene_name, .data$strand1))

        if (exportAntisense) {
            utils::write.table(
                antisenseCircRNAs,
                "antisenseCircRNAs.txt",
                quote = FALSE,
                row.names = FALSE,
                col.names = TRUE,
                sep = "\t"
            )
        }

        return(antisenseCircRNAs)
    }


# If the function you are looking for is not here check supportFunction.R
# Functions in supportFunction.R are used by multiple functions.
