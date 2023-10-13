#' @title Annotate circRNA features.
#'
#' @description The function annotateBSJs() annotates the circRNA structure
#' and the introns flanking the corresponding back-spliced junctions.
#' The genomic features are extracted from the user provided gene annotation.
#'
#' @param backSplicedJunctions A data frame containing the back-spliced junction
#' coordinates (e.g. detected or randomly selected).
#' For detected back-spliced junctions see \code{\link{getBackSplicedJunctions}}
#' and \code{\link{mergeBSJunctions}} (to group circRNAs detected by multiple
#' detection tools). For randomly selected back-spliced junctions see
#' \code{\link{getRandomBSJunctions}}.
#'
#' @param gtf A data frame containing genome annotation information. It can be
#' generated with \code{\link{formatGTF}}.
#'
#' @param isRandom A logical indicating whether the back-spliced junctions have
#' been randomly generated with \code{\link{getRandomBSJunctions}}.
#' Deafult value is FALSE.
#'
#' @param pathToTranscripts A string containing the path to the transcripts.txt
#' file. The file transcripts.txt contains the transcript ids of the
#' circRNA host gene to analyze. It must have one column with header id.
#' By default pathToTranscripts is set to NULL and the file it is searched in
#' the working directory. If transcripts.txt is located in a different directory
#' then the path needs to be specified. If this file is empty or absent the
#' longest transcript of the circRNA host gene containing the back-spliced
#' junctions are considered in the annotation analysis.
#'
#' @return A data frame.
#'
#' @examples
#' # Load a data frame containing detected back-spliced junctions
#' data("mergedBSJunctions")
#'
#' # Load short version of the gencode v19 annotation file
#' data("gtf")
#'
#' # Annotate the first back-spliced junctions
#' annotatedBSJs <- annotateBSJs(mergedBSJunctions[1, ], gtf)
#'
#' @importFrom magrittr %>%
#' @importFrom stringr str_detect
#' @importFrom rlang .data
#' @importFrom utils read.table
#' @import dplyr
#' @export
annotateBSJs <- function(backSplicedJunctions,
    gtf,
    isRandom = FALSE,
    pathToTranscripts = NULL) {
    # Read trancripts.txt
    transcriptsFromFile <- .readTranscripts(pathToTranscripts)
    if (nrow(transcriptsFromFile) == 0 & !isRandom) {
        cat("transcripts.txt is empty or absent. The longest
                transcripts for all circRNAs will be analyzed") }
    # Check the validity and the content
    backSplicedJunctions <- .checkBSJsDF(backSplicedJunctions)
    # Create an empty data frame to store the extracted information
    annotatedBSJs <- .createAnnotatedBSJsDF(backSplicedJunctions)

    # Since the coordinates of the detected back-spliced junctions might
    # not exactly correspond to annotated exonic coordinates, a gap of 10
    # nucleotides is allowed.
    for (i in seq_along(backSplicedJunctions$id)) {
        exJ1 <- substr(backSplicedJunctions[i, "startUpBSE"],
            1,
            nchar(backSplicedJunctions[i, "startUpBSE"]) - 1)
        exJ2 <- substr(backSplicedJunctions[i, "endDownBSE"],
            1,
            nchar(backSplicedJunctions[i, "endDownBSE"]) - 1)
        gene <- backSplicedJunctions$gene[i]

        annotatedBSJs$id[i] <- backSplicedJunctions$id[i]
        annotatedBSJs$gene[i] <- backSplicedJunctions$gene[i]
        annotatedBSJs$strand[i] <- backSplicedJunctions$strand[i]
        annotatedBSJs$chrom[i] <- backSplicedJunctions$chrom[i]

        # Retrieves all transcripts whose exon coordinates overlap that of
        # the given back-spliced junction coordinates (exJ1, exJ2)
        allTranscripts <-  .getAllTranscripts(gtf, gene, exJ1, exJ2)
        if (!is.na(allTranscripts[1])) {
            annotatedBSJs$allTranscripts[i] <- paste(allTranscripts, collapse = ",")
            # Get transcript to analyze
            transcriptToAnalyze <- .getTranscriptToAnalyze(transcriptsFromFile, allTranscripts, gtf)
            annotatedBSJs$transcript[i] <- transcriptToAnalyze$transcript_id[1]
            annotatedBSJs$totExons[i] <- nrow(transcriptToAnalyze)

            # get BSEs from transcriptToAnalyze
            bsExons <- .getBSEsFromTranscript(transcriptToAnalyze, exJ1, exJ2)
            annotatedBSJs$exNumUpBSE[i] <- bsExons$exon_number[1]
            annotatedBSJs$exNumDownBSE[i] <- bsExons$exon_number[2]
            annotatedBSJs$numOfExons[i] <- length(c(bsExons$exon_number[1]:bsExons$exon_number[2]))
            if (bsExons$strand[1] == "-") {
                annotatedBSJs$startUpBSE[i] <- bsExons$end[1]
                annotatedBSJs$endUpBSE[i] <- bsExons$start[1]
                annotatedBSJs$startDownBSE[i] <- bsExons$end[2]
                annotatedBSJs$endDownBSE[i] <- bsExons$start[2]
                # For negative strand
            } else if (bsExons$strand[1] == "+") {
                annotatedBSJs$startUpBSE[i] <- bsExons$start[1]
                annotatedBSJs$endUpBSE[i] <- bsExons$end[1]
                annotatedBSJs$startDownBSE[i] <- bsExons$start[2]
                annotatedBSJs$endDownBSE[i] <- bsExons$end[2]
            }
            # Retrive flanking intron coordinates
            intronCoords <- .getIntronCoords(transcriptToAnalyze, bsExons)
            annotatedBSJs$startUpIntron[i] <- intronCoords$startUpIntron[1]
            annotatedBSJs$endUpIntron[i] <- intronCoords$endUpIntron[1]
            annotatedBSJs$startDownIntron[i] <- intronCoords$startDownIntron[1]
            annotatedBSJs$endDownIntron[i] <- intronCoords$endDownIntron[1]

            # Calculate circRNA length (exon only)
            lenCircRNA <- .getLengthCirc(transcriptToAnalyze, bsExons)
            annotatedBSJs$lenCircRNA[i] <- lenCircRNA[[1]]
        } else if (is.na(allTranscripts[1])) {
            annotatedBSJs$startUpBSE[i] <-  backSplicedJunctions$startUpBSE[i]
            annotatedBSJs$endDownBSE[i] <- backSplicedJunctions$endDownBSE[i]
        }
    }
    # Retrieve the length (bp) of the bse and flanking introns
    lenBSEfi <- .getLengthBSEfi(annotatedBSJs)
    annotatedBSJs$lenUpIntron <- lenBSEfi$lenUpIntron
    annotatedBSJs$lenUpBSE <- lenBSEfi$lenUpBSE
    annotatedBSJs$lenDownBSE <- lenBSEfi$lenDownBSE
    annotatedBSJs$lenDownIntron <- lenBSEfi$lenDownIntron
    annotatedBSJs$meanLengthBSEs <- lenBSEfi$meanLengthBSEs
    annotatedBSJs$meanLengthIntrons <- lenBSEfi$meanLengthIntrons
    return(annotatedBSJs)
}


# The function getAnnotatedBSJsColNames() returns the column
.getAnnotatedBSJsColNames <- function() {
    basicColumns <- .getBasicColNames()

    annotatedBSJsColNames <-
        c(
            basicColumns[1],
            basicColumns[2],
            "allTranscripts",
            "transcript",
            "totExons",
            basicColumns[3],
            basicColumns[4],
            c("startUpIntron",
                "endUpIntron"),
            basicColumns[5],
            c("endUpBSE",
                "startDownBSE"),
            basicColumns[6],
            c(
                "startDownIntron",
                "endDownIntron",
                "exNumUpBSE",
                "exNumDownBSE",
                "numOfExons",
                "lenUpIntron",
                "lenUpBSE",
                "lenDownBSE",
                "lenDownIntron",
                "meanLengthBSEs",
                "meanLengthIntrons",
                "lenCircRNA"
            )
        )

    return(annotatedBSJsColNames)

}


# The function getAllTranscripts() retrieves all transcripts whose exon
# coordinates overlap that of the given back-spliced junctions coordinates
# (exJ1 and exJ2). If the trancripts are not found it
# might be that the annotation file is different from the one used
# for mapping or it can be due to the circRNA prediction tool
# analysis. In both cases since the genomic features can
# not be extrated an NA value is reported.
.getAllTranscripts <- function(gtf, gene, exJ1, exJ2) {
    # Find the isofom containing the back-spliced exon junction coordinates
    #  given in input
    start <- gtf %>%
        dplyr::filter(gene_name == gene,
            type == "exon") %>%
        dplyr::group_by(transcript_id) %>%
        dplyr::filter(
            stringr::str_detect(start, exJ1) |
                stringr::str_detect(start, exJ2)
        )

    # Find the isofom containing the back-spliced exon junction coordinates
    # given in input
    end <- gtf %>%
        dplyr::filter(gene_name == gene,
            type == "exon") %>%
        dplyr::group_by(transcript_id) %>%
        dplyr::filter(stringr::str_detect(end, exJ1) |
                stringr::str_detect(end, exJ2))

    # Check whether the end and the start coordinates are present and
    # if so take the the transcripts and sort by lenght. The first transcript
    # will be the longest one
    inter <- intersect(start$transcript_id, end$transcript_id)
    if (length(inter) != 0) {
        allTranscripts <- gtf %>%
            dplyr::filter(transcript_id %in% inter) %>%
            dplyr::group_by(transcript_id) %>%
            dplyr::summarise(len = sum(width)) %>%
            dplyr::arrange(desc(len)) %>%
            .$transcript_id

    } else{
        # If one of the 2 back-spliced junction is not among the exons annotated
        # in the gtf file given in input, it can be that a different gtf file
        # has been used during the mapping procedure to idetified the circRNAs.
        allTranscripts <-  NA_character_

    }

    # Return a data frame with the genomic coordinates of all the exons of
    # the longest isoform containing the given back-spliced exon junction
    # coordinates (exJ1 and exJ2)
    return(allTranscripts)
}


# The function getFlankIntrons() retrieves the genomic coordinates
# of the introns flanking the back-spliced exons given in input.
# This function is called when the back-spliced exons are not the first and
# last exon of the transcript.
.getFlankIntrons <- function(transcriptToAnalyze, bsExons) {
    # Create an empty data frame
    intronCoords <- data.frame(matrix(
        nrow = 1,
        ncol = length(.getAnnotatedBSJsColNames())
    ))
    colnames(intronCoords) <- .getAnnotatedBSJsColNames()

    # Retrieve the upstream and the downstream exons (adiacent exons) flanking
    # the back-spliced given in input.
    adiacentExons <- transcriptToAnalyze %>%
        dplyr::filter(exon_number %in% c(bsExons$exon_number[1] - 1,
            bsExons$exon_number[2] + 1))

    # Retrive the intron coordinates flanking the back-spliced exons given in
    # input
    # for positive strand
    if (bsExons$strand[1] == "+") {
        intronCoords$startUpIntron <- adiacentExons$end[1] + 1
        intronCoords$endUpIntron <- bsExons$start[1] - 1
        intronCoords$startDownIntron <- bsExons$end[2] + 1
        intronCoords$endDownIntron <- adiacentExons$start[2] - 1

        # for negative strand
    } else if (bsExons$strand[1] == "-") {
        intronCoords$startUpIntron <- adiacentExons$start[1] - 1
        intronCoords$endUpIntron <- bsExons$end[1] + 1
        intronCoords$startDownIntron <- bsExons$start[2] - 1
        intronCoords$endDownIntron <- adiacentExons$end[2] + 1
    }

    return(intronCoords)

}


# The function getFlankIntronFirst() retrieves the genomic
# coordinates of the introns flanking the back-spliced exons given in input.
# This function is called when one of the two back-spliced exons is the
# first of the transcript.
.getFlankIntronFirst <- function(transcriptToAnalyze, bsExons) {
    # Create an empty data frame with 1 row
    intronCoords <- data.frame(matrix(
        nrow = 1,
        ncol = length(.getAnnotatedBSJsColNames())
    ))
    colnames(intronCoords) <- .getAnnotatedBSJsColNames()

    # Retrieve the upstream and the downstream exons (adiacent exons) flanking
    # the back-spliced exons given in input.
    adiacentExons <- transcriptToAnalyze %>%
        dplyr::filter(exon_number ==
                bsExons$exon_number[2] + 1)

    # Retrive the coordinates of the downstream intron.
    # Since one of the exon is the first exon of the transcritpt the upstream
    # introns is not available "NA"

    # for positive strand
    if (bsExons$strand[1] == "+") {
        intronCoords$startDownIntron <- bsExons$end[2] + 1
        intronCoords$endDownIntron <- adiacentExons$start[1] - 1


        # for negative strand
    } else if (bsExons$strand[1] == "-") {
        intronCoords$startDownIntron <- bsExons$start[2] - 1
        intronCoords$endDownIntron <- adiacentExons$end[1] + 1
    }

    return(intronCoords)

}


# The function getFlankIntronLast() retrieves the coordinates of
# the introns flanking the back-spliced exons given in input.
# This function is called when one of the two back-spliced exons is the last
# of the transcript.
.getFlankIntronLast <- function(transcriptToAnalyze, bsExons) {
    # Create an empty data frame with 1 row
    intronCoords <- data.frame(matrix(
        nrow = 1,
        ncol = length(.getAnnotatedBSJsColNames())
    ))
    colnames(intronCoords) <- .getAnnotatedBSJsColNames()

    # Retrieve the upstream and the downstream exons (adiacent exons) flanking
    # the back-spliced exons given in input.
    adiacentExons <- transcriptToAnalyze %>%
        dplyr::filter(exon_number ==
                bsExons$exon_number[1] - 1)

    # Retrive the coordinates of the upstream intron
    # Since one the exon is the last exon of the transcritpt the downstream
    # intron is not available "NA"

    # For positive strand
    if (bsExons$strand[1] == "+") {
        intronCoords$startUpIntron <- adiacentExons$end[1] + 1
        intronCoords$endUpIntron <- bsExons$start[1] - 1

        # For negative strand
    } else if (bsExons$strand[1] == "-") {
        intronCoords$startUpIntron <- adiacentExons$start[1] - 1
        intronCoords$endUpIntron <- bsExons$end[1] + 1
    }
    return(intronCoords)

}


# Get coordinates of flanking introns
.getIntronCoords <- function(transcriptToAnalyze, bsExons) {

    # Create an empty data frame with 1 row
    intronCoords <- data.frame(matrix(
        nrow = 1,
        ncol = length(.getAnnotatedBSJsColNames())
    ))
    colnames(intronCoords) <- .getAnnotatedBSJsColNames()

    # Retrieve the coordinates of the flanking introns
    if (bsExons$exon_number[1] != 1 &
            bsExons$exon_number[2] != nrow(transcriptToAnalyze)) {
        # If the BSEs are not the first and last of the transcript
        intronCoords <-
            .getFlankIntrons(transcriptToAnalyze, bsExons)
        # If one the BSE is the first exon of the transcript
    } else if (bsExons$exon_number[1] == 1 &
            bsExons$exon_number[2] != nrow(transcriptToAnalyze)) {
        intronCoords <- .getFlankIntronFirst(transcriptToAnalyze, bsExons)
        # If one the BSE is the last exon of the transcript
    } else if (bsExons$exon_number[2] == nrow(transcriptToAnalyze) &
            bsExons$exon_number[1] != 1) {
        intronCoords <- .getFlankIntronLast(transcriptToAnalyze, bsExons)
    }

    return(intronCoords)
}


# The function getLengthBSEfi() calculates the length (nt) of the
# back-spliced exons and the corresponding flanking introns.
.getLengthBSEfi <- function(annotatedBSJs) {
    colNames <- c(
        "id",
        "lenUpIntron",
        "lenUpBSE",
        "lenDownBSE",
        "lenDownIntron",
        "meanLengthBSEs",
        "meanLengthIntrons"
    )

    # Create an empty dataframe
    lenBSEfi <-
        data.frame(matrix(
            nrow = nrow(annotatedBSJs),
            ncol = length(colNames)
        ))
    colnames(lenBSEfi) <- colNames


    # Calculate the back-spliced exons and introns length and select the needed
    # columns
    lenBSEfi <- annotatedBSJs %>%
        dplyr::mutate(
            lenUpIntron = abs(endUpIntron - startUpIntron),
            lenUpBSE = abs(endUpBSE - startUpBSE),
            lenDownBSE = abs(endDownBSE - startDownBSE),
            lenDownIntron = abs(endDownIntron - startDownIntron),
            meanLengthBSEs = base::rowMeans(base::cbind(lenUpBSE, lenDownBSE), na.rm =
                    TRUE),
            meanLengthIntrons = base::rowMeans(
                base::cbind(lenUpIntron, lenDownIntron),
                na.rm =
                    TRUE
            )
        ) %>%
        dplyr::select(
            id,
            lenUpIntron,
            lenUpBSE,
            lenDownBSE,
            lenDownIntron,
            meanLengthBSEs,
            meanLengthIntrons
        )
    # Return a data frame containing the length (bp) of the back-spliced exons
    # and the corresponding falnking introns.

    # Repalce NaN values with NA
    lenBSEfi[is.na(lenBSEfi)] <- NA_character_
    return(lenBSEfi)

}



# Create annotatedBSJs data frame
.createAnnotatedBSJsDF <- function(backSplicedJunctions) {
    annotatedBSJs <-
        data.frame(matrix(
            nrow = nrow(backSplicedJunctions),
            ncol = length(.getAnnotatedBSJsColNames())
        ))
    colnames(annotatedBSJs) <- .getAnnotatedBSJsColNames()

    return(annotatedBSJs)
}


# Get transxript to analyze
.getTranscriptToAnalyze <-
    function(transcriptsFromFile, allTranscripts, gtf) {
        if (nrow(transcriptsFromFile) > 0 &
                any(transcriptsFromFile$id %in% allTranscripts)) {
            # Analyze the transcript specified by the user in
            # transcripts.txt
            index <-
                which(allTranscripts %in%  transcriptsFromFile$id)
            # In case multiple transcripts are given in input for
            # the same circRNA, the first is taken. Only one isoform
            # at the time can be analyzed.
            # For this reason index[1]
            transcript <-  allTranscripts[index[1]]
        } else{
            # The first transcript in allTranscript is the longest
            transcript <- allTranscripts[1]
        }

        # The isoform that is analyzed can be the longest isoform
        # containing the BSJs or the isoform specified by the user
        # in transcripts.txt
        # TODO consider the exons to keep if this is known
        transcriptToAnalyze <- gtf  %>%
            dplyr::filter(transcript_id ==
                    transcript) %>%
            dplyr::arrange(exon_number)

        return(transcriptToAnalyze)
    }


# Get back-Splice exons from transcriptToAnalyze
.getBSEsFromTranscript <- function(transcriptToAnalyze, exJ1, exJ2) {
    # If the isoform is present then it retrieves the missing
    # coordinates of the back-spliced exons
    bsExons <- transcriptToAnalyze %>%
        dplyr::filter(
            stringr::str_detect(start, exJ1) |
                stringr::str_detect(end, exJ1) |
                stringr::str_detect(start, exJ2) |
                stringr::str_detect(end, exJ2)
        ) %>%
        dplyr::arrange(exon_number)

    # If only one exon is present then the row is duplicated.
    # This step is only made to simplify the code below without
    # introducing an additional if statment when a single
    # back-spliced exon is present.
    if (nrow(bsExons) == 1) {
        bsExons[2,] <- bsExons[1,]
    }

    return(bsExons)
}

# Get length circRNA
.getLengthCirc <- function(transcriptToAnalyze, bsExons) {
    # Calculate circRNA length considering only the lenght of the
    # exons in between the back-spliced junctions
    lenCircRNA <-
        transcriptToAnalyze %>%
        dplyr::filter(exon_number %in%
                c(bsExons$exon_number[1]:bsExons$exon_number[2])) %>%
        dplyr::summarise(lenCircRNA = sum(width))
    return(lenCircRNA)
}


# If the function you are looking for is not here check supportFunction.R
# Functions in supportFunction.R are used by multiple functions.
