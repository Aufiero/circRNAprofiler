#' @title Plot length introns flanking back-spliced junctions
#'
#' @description The function plotLenIntrons() generates vertical boxplots for
#' comparison of length of introns flanking the back-spliced junctions
#' (e.g. detected Vs randomly selected).
#'
#' @param annotatedFBSJs A data frame with the annotated back-spliced junctions
#' (e.g. detected). It can be generated with \code{\link{annotateBSJs}}.
#' These act as foreground back-spliced junctions.
#'
#' @param annotatedBBSJs A data frame with the annotated back-spliced junctions
#' (e.g. randomly selected). It can generated with \code{\link{annotateBSJs}}.
#' These act as background back-spliced junctions.
#'
#' @param df1Name A string specifying the name of the first data frame. This
#' will be displayed in the legend of the plot. Deafult value is "foreground".
#'
#' @param df2Name A string specifying the name of the first data frame. This
#' will be displayed in the legend of the plot. Deafult value is "background".
#'
#' @param title A character string specifying the title of the plot.
#' 
#' @param setyLim A logical specifying whether to set y scale limits.
#' If TRUE the value in ylim will be used. Deafult value is FALSE.
#'
#' @param ylim An integer specifying the lower and upper y axis limits
#' Deafult values are c(0, 8).
#'
#' @return A ggplot object.
#'
#' @examples
#' # Load data frame containing detected back-spliced junctions
#' data("mergedBSJunctions")
#'
#' # Load short version of the gencode v19 annotation file
#' data("gtf")
#'
#' # Annotate the first 10 back-spliced junctions
#' annotatedFBSJs <- annotateBSJs(mergedBSJunctions[1:10, ], gtf)
#'
#' # Get random back-spliced junctions
#' randomBSJunctions <- getRandomBSJunctions( gtf, n = 10, f = 10)
#'
#' # Annotate random back-spliced junctions
#' annotatedBBSJs <- annotateBSJs(randomBSJunctions, gtf, isRandom = TRUE)
#'
#' # Plot
#' p <- plotLenIntrons(
#'     annotatedFBSJs,
#'     annotatedBBSJs,
#'     df1Name = "foreground",
#'     df2Name = "background",
#'     title = "")
#' p
#'
#' @import ggplot2
#' @import dplyr
#' @import magrittr
#' @import reshape2
#' @importFrom rlang .data
#' @export
plotLenIntrons <-
    function(annotatedFBSJs,
        annotatedBBSJs,
        df1Name = "foreground",
        df2Name = "background",
        title = "",
        setyLim = FALSE,
        ylim = c(0, 8)) {
        # Reshape the data frame
        reshForeground <- annotatedFBSJs %>%
            dplyr::mutate(circRNA = rep("foreground", nrow(annotatedFBSJs))) %>%
            dplyr::select(
                .data$id,
                .data$gene,
                .data$circRNA,
                .data$lenUpIntron,
                .data$lenDownIntron,
                .data$meanLengthIntrons
            ) %>%
            reshape2::melt(
                id.vars = c("id", "gene", "circRNA"),
                variable.name = "feature",
                value.name = "length"
            )
        # Reshape the data frame
        reshBackground <- annotatedBBSJs %>%
            dplyr::mutate(circRNA = rep("background", nrow(annotatedBBSJs))) %>%
            dplyr::select(
                .data$id,
                .data$gene,
                .data$circRNA,
                .data$lenUpIntron,
                .data$lenDownIntron,
                .data$meanLengthIntrons
            ) %>%
            reshape2::melt(
                id.vars = c("id", "gene", "circRNA"),
                variable.name = "feature",
                value.name = "length"
            )
        combinedFB <- rbind(reshForeground, reshBackground)
        combinedFB$length <- as.numeric(combinedFB$length)
        # Plot
        if(setyLim){
            ymin <- ylim[1]
            ymax <- ylim[2]
        }else{
            ymin<- 0
            ymax<- base::max(log10(combinedFB$length), na.rm = TRUE)
        }
        
        p <-
            ggplot(combinedFB, aes(x = .data$feature, y = log10(.data$length)))+
            geom_boxplot(aes(fill = .data$circRNA), na.rm = TRUE) +
            ylim(ymin,ymax )+
            labs(title = title, x = "", y = "Log10 length (nt)") +
            theme_classic() +
            theme(plot.title = element_text(hjust = 0.5),
                axis.text.x = element_text(hjust = 0.5)) +
            scale_fill_grey(name = "circRNA", labels = c(df2Name, df1Name))
        return(p)
    }


#' @title Plot length back-spliced exons
#'
#' @description The function plotLenBSEs() generates vertical boxplots for
#' comparison of length of back-spliced exons (e.g. detected Vs randomly
#' selected).
#'
#' @param annotatedFBSJs A data frame with the annotated back-spliced junctions
#' (e.g. detected). It can be generated with \code{\link{annotateBSJs}}.
#' These act as foreground back-spliced junctions.
#'
#' @param annotatedBBSJs A data frame with the annotated back-spliced junctions
#' (e.g. randomly selected). It can generated with \code{\link{annotateBSJs}}.
#' These act as background back-spliced junctions.
#'
#' @param df1Name A string specifying the name of the first data frame. This
#' will be displayed in the legend of the plot.
#'
#' @param df2Name A string specifying the name of the first data frame. This
#' will be displayed in the legend of the plot.
#'
#' @param title A character string specifying the title of the plot
#' 
#' @param setyLim A logical specifying whether to set y scale limits.
#' If TRUE the value in ylim will be used. Deafult value is FALSE.
#'
#' @param ylim An integer specifying the lower and upper y axis limits
#' Deafult values are c(0, 8).
#'
#' @return A ggplot object.
#'
#' @examples
#' # Load data frame containing detected back-spliced junctions
#' data("mergedBSJunctions")
#'
#' # Load short version of the gencode v19 annotation file
#' data("gtf")
#'
#' # Annotate the first 10 back-spliced junctions
#' annotatedFBSJs <- annotateBSJs(mergedBSJunctions[1:10, ], gtf)
#'
#' # Get random back-spliced junctions
#' randomBSJunctions <- getRandomBSJunctions(n = 10, f = 10, gtf)
#'
#' # Annotate random back-spliced junctions
#' annotatedBBSJs <- annotateBSJs(randomBSJunctions, gtf, isRandom = TRUE)
#'
#' # Plot
#' p <- plotLenBSEs(
#'     annotatedFBSJs,
#'     annotatedBBSJs,
#'     df1Name = "foreground",
#'     df2Name = "background",
#'     title = "")
#' p
#'
#'
#' @import ggplot2
#' @import dplyr
#' @import magrittr
#' @import reshape2
#' @importFrom rlang .data
#' @export
plotLenBSEs <-
    function(annotatedFBSJs,
        annotatedBBSJs,
        df1Name = "foreground",
        df2Name = "background",
        title = "",
        setyLim = FALSE,
        ylim = c(0, 8)) {
        # Reshape the data frame
        reshForeground <- annotatedFBSJs %>%
            dplyr::mutate(circRNA = rep("foreground", nrow(annotatedFBSJs))) %>%
            dplyr::select(
                .data$id,
                .data$gene,
                .data$circRNA,
                .data$lenUpBSE,
                .data$lenDownBSE,
                .data$meanLengthBSEs
            ) %>%
            reshape2::melt(
                id.vars = c("id", "gene", "circRNA"),
                variable.name = "feature",
                value.name = "length"
            )
        # Reshape the data frame
        reshBackground <- annotatedBBSJs %>%
            dplyr::mutate(circRNA = rep("background", nrow(annotatedBBSJs))) %>%
            dplyr::select(
                .data$id,
                .data$gene,
                .data$circRNA,
                .data$lenUpBSE,
                .data$lenDownBSE,
                .data$meanLengthBSEs
            ) %>%
            reshape2::melt(
                id.vars = c("id", "gene", "circRNA"),
                variable.name = "feature",
                value.name = "length"
            )
        combinedFB <- rbind(reshForeground, reshBackground)
        combinedFB$length <- as.numeric(combinedFB$length)
        # Plot
        if(setyLim){
            ymin <- ylim[1]
            ymax <- ylim[2]
        }else{
            ymin<- 0
            ymax<- base::max(log10(combinedFB$length), na.rm = TRUE)
        }
        p <-
            ggplot(combinedFB, aes(x = .data$feature, y = log10(.data$length))) +
            geom_boxplot(aes(fill = .data$circRNA), na.rm = TRUE) +
            ylim(ymin,ymax)+
            labs(title = title, x = "", y = "Log10 length (nt)") +
            theme_classic() +
            theme(plot.title = element_text(hjust = 0.5),
                axis.text.x = element_text(hjust = 0.5)) +
            scale_fill_grey(name = "circRNA", labels = c(df2Name, df1Name))
        return(p)
    }

#' @title Plot circRNA host genes
#'
#' @description The function plotHostGenes() generates a bar chart showing the
#' no. of circRNAs produced from each the circRNA host gene.
#'
#' @param annotatedBSJs A data frame with the annotated back-spliced junctions.
#' This data frame can be generated with \code{\link{annotateBSJs}}.
#'
#' @param title A character string specifying the title of the plot.
#'
#' @return A ggplot object.
#'
#' @examples
#' # Load data frame containing detected back-spliced junctions
#' data("mergedBSJunctions")
#'
#' # Load short version of the gencode v19 annotation file
#' data("gtf")
#'
#' # Annotate the first 10 back-spliced junctions
#' annotatedBSJs <- annotateBSJs(mergedBSJunctions[1:10, ], gtf)
#'
#' # Plot
#' p <- plotHostGenes(annotatedBSJs, title = "")
#' p
#'
#' @import ggplot2
#' @import dplyr
#' @import magrittr
#' @importFrom rlang .data
#' @export
plotHostGenes <-
    function(annotatedBSJs, title = "") {
        # Plot
        p <-  annotatedBSJs %>%
            dplyr::group_by(.data$gene) %>%
            dplyr::summarise(n1 = n()) %>%
            dplyr::group_by(.data$n1) %>%
            dplyr::summarise(n2 = n()) %>%
            ggplot(aes(x = factor(.data$n1), y = factor(.data$n2))) +
            geom_bar(position = "dodge", stat = "identity") +
            labs(title = title, x = "No. of circRNAs", y = "No. of genes") +
            coord_flip() +
            theme_classic() +
            theme(plot.title = element_text(hjust = 0.5))

        return(p)
    }



#' @title Plot back-spliced exon positions
#'
#' @description The function plotExPosition() generates a bar chart showing
#' the position of the back-spliced exons within the transcript.
#'
#' @param annotatedBSJs A data frame with the annotated back-spliced junctions.
#' This data frame can be generated with \code{\link{annotateBSJs}}.
#'
#' @param title A character string specifying the title of the plot.
#'
#' @param n An integer specyfing the position counts cut-off. If 0 is specified
#' all position are plotted. Deafult value is 0.
#'
#' @param flip A logical specifying whether to flip the transcripts. If TRUE all
#' transcripts are flipped and the last exons will correspond to the first ones,
#' the second last exons to the second ones etc. Default value is FALSE.
#'
#' @return A ggplot object.
#'
#' @examples
#' # Load data frame containing detected back-spliced junctions
#' data("mergedBSJunctions")
#'
#' # Load short version of the gencode v19 annotation file
#' data("gtf")
#'
#' # Annotate the first 10 back-spliced junctions
#' annotatedBSJs <- annotateBSJs(mergedBSJunctions[1:10, ], gtf)
#'
#' # Plot
#' p <- plotExPosition(annotatedBSJs, title = "", n = 0, flip = FALSE)
#' p
#'
#'
#' @import ggplot2
#' @import dplyr
#' @import magrittr
#' @import reshape2
#' @importFrom rlang .data
#' @export
plotExPosition <-
    function(annotatedBSJs,
        title = "",
        n = 0,
        flip = FALSE) {
        if (flip) {
            annotatedBSJs <- annotatedBSJs %>%
                dplyr::mutate(
                    exNumUpBSE = (.data$totExons - .data$exNumUpBSE) + 1,
                    exNumDownBSE = (.data$totExons - .data$exNumDownBSE) + 1
                ) %>%
                dplyr::select(.data$id,
                    .data$gene,
                    .data$exNumUpBSE,
                    .data$exNumDownBSE)
        } else{
            annotatedBSJs <- annotatedBSJs %>%
                dplyr::select(.data$id,
                    .data$gene,
                    .data$exNumUpBSE,
                    .data$exNumDownBSE)
        }
        # Plot
        p <-  annotatedBSJs %>%
            reshape2::melt(
                id.vars = c("id", "gene"),
                variable.name = "feature",
                value.name = "exNum"
            ) %>%
            dplyr::filter(!is.na(.data$exNum)) %>%
            dplyr::group_by(.data$exNum) %>%
            dplyr::summarise(n1 = n()) %>%
            dplyr::filter(.data$n1 > n) %>%
            ggplot(aes(x = factor(.data$exNum), y = .data$n1)) +
            geom_bar(position = "dodge", stat = "identity") +
            labs(title = title, x = "Exon position", y = "Frequency") +
            theme_classic() +
            theme(
                axis.text.x = element_text(
                    angle = 90,
                    hjust = 0.5,
                    vjust = 0.5
                ),
                plot.title = element_text(hjust = 0.5)
            )
        return(p)
    }

#' @title Plot exons between back-spliced junctions
#'
#' @description The function plotExBetweenBSEs() generates a bar chart showing
#' the no. of exons in between the back-spliced junctions.
#'
#' @param annotatedBSJs A data frame with the annotated back-spliced junctions.
#' This data frame can be generated with \code{\link{annotateBSJs}}.
#'
#' @param title A character string specifying the title of the plot.
#'
#' @return A ggplot object.
#'
#' @examples
#' # Load data frame containing detected back-spliced junctions
#' data("mergedBSJunctions")
#'
#' # Load short version of the gencode v19 annotation file
#' data("gtf")
#'
#' # Annotate the first 10 back-spliced junctions
#' annotatedBSJs <- annotateBSJs(mergedBSJunctions[1:10, ], gtf)
#'
#' # Plot
#' p <- plotExBetweenBSEs(annotatedBSJs, title = "")
#' p
#'
#'
#' @import ggplot2
#' @import dplyr
#' @import magrittr
#' @importFrom rlang .data
#' @export
plotExBetweenBSEs <-
    function(annotatedBSJs, title = "") {
        # Plot
        p <-  annotatedBSJs %>%
            dplyr::filter(!is.na(.data$numOfExons)) %>%
            group_by(.data$numOfExons) %>%
            summarise(n1 = n()) %>%
            ggplot(aes(x = factor(.data$numOfExons), y = .data$n1)) +
            geom_bar(position = "dodge", stat = "identity") +
            labs(title = title, x = "No. of exons", y =
                    "Frequency") +
            theme_classic() +
            theme(
                axis.text.x = element_text(
                    angle = 90,
                    hjust = 0.5,
                    vjust = 0.5
                ),
                plot.title = element_text(hjust = 0.5)
            )

        return(p)

    }

#' @title Plot exons in the circRNA host transcript
#'
#' @description The function plotTotExons() generates a bar chart showing the
#' total number of exons (totExon column) in the transcripts selected for
#' the downstream analysis.
#'
#' @param annotatedBSJs A data frame with the annotated back-spliced junctions.
#' This data frame can be generated with \code{\link{annotateBSJs}}.
#'
#' @param title A character string specifying the title of the plot
#'
#' @return A ggplot object.
#'
#' @examples
#' # Load data frame containing detected back-spliced junctions
#' data("mergedBSJunctions")
#'
#' # Load short version of the gencode v19 annotation file
#' data("gtf")
#'
#' # Annotate the first 10 back-spliced junctions
#' annotatedBSJs <- annotateBSJs(mergedBSJunctions[1:10, ], gtf)
#'
#' # Plot
#' p <- plotTotExons(annotatedBSJs, title = "")
#' p
#'
#'
#' @import ggplot2
#' @import dplyr
#' @import magrittr
#' @importFrom rlang .data
#' @export
plotTotExons <-
    function(annotatedBSJs, title = "") {
        # Keep unique transcript
        annotatedBSJsNoDup <-
            annotatedBSJs[!duplicated(annotatedBSJs$transcript), ]

        p <- annotatedBSJsNoDup %>%
            dplyr::select(.data$totExons) %>%
            dplyr::filter(!is.na(.data$totExons)) %>%
            dplyr::group_by(.data$totExons) %>%
            dplyr::summarise(n1 = n()) %>%
            ggplot(aes(x = factor(.data$totExons), y = .data$n1)) +
            geom_bar(position = "dodge", stat = "identity") +
            labs(title = title, x = "No. of exons", y =
                    "Frequency") +
            theme_classic() +
            theme(
                axis.text.x = element_text(
                    angle = 90,
                    hjust = 0.5,
                    vjust = 0.5
                ),
                plot.title = element_text(hjust = 0.5)
            )
        return(p)
    }



#' @title Plot differential circRNA expression results
#'
#' @description The function volcanoPlot() generates a volcano plot with the
#' results of the differential expression analysis.
#'
#' @param res A data frame containing the the differential expression resuls.
#' It can be generated with \code{\link{getDeseqRes}} or
#' \code{\link{getEdgerRes}}.
#'
#' @param log2FC An integer specifying the log2FC cut-off. Deafult value is 1.
#'
#' @param padj An integer specifying the adjusted P value cut-off.
#' Deafult value is 0.05.
#'
#' @param title A character string specifying the title of the plot.
#'
#' @param gene A logical specifying whether to show all the host gene names of
#' the differentially expressed circRNAs to the plot. Deafult value is FALSE.
#'
#' @param geneSet A character vector specifying which host gene name of
#' the differentially expressed circRNAs to show in the plot. Multiple host
#' gene names can be specofied. E.g. c('TTN, RyR2')
#'
#' @param setxLim A logical specifying whether to set x scale limits.
#' If TRUE the value in xlim will be used. Deafult value is FALSE.
#'
#' @param setyLim A logical specifying whether to set y scale limits.
#' If TRUE the value in ylim will be used. Deafult value is FALSE.
#'
#' @param xlim A numeric vector specifying the lower and upper x axis limits.
#' Deafult values are c(-8 , 8).
#'
#' @param ylim An integer specifying the lower and upper y axis limits
#' Deafult values are c(0, 5).
#'
#' @param color A string specifying the color of the differentially expressed
#' circRNAs. Default value is "blue".
#'
#' @return A ggplot object.
#'
#' @examples
#' # Load data frame containing detected back-spliced junctions
#' data("mergedBSJunctions")
#'
#' pathToExperiment <- system.file("extdata", "experiment.txt",
#'     package ="circRNAprofiler")
#'
#' # Filter circRNAs
#' filterdCirc <- filterCirc(
#'     mergedBSJunctions,
#'     allSamples = FALSE,
#'     min = 5,
#'     pathToExperiment)
#'
#' # Find differentially expressed circRNAs
#'  deseqResBvsA <- getDeseqRes(
#'     filterdCirc,
#'     condition = "A-B",
#'     pAdjustMethod = "BH",
#'     pathToExperiment)
#'
#' # Plot
#' p <- volcanoPlot(
#'     deseqResBvsA,
#'     log2FC = 1,
#'     padj = 0.05,
#'     title = "",
#'     setxLim = TRUE,
#'     xlim = c(-8 , 7.5),
#'     setyLim = FALSE,
#'     ylim = c(0 , 4),
#'     gene = FALSE)
#' p
#'
#' @import ggplot2
#' @importFrom stats na.omit
#' @export
volcanoPlot <- function(res,
    log2FC = 1,
    padj = 0.05,
    title = "",
    gene = FALSE,
    geneSet = c(''),
    setxLim = FALSE,
    xlim = c(-8 , 8),
    setyLim = FALSE,
    ylim = c(0, 5),
    color = "blue") {
    res <- stats::na.omit(res)
    diffExpCirc <- stats::na.omit(res[abs(res$log2FC) >= log2FC &
            res$padj <= padj, ])

    xyLim <- .setxyLim(setxLim, xlim, setyLim, ylim, res)
    xmin <- xyLim$xmin
    xmax <- xyLim$xmax
    ymin <- xyLim$ymin
    ymax <- xyLim$ymax

    p <- ggplot(data = res, aes(x = log2FC, y = -log10(padj))) +
        geom_point(colour = "black",
            size = 3,
            na.rm = TRUE) +
        xlim(xmin, xmax) +
        ylim(ymin, ymax) +
        labs(title = title, x = "log2 FC", y = "-log10 padj") +
        geom_point(
            data = diffExpCirc,
            aes(x = log2FC, y = -log10(padj)),
            color = color,
            size = 3,
            na.rm = TRUE
        ) +
        theme(
            panel.background = element_blank(),
            plot.title = element_text(hjust = 0.5),
            panel.border = element_rect(
                colour = "black",
                fill = NA,
                size = 1.4
            )
        ) +
        geom_hline(
            yintercept = 1.3,
            linetype = "dashed",
            color = "black",
            size = 0.5
        ) +
        geom_vline(
            xintercept = -log2FC,
            linetype = "dashed",
            color = "black",
            size = 0.5
        ) +
        geom_vline(
            xintercept = log2FC,
            linetype = "dashed",
            color = "black",
            size = 0.5
        )

    if (gene) {
        p <- .addGeneName(p, diffExpCirc, log2FC, padj)
    }

    if (length(geneSet)>0){
        diffExpCircFiltered <- diffExpCirc %>%
            dplyr::filter(.data$gene %in% geneSet)
        p <- .addGeneName(p, diffExpCircFiltered, log2FC, padj)
    }

    return(p)
}

# Set xyLim
.setxyLim <- function(setxLim = FALSE,
    xlim = c(-8 , 8),
    setyLim = FALSE,
    ylim = c(0, 5),
    res) {
    if (setxLim) {
        xmin <- xlim[1]
        xmax <- xlim[2]
    } else{
        xmin <- min(res$log2FC)
        xmax <- max(res$log2FC)
    }

    if (setyLim) {
        ymin <- ylim[1]
        ymax <- ylim[2]
    } else{
        ymin <- min(-log10(res$padj))
        ymax <- max(-log10(res$padj))
    }

    xyLim <- vector("list", 4)
    names(xyLim)[1] <- "xmin"
    names(xyLim)[2] <- "xmax"
    names(xyLim)[1] <- "ymin"
    names(xyLim)[2] <- "ymax"

    xyLim$xmin <- xmin
    xyLim$xmax <- xmax

    xyLim$ymin <- ymin
    xyLim$ymax <- ymax

    return(xyLim)
}

# Add gene names to the volcano plot
.addGeneName <- function(p, diffExpCirc, log2FC, padj){
    p <- p +
        geom_text(
            data = diffExpCirc,
            aes(
                x = log2FC,
                y = -log10(padj),
                label = .data$gene
            ),
            size = 3,
            colour = "black",
            hjust = 0.5,
            vjust = -0.3
        )
    return(p)
}


#' @title Plot motifs analysis results
#'
#' @description The function plotMotifs() generates 2 bar charts showing the
#' log2FC and the number of occurences of each motif found in the target
#' sequences (e.g detected Vs randomly selected).
#'
#'
#' @param mergedMotifsFTS A data frame containing the number of occurences
#' of each motif found in foreground target sequences (e.g from detected
#' back-spliced junctions). It can be generated with the
#' \code{\link{mergeMotifs}}.
#'
#' @param mergedMotifsBTS A data frame containing the number of occurences
#' of each motif found in the background target sequences (e.g. from
#' random back-spliced junctions). It can be generated with the
#' \code{\link{mergeMotifs}}.
#'
#' @param log2FC An integer specifying the log2FC cut-off. Default value is 1.
#'
#' NOTE: log2FC is calculated as follow: normalized number of occurences of
#' each motif found in the foreground target sequences / normalized number of
#' occurences of each motif found in the background target sequences.
#'
#' To avoid infinity as a value for fold change, 1 was added to number of occurences
#' of each motif found in the foreground and background target sequences before
#' the normalization.
#'
#' @param n An integer specifying the number of motifs cutoff.
#' E.g. if 3 is specifiyed only motifs that are found at least 3 times in the
#' foreground or background target sequences are retained.
#' Deafaut value is 0.
#'
#' @param removeNegLog2FC A logical specifying whether to remove the RBPs having
#' a negative log2FC. If TRUE then only positive log2FC will be visualized.
#' Default value is FALSE.
#'
#' @param nf1 An integer specifying the normalization factor for the
#' first data frame mergedMotifsFTS. The occurrences of each motif plus 1 are
#' divided by nf1. The normalized values are then used for fold-change calculation.
#' Set this to the number or length of target sequences (e.g from detected
#' back-spliced junctions) where the motifs were extracted from.
#' Default value is 1.
#'
#' @param nf2 An integer specifying the normalization factor for the
#' second data frame mergedMotifsBTS. The occurrences of each motif plus 1 are
#' divided by nf2. The The normalized values are then used for fold-change calculation.
#' Set this to the number of target sequences (e.g from random
#' back-spliced junctions) where the motifs were extracted from.
#' Default value is 1.
#'
#' NOTE: If nf1 and nf2 is set equals to 1, the number or length of target sequences
#' (e.g detected Vs randomly selected) where the motifs were extrated from,
#' is supposed to be the same.
#'
#' @param df1Name A string specifying the name of the first data frame. This
#' will be displayed in the legend of the plot. Deafult value is "foreground".
#'
#' @param df2Name A string specifying the name of the first data frame. This
#' will be displayed in the legend of the plot. Deafult value is "background".
#'
#' @param angle An integer specifying the rotation angle of the axis labels.
#' Default value is 0.
#'
#' @return A ggplot object.
#'
#' @examples
#' # Load data frame containing detected back-spliced junctions
#' data("mergedBSJunctions")
#'
#' # Load short version of the gencode v19 annotation file
#' data("gtf")
#'
#' # Annotate the first back-spliced junctions
#' annotatedFBSJs <- annotateBSJs(mergedBSJunctions[1, ], gtf)
#'
#' # Get random back-spliced junctions
#' randomBSJunctions <- getRandomBSJunctions(gtf, n = 1, f = 10)
#'
#' # Annotate random back-spliced junctions
#' annotatedBBSJs <- annotateBSJs(randomBSJunctions, gtf, isRandom = TRUE)
#'
#' # Get genome
#' genome <- BSgenome::getBSgenome("BSgenome.Hsapiens.UCSC.hg19")
#'
#' # Retrieve target sequences from detected back-spliced junctions
#' targetsFTS <- getSeqsFromGRs(
#'    annotatedFBSJs,
#'    genome,
#'    lIntron = 200,
#'    lExon = 10,
#'    type = "ie"
#'    )
#'
#' # Retrieve target sequences from random back-spliced junctions
#' targetsBTS <- getSeqsFromGRs(
#'    annotatedBBSJs,
#'    genome,
#'    lIntron = 200,
#'    lExon = 10,
#'    type = "ie"
#'    )
#'
#' # Get motifs
#'  motifsFTS <- getMotifs(
#'      targetsFTS,
#'      width = 6,
#'      database = 'ATtRACT',
#'      species = "Hsapiens",
#'      rbp = TRUE,
#'      reverse = FALSE)
#'
#' motifsBTS <- getMotifs(
#'      targetsBTS,
#'      width = 6,
#'      database = 'ATtRACT',
#'      species = "Hsapiens",
#'      rbp = TRUE,
#'      reverse = FALSE)
#'
#' # Merge motifs
#'  mergedMotifsFTS <- mergeMotifs(motifsFTS)
#'  mergedMotifsBTS <- mergeMotifs(motifsBTS)
#'
#' # Plot
#'  p <- plotMotifs(
#'      mergedMotifsFTS,
#'      mergedMotifsBTS,
#'      log2FC = 2,
#'      nf1 = nrow(annotatedFBSJs),
#'      nf2 = nrow(annotatedBBSJs),
#'      df1Name = "foreground",
#'      df2Name = "background")
#'
#' @importFrom rlang .data
#' @import ggplot2
#' @import dplyr
#' @import magrittr
#' @export
plotMotifs <-
    function(mergedMotifsFTS,
        mergedMotifsBTS,
        log2FC = 1,
        n = 0,
        removeNegLog2FC = FALSE,
        nf1 = 1,
        nf2 = 1,
        df1Name = "foreground",
        df2Name = "background",
        angle = 0) {
        mergedMotifsAll <-
            .getMergedMotifsAll(
                mergedMotifsFTS,
                mergedMotifsBTS,
                log2FC,
                n,
                nf1,
                nf2
            )

        if(removeNegLog2FC){
            mergedMotifsAll<- mergedMotifsAll %>%
                dplyr::filter(log2FC >=0)
        }

        p <- list()
        p[[1]] <-  mergedMotifsAll %>%
            ggplot(aes(x = .data$id,
                y = .data$log2FC)) +
            geom_bar(position = "dodge", stat = "identity") +
            labs(title = "", x = "id", y = "log2 FC") +
            coord_flip() +
            theme_classic() +
            theme(plot.title = element_text(angle = 90, hjust = 0.5))


        p[[2]] <- mergedMotifsAll %>%
            dplyr::select(.data$id, .data$foregroundNorm, .data$backgroundNorm) %>%
            reshape2::melt(
                id.vars = "id",
                variable.name = "circRNA",
                value.name = "count"
            ) %>%
            ggplot(aes(
                x = .data$id,
                y = .data$count,
                fill = .data$circRNA
            )) +
            geom_bar(position = "dodge", stat = "identity") +
            labs(title = "", x = "id", y = "normalized count") +
            coord_flip() +
            theme_classic() +
            theme(plot.title = element_text(angle = 90, hjust = 0.5),
                  axis.text.x = element_text(angle = angle, hjust = 1)) +
            scale_fill_grey(name = "circRNA", labels = c(df1Name, df2Name))

        p[[3]] <- mergedMotifsAll %>%
            dplyr::arrange(desc(.data$log2FC))
        return(p)
    }


# get mergedMotifsAll data frame
.getMergedMotifsAll <- function(mergedMotifsFTS,
    mergedMotifsBTS,
    log2FC = 1,
    n = 0,
    nf1 = 1,
    nf2 = 1) {
    mergedMotifsAll <-
        base::merge(mergedMotifsFTS,
            mergedMotifsBTS,
            by = "id",
            all = TRUE) %>%
        dplyr::rename(foreground = .data$count.x,
            background = .data$count.y) %>%
        dplyr::mutate(
            foreground = ifelse(is.na(.data$foreground), 0, .data$foreground),
            background = ifelse(is.na(.data$background), 0, .data$background)
        ) %>%
        dplyr::filter(.data$foreground >= n | .data$background >= n )%>%
        dplyr::mutate(
            foregroundNorm = (.data$foreground+1) / nf1,
            backgroundNorm = (.data$background+1) / nf2
        ) %>%
        dplyr::mutate(log2fc = log2(.data$foregroundNorm / .data$backgroundNorm)) %>%
        dplyr::filter(abs(.data$log2fc) >= log2FC) %>%
        dplyr::arrange(.data$log2fc) %>%
        dplyr::rename(log2FC = .data$log2fc) %>%
        dplyr::mutate(id = factor(.data$id, levels = .data$id)) %>%
        dplyr::mutate(motif.x = ifelse(is.na(.data$motif.x), .data$motif.y, .data$motif.x)) %>%
        dplyr::rename(motifF = .data$motif.x) %>%
        dplyr::rename(motifB = .data$motif.y) %>%
        dplyr::select(
            .data$id,
            .data$foreground,
            .data$background,
            .data$foregroundNorm,
            .data$backgroundNorm,
            .data$log2FC,
            .data$motifF,
            .data$motifB
        )
    return(mergedMotifsAll)
}


#' @title Plot miRNA analysis results
#'
#' @description The function plotMiR() generates a scatter plot showing the
#' number of miRNA binding sites for each miR found in the target sequence.
#'
#' @param rearragedMiRres A list containing containing rearranged
#' miRNA analysis results. See \code{\link{getMiRsites}} and then
#' \code{\link{rearrangeMiRres}}.
#'
#' @param n An integer specifying the miRNA binding sites cut-off. The miRNA
#' with a number of binding sites equals or higher to the cut-off will be
#' colored. Deafaut value is 40.
#'
#' @param color A string specifying the color of the top n miRs. Default value
#' is "blue".
#'
#' @param  miRid A logical specifying whether or not to show the miR ids in
#' the plot. default value is FALSE.
#'
#' @param id An integer specifying which element of the list
#' rearragedMiRres to plot. Each element of the list contains
#' the miR resutls relative to one circRNA. Deafult value is 1.
#'
#' @return A ggplot object.
#'
#'
#' @examples
#' # Load data frame containing detected back-spliced junctions
#' data("mergedBSJunctions")
#'
#' # Load short version of the gencode v19 annotation file
#' data("gtf")
#'
#' # Annotate the first back-spliced junctions
#' annotatedBSJs <- annotateBSJs(mergedBSJunctions[1, ], gtf)
#'
#' # Get genome
#' genome <- BSgenome::getBSgenome("BSgenome.Hsapiens.UCSC.hg19")
#'
#' # Retrieve target sequences.
#' targets <- getCircSeqs(
#'     annotatedBSJs,
#'     gtf,
#'     genome)
#'
#' # Screen target sequence for miR binding sites.
#' pathToMiRs <- system.file("extdata", "miRs.txt", package="circRNAprofiler")
#'
#' # miRsites <- getMiRsites(
#' #     targets,
#' #     miRspeciesCode = "hsa",
#' #     miRBaseLatestRelease = TRUE,
#' #     totalMatches = 6,
#' #     maxNonCanonicalMatches = 1,
#' #     pathToMiRs)
#'
#' # Rearrange miR results
#' # rearragedMiRres <- rearrangeMiRres(miRsites)
#'
#' # Plot
#' # p <- plotMiR(
#' #     rearragedMiRres,
#' #     n = 20,
#' #     color = "blue",
#' #     miRid = TRUE,
#' #     id = 3)
#' # p
#'
#' @importFrom rlang .data
#' @import ggplot2
#' @import dplyr
#' @import magrittr
#' @export
plotMiR <-
    function(rearragedMiRres,
        n = 40,
        color = "blue",
        miRid = FALSE,
        id = 1) {
        topMir <- .getTopMir(rearragedMiRres, id, n)
        p <-  rearragedMiRres[[id]][[2]] %>%
            dplyr::mutate(miRid = stringr::str_replace(.data$miRid, ">", "")) %>%
            dplyr::filter(!is.na(.data$counts)) %>%
            dplyr::arrange(counts) %>%
            ggplot(aes(x = .data$miRid, y = .data$counts)) +
            geom_point(colour = "black", size = 4) +
            labs(title = "", x = "miRNA", y = "No. of binding sites") +
            geom_point(
                data = topMir,
                aes(x = .data$miRid, y = .data$counts),
                color = color,
                size = 4
            ) +
            theme(
                panel.background = element_blank(),
                plot.title = element_text(hjust = 0.5),
                axis.text.x = element_blank(),
                # element_text(angle = 90, hjust = 1),
                axis.ticks.x = element_blank(),
                panel.border = element_rect(
                    colour = "black",
                    fill = NA,
                    size = 1.4
                )
            ) +
            expand_limits(y = 0)
        if (miRid) {
            p <- .addMiRid(p, topMir )
        }
        return(p)
    }


# add miR ids to the miR plot
.addMiRid <- function(p, topMir){
    p <- p +
        geom_text(
            data = topMir,
            aes(
                x = .data$miRid,
                y = .data$counts,
                label = .data$miRid
            ),
            size = 4,
            colour = "black",
            hjust = 0.5,
            vjust = -0.3
            #angle = 90
        )
    return(p)
}



# Get miRNAs with a number of binding sites higher than n
.getTopMir <- function(rearragedMiRres,
    id = 1,
    n = 40) {
    topMir <- rearragedMiRres[[id]][[2]] %>%
        dplyr::filter(.data$counts >= n) %>%
        dplyr::mutate(miRid = stringr::str_replace(.data$miRid, ">", ""))
    return(topMir)
}

# If the function you are looking for is not here check supportFunction.R
# Functions in supportFunction.R are used by multiple functions.
