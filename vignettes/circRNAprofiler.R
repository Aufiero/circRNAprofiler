## ---- global_options, include = FALSE------------------------------------
library(knitr)
library(citr)
knitr::opts_chunk$set(fig.path='figs/', warning=FALSE, message=FALSE, collapse=TRUE)
source("render_toc.R")

## ---- toc, echo = FALSE--------------------------------------------------
render_toc("circRNAprofiler.Rmd")

## ---- echo=FALSE, out.width='90%', fig.align="center", fig.cap="\\label{fig:figs} Schematic representation of the circRNA analysis workflow implemented by circRNAprofiler. In the grey box are reported the modules with the corresponding main functions."----

knitr::include_graphics("./images/image1.png")

## ---- eval = FALSE-------------------------------------------------------
#  install.packages("devtools")
#  library(devtools)
#  install_github("Aufiero/circRNAprofiler")

## ------------------------------------------------------------------------
library(circRNAprofiler)

# Packages needed for the vignettes
library(ggpubr)
library(ggplot2)
library(gridExtra)

## ---- echo=FALSE, out.width='40%', fig.align="center", fig.cap="\\label{fig:figs} Project folder structure"----
knitr::include_graphics("./images/image2.png")

## ---- eval = FALSE-------------------------------------------------------
#  initCircRNAprofiler(projectFolderName = "circRNAprofiler", predictionTools =
#                        "mapsplice")

## ---- eval = FALSE-------------------------------------------------------
#  initCircRNAprofiler(
#      projectFolderName = "circRNAprofiler",
#      predictionTools = c("mapsplice", "nclscan", "circmarker")
#  )

## ---- echo=FALSE---------------------------------------------------------
experiment <-
    read.table(
        "experiment.txt",
        header = TRUE,
        stringsAsFactors = FALSE,
        sep = "\t"
    )
knitr::kable(experiment)

## ---- echo=FALSE---------------------------------------------------------
motifs <-
    read.table(
        "motifs.txt",
        stringsAsFactors = FALSE,
        header = TRUE,
        sep = "\t"
    )
knitr::kable(motifs)

## ---- echo=FALSE---------------------------------------------------------
traits <-
    read.table(
        "traits.txt",
        stringsAsFactors = FALSE,
        header = TRUE,
        sep = "\t"
    )
knitr::kable(head(traits))

## ---- echo=FALSE---------------------------------------------------------
miRs <-
    read.table(
        "miRs.txt",
        header = TRUE,
        stringsAsFactors = FALSE,
        sep = "\t"
    )
knitr::kable(head(miRs))

## ---- echo=FALSE---------------------------------------------------------
transcripts <-
    read.table(
        "transcripts.txt",
        header = TRUE,
        stringsAsFactors = FALSE,
        sep = "\t"
    )
knitr::kable(transcripts)

## ---- echo=FALSE---------------------------------------------------------
circRNApredictions <-
    read.table(
        "circRNAs_test.txt",
        header = TRUE,
        stringsAsFactors = FALSE,
        sep = "\t"
    )
knitr::kable(circRNApredictions)

## ---- eval = FALSE-------------------------------------------------------
#  check <- checkProjectFolder()
#  check

## ------------------------------------------------------------------------
# For example purpose load a short version of the already formatted annotation
# file gencode.V19.annotation.gtf (downloaded from https://www.gencodegenes.org/)
data("gtf")
head(gtf)

# Alternatively put the gtf file in the project folder, then run:
# gtf <- formatGTF(pathToGTF = "gencode.V19.annotation.gtf")

## ---- eval= FALSE--------------------------------------------------------
#  # Load the object containing the detected circRNAs
#  data("backSplicedJunctions")
#  head(backSplicedJunctions)
#  
#  # Alternatively run:
#  # backSplicedJunctions <- getBackSplicedJunctions(gtf)

## ---- fig.align="center", fig.width = 10, fig.height = 3-----------------
# Plot
p <- ggplot(backSplicedJunctions, aes(x = tool)) +
    geom_bar() +
    labs(title = "", x = "Detection tool", y = "No. of circRNAs") +
    theme_classic()

dt <- getDetectionTools() %>%
    dplyr::filter( name %in% c("mapsplice","nclscan", "circmarker"))%>%
    gridExtra::tableGrob(rows=NULL)

# Merge plots
gridExtra::grid.arrange(p, dt, nrow=1)

## ------------------------------------------------------------------------
# Load object containing the merged back-spliced junctions
data("mergedBSJunctions")
head(mergedBSJunctions)

# Alternatively run:
# mergedBSJunctions <-
# mergeBSJunctions(backSplicedJunctions, gtf)

## ---- fig.align = "center", fig.width = 10, fig.height = 4---------------
# Plot
p <- ggplot(mergedBSJunctions, aes(x = tool)) +
    geom_bar() +
    labs(title = "", x = "Detection tool", y = "No. of circRNAs") +
    theme_classic()

# Run getDetectionTools() to get the code corresponding to the circRNA 
# detection tools.
dt <-getDetectionTools() %>%
    dplyr::filter( name %in% c("mapsplice","nclscan", "circmarker"))%>%
    gridExtra::tableGrob(rows=NULL)

# Merge plots
gridExtra::grid.arrange(p, dt, nrow=1)


## ------------------------------------------------------------------------
filteredCirc <-
filterCirc(mergedBSJunctions, allSamples = FALSE, min = 5)

## ---- fig.align="center", fig.width = 10, fig.height = 4-----------------
# Plot
p <- ggplot(filteredCirc, aes(x = tool)) +
    geom_bar() +
    labs(title = "", x = "Detection tool", y = "No. of circRNAs") +
    theme_classic()

# Run getDetectionTools() to get the code corresponding to the circRNA
# detection tools.
dt <-getDetectionTools() %>%
    dplyr::filter( name %in% c("mapsplice","nclscan", "circmarker"))%>%
    gridExtra::tableGrob(rows=NULL)

# Merge plots
gridExtra::grid.arrange(p, dt, nrow=1)

## ------------------------------------------------------------------------
# Compare condition B Vs A
deseqResBvsA <-
    getDeseqRes(
        filteredCirc,
        condition = "A-B",
        fitType = "local",
        pAdjustMethod = "BH"
    )
head(deseqResBvsA)

## ------------------------------------------------------------------------
# Compare condition C Vs A
deseqResCvsA <-
    getDeseqRes(
        filteredCirc,
        condition = "A-C",
        fitType = "local",
        pAdjustMethod = "BH"
    )
head(deseqResCvsA)

## ---- fig.align="center", fig.height= 8, fig.width = 8-------------------
# We set the xlim and ylim to the same values for both plots to make them
# comparable. Before setting the axis limits, you should visualize the 
# plots with the default values to be able to define the correct limits 
p1 <-
    volcanoPlot(
        deseqResBvsA,
        log2FC = 1,
        padj = 0.05,
        title = "DCMs Vs. Con",
        setxLim = TRUE,
        xlim = c(-8 , 7.5),
        setyLim = FALSE,
        ylim = c(0 , 4),
        gene = FALSE
    )
p2 <-
    volcanoPlot(
        deseqResCvsA,
        log2FC = 1,
        padj = 0.05,
        title = "HCMs Vs. Con",
        setxLim = TRUE,
        xlim = c(-8 , 7.5),
        setyLim = TRUE,
        ylim = c(0 , 4),
        gene = FALSE
    )
ggarrange(p1, 
          p2, 
          ncol = 1, 
          nrow = 2)

## ----eval = FALSE--------------------------------------------------------
#  # Compare condition B Vs A
#  edgerResBvsA <-
#      getEdgerRes(
#          filteredCirc,
#          condition = "A-B",
#          normMethod = "TMM",
#          pAdjustMethod = "BH"
#      )
#  head(edgerResBvsA)

## ----eval = FALSE--------------------------------------------------------
#  # Compare condition C Vs A
#  edgerResCvsA <-
#      getEdgerRes(
#          filteredCirc,
#          condition = "A-C",
#          normMethod = "TMM",
#          pAdjustMethod = "BH"
#      )
#  head(edgerResCvsA)

## ---- eval = FALSE-------------------------------------------------------
#  liftedBSJCoords <- liftBSJCoords(filteredCirc, map = "hg19ToMm9",
#                                   annotationHubID = "AH14155")

## ------------------------------------------------------------------------
# As an example of the 1458 filtered circRNAs we annotate only the firt 30 
# circRNAs
annotatedBSJs <- annotateBSJs(filteredCirc[1:30,], gtf) 
head(annotatedBSJs)

## ------------------------------------------------------------------------
# First find frequency of single exon circRNAs
f <-
    sum((annotatedBSJs$exNumUpBSE == 1 |
             annotatedBSJs$exNumDownBSE == 1) ,
        na.rm = TRUE) / (nrow(annotatedBSJs) * 2)

# Retrieve random back-spliced junctions
randomBSJunctions <-
    getRandomBSJunctions(n = nrow(annotatedBSJs), f = f, gtf)
head(randomBSJunctions)

## ---- eval = FALSE-------------------------------------------------------
#  annotatedRBSJs <- annotateBSJs(randomBSJunctions, gtf, isRandom = TRUE)

## ---- fig.align="center", fig.width = 10, fig.height = 7, eval = FALSE----
#  # annotatedBSJs act as foreground data set
#  # annotatedRBSJs act as background data set
#  
#  # Length of flanking introns
#  p1 <- plotLenIntrons(
#      annotatedBSJs,
#      annotatedRBSJs,
#      title = "Length flanking introns",
#      df1Name = "detected",
#      df2Name = "random"
#  )
#  
#  # Length of back-splided exons
#  p2 <- plotLenBSEs(
#      annotatedBSJs,
#      annotatedRBSJs,
#      title = "Length back-splided exons",
#      df1Name = "detected",
#      df2Name = "random"
#  )
#  
#  # No. of circRNAs produced from the host genes
#  p3 <-
#      plotHostGenes(annotatedBSJs, title = "# CircRNAs produced from the host genes")
#  
#  # No. of exons in between the back-spliced junctions
#  p4 <-
#      plotExBetweenBSEs(annotatedBSJs, title = "# Exons between back-spliced junctions")
#  
#  # Position of back-spliced exons within the host transcripts
#  p5 <-
#      plotExPosition(annotatedBSJs,
#                     flip = FALSE,
#                     n = 1,
#                     title = "Position back-spliced exons in the transcripts")
#  
#  # Additional plotting functions
#  # Position of back-spliced exons within the flipped host transcripts
#  p6 <-
#      plotExPosition(annotatedBSJs,
#                     flip = TRUE,
#                     n = 1,
#                     title = "Position back-spliced exons in the flipped transcripts")
#  
#  # Total no. of exons within the host transcripts
#  p7 <-
#      plotTotExons(annotatedBSJs, title = " Total number of exons in the transcripts")
#  
#  # Combine plots
#  ggarrange(p1,
#            p2,
#            p3,
#            p4,
#            p5,
#            ncol = 2,
#            nrow = 3)
#  

## ---- echo=FALSE, out.width='100%', fig.align="center", fig.cap="\\label{fig:figs} Comparison of structural features extracted from the subset of 1458 filtered  back-spliced junctions compared to an equal number of randomly generated back-spliced junctions."----
knitr::include_graphics("./images/image3.png")

## ---- eval = FALSE-------------------------------------------------------
#  # Select ALPK2:-:chr18:56247780:56246046 circRNA
#  annotatedCirc <-
#  annotatedBSJs[annotatedBSJs$id == "ALPK2:-:chr18:56247780:56246046", ]
#  
#  # As background data set we used all the remaining 1457 filered circRNAs.
#  # Alternatively the subset of randomly generated back-spliced junctions can be used.
#  annotatedBackgroundCircs <-
#  annotatedBSJs[which(annotatedBSJs$id != "ALPK2:-:chr18:56247780:56246046"), ]

## ----eval = FALSE--------------------------------------------------------
#  # Foreground target sequences
#  targetsFTS_circ <-
#      getCircSeqs(annotatedCirc, gtf, species = "Hsapiens", genome = "hg19")
#  

## ---- eval = FALSE-------------------------------------------------------
#  # Foreground target sequences
#  targetsFTS_bsj <-
#      getSeqsAcrossBSJs(annotatedCirc, gtf, species = "Hsapiens", genome = "hg19")
#  

## ---- eval = FALSE-------------------------------------------------------
#  # Foreground target sequences
#  targetsFTS_gr <-
#      getSeqsFromGRs(
#          annotatedCirc,
#          lIntron = 200,
#          lExon = 9,
#          type = "ie",
#          species = "Hsapiens",
#          genome = "hg19"
#      )
#  # Background target sequences.
#  targetsBTS_gr <-
#      getSeqsFromGRs(
#          annotatedBackgroundCircs,
#          lIntron = 200,
#          lExon = 9,
#          type = "ie",
#          species = "Hsapiens",
#          genome = "hg19"
#      )

## ---- eval = FALSE-------------------------------------------------------
#  # Find motifs in the foreground target sequences
#  motifsFTS_gr <-
#      getMotifs(targetsFTS_gr,
#                species = "Hsapiens",
#                width = 6,
#                rbp = TRUE,
#                reverse = FALSE)
#  # Find motifs in the background target sequences
#  motifsBTS_gr <-
#      getMotifs(targetsBTS_gr,
#                species = "Hsapiens",
#                width = 6,
#                rbp = TRUE,
#                reverse = FALSE)

## ---- eval = FALSE-------------------------------------------------------
#  mergedMotifsFTS_gr <- mergeMotifs(motifsFTS_gr)
#  mergedMotifsBTS_gr <- mergeMotifs(motifsBTS_gr)

## ---- fig.align="center", fig.width = 7, fig.height = 7, eval = FALSE----
#  # Plot
#  p <-
#      plotMotifs(
#          mergedMotifsFTS_gr,
#          mergedMotifsBTS_gr,
#          nf1 = nrow(annotatedCirc),
#          nf2 = nrow(annotatedBackgroundCircs),
#          log2FC = 1,
#          df1Name = "circALPK2",
#          df2Name = "Other circRNAs"
#      )
#  ggarrange(p[[1]],
#            p[[2]],
#            labels = c("", ""),
#            ncol = 2,
#            nrow = 1)
#  

## ---- echo=FALSE, out.width='80%', fig.align="center", fig.cap="\\label{fig:figs} Bar chart showing the log2FC (cut-off = 1) and the counts of the RBP motifs found in the region flanking the predicted back-spliced junctions of circALPK2 compared to the remaining 1457 detected circRNAs."----
knitr::include_graphics("./images/image4.png")

## ---- eval = FALSE-------------------------------------------------------
#  # Type p[[3]] to get table
#  head(p[[3]])

## ----eval = FALSE--------------------------------------------------------
#  miRsites <-
#      getMiRsites(
#          targetsFTS_circ,
#          species = "Hsapiens",
#          genome = "hg19",
#          miRspeciesCode = "hsa",
#          miRBaseLatestRelease = TRUE,
#          totalMatches = 6,
#          maxNonCanonicalMatches = 1
#      )
#  

## ----eval = FALSE--------------------------------------------------------
#  rearragedMiRres <-
#      rearrangeMiRres(
#          miRsites,
#          n = 40,
#          color = "blue",
#          miRid = TRUE,
#          id = 1
#      )

## ---- eval=FALSE---------------------------------------------------------
#  # If multiple circRNAs have been anlyzed for the presence of miR binding sites
#  # the following code can store the predictions for each circRNA in a
#  # different sheet of an xls file for a better user consulation.
#  i <- 1
#  j <- 1
#  while (i <= (length(rearragedMiRres))) {
#      write.xlsx2(
#          rearragedMiRres[[i]][[1]],
#          "miRsites_TM7_NCM0.xlsx",
#          paste("sheet", j, sep = ""),
#          append = TRUE
#      )
#      j <- j + 1
#      write.xlsx2(
#          rearragedMiRres[[i]][[2]],
#          "miRsites_TM7_NCM0.xlsx",
#          paste("sheet", j, sep = ""),
#          append = TRUE
#      )
#      i <- i + 1
#      j <- j + 1
#  }
#  

## ---- eval = FALSE-------------------------------------------------------
#  p <- plotMiR(rearragedMiRres,
#               n = 40,
#               color = "blue",
#               miRid = TRUE,
#               id = 1)
#  p

## ---- eval = FALSE-------------------------------------------------------
#  snpsGWAS <- annotateSNPsGWAS(targetsFTS_gr, genome = "hg19")
#  

## ---- eval = FALSE-------------------------------------------------------
#  repeats <-
#      annotateRepeats(targetsFTS_gr, annotationHubID = "AH5122",
#                      complementary = TRUE)

## ------------------------------------------------------------------------
sessionInfo()

