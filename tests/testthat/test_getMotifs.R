setwd(paste(getwd(), "testdata", sep = "/"))

context("Test that getMotifs() function works correctly")
test_that("getMotifs() and mergeMotifs() generate the correct data structure
    with GR seqs",
    {
        gtf <- formatGTF(pathToGTF = "gencodeVM16.gtf")

        # Create the backSplicedJunctions without retriving the missing coordinates
        backSplicedJunctions <- getBackSplicedJunctions(gtf)
        mergedBSJunctions <-
            mergeBSJunctions(backSplicedJunctions, gtf)

        # Annonate BSJs
        annotatedBSJs <- annotateBSJs(mergedBSJunctions, gtf)

        if (!requireNamespace("BSgenome.Mmusculus.UCSC.mm10", quietly = TRUE))
            install.packages("BSgenome.Mmusculus.UCSC.mm10")

        # Get BSgenome object
        genome <- BSgenome::getBSgenome("BSgenome.Mmusculus.UCSC.mm10")

        # Retrieve sequences
        targets <-
            getSeqsFromGRs(
                annotatedBSJs,
                genome,
                lIntron = 101,
                lExon = 9,
                type = "ie"
            )

        # Retrieve motifs
        # rbp = TRUE reverse = TRUE
        motifs <-
            getMotifs(
                targets,
                width = 6,
                species = "Mmusculus",
                rbp = TRUE,
                reverse = TRUE
            )
        expect_is(motifs, "list")
        expect_equal(length(motifs), 2)

        # check nrow all data frame
        expect_equal(nrow(motifs$upGR$targets), nrow(annotatedBSJs))
        expect_equal(nrow(motifs$downGR$targets), nrow(annotatedBSJs))
        expect_equal(ncol(motifs$upGR$counts) - 1, nrow(motifs$upGR$motif))

        # Merge motifs
        mergedMotifs <- mergeMotifs(motifs)
        expect_is(mergedMotifs, "data.frame")

        # Retrieve motifs
        # rbp = FALSE reverse = FALSE
        motifs <-
            getMotifs(
                targets,
                width = 6,
                species = "Mmusculus",
                rbp = FALSE,
                reverse = FALSE
            )
        expect_is(motifs, "list")
        expect_equal(length(motifs), 2)

        # check nrow all data frame
        expect_equal(nrow(motifs$upGR$targets), nrow(annotatedBSJs))
        expect_equal(nrow(motifs$downGR$targets), nrow(annotatedBSJs))
        expect_equal(ncol(motifs$upGR$counts) - 1, nrow(motifs$upGR$motif))

        # Merge motifs
        mergedMotifs <- mergeMotifs(motifs)
        expect_is(mergedMotifs, "data.frame")

    })


test_that("getMotifs() generates a list with the correct
    content with GR seqs",
    {
        gtf <- formatGTF(pathToGTF = "gencodeVM16.gtf")

        # Create the backSplicedJunctions without retriving the missing coordinates
        backSplicedJunctions <- getBackSplicedJunctions(gtf)
        mergedBSJunctions <-
            mergeBSJunctions(backSplicedJunctions, gtf)

        # Annonate BSJs
        annotatedBSJs <- annotateBSJs(mergedBSJunctions, gtf)

        if (!requireNamespace("BSgenome.Mmusculus.UCSC.mm10", quietly = TRUE))
            install.packages("BSgenome.Mmusculus.UCSC.mm10")

        # Get BSgenome object
        genome <- BSgenome::getBSgenome("BSgenome.Mmusculus.UCSC.mm10")

        # Retrieve sequences
        targets <-
            getSeqsFromGRs(
                annotatedBSJs,
                genome,
                lIntron = 101,
                lExon = 9,
                type = "ie"
            )

        # Retrieve motifs
        motifs <-
            getMotifs(
                targets,
                width = 4,
                species = "Mmusculus",
                rbp = TRUE,
                reverse = TRUE
            )

        expect_equal(motifs$upGR$counts$UCUU[4], 0)
        expect_equal(motifs$downGR$counts$UCUU[4], 1)

        expect_equal(motifs$upGR$counts$UUCU[4], 0)
        expect_equal(motifs$downGR$counts$UUCU[4], 2)


        expect_equal(motifs$upGR$location$UUCU[4], "NA")
        expect_equal(motifs$downGR$location$UUCU[4], "26,102")

        expect_equal(motifs$upGR$location$UCUU[4], "NA")
        expect_equal(motifs$downGR$location$UCUU[4], "103")

        # check 2nd pattern
        expect_equal(motifs$upGR$counts$AGAG[4], 2)
        expect_equal(motifs$downGR$counts$AGAG[4], 2)

        expect_equal(motifs$upGR$counts$GAGA[4], 1)
        expect_equal(motifs$downGR$counts$GAGA[4], 2)

        expect_equal(motifs$upGR$location$GAGA[4], "13")
        expect_equal(motifs$downGR$location$GAGA[4], "18,56")

        expect_equal(motifs$upGR$location$AGAG[4], "12,14")
        expect_equal(motifs$downGR$location$AGAG[4], "17,55")

    })




test_that(
    "getMotifs() and mergeMotifs() generates the correct data structure
    with circRNA and BSJ seqs",
    {
        gtf <- formatGTF(pathToGTF = "gencodeVM16.gtf")

        # Create the backSplicedJunctions without retriving the missing coordinates
        backSplicedJunctions <- getBackSplicedJunctions(gtf)
        mergedBSJunctions <-
            mergeBSJunctions(backSplicedJunctions, gtf)

        # Annonate BSJs
        annotatedBSJs <-
            annotateBSJs(mergedBSJunctions[1:3, ], gtf)

        if (!requireNamespace("BSgenome.Mmusculus.UCSC.mm10", quietly = TRUE))
            install.packages("BSgenome.Mmusculus.UCSC.mm10")

        # Get BSgenome object
        genome <- BSgenome::getBSgenome("BSgenome.Mmusculus.UCSC.mm10")

        # Retrieve sequences (BSJs)
        targets <-
            getSeqsAcrossBSJs(annotatedBSJs,
                gtf,
                genome
                )

        # Retrieve motifs
        # rbp = TRUE reverse = TRUE
        motifs <-
            getMotifs(
                targets,
                width = 6,
                species = "Mmusculus",
                rbp = TRUE,
                reverse = TRUE
            )
        expect_is(motifs, "list")
        expect_equal(length(motifs), 1)

        # check nrow all data frame
        expect_equal(nrow(motifs$bsj$targets), nrow(annotatedBSJs))
        expect_equal(ncol(motifs$bsj$counts) - 1, nrow(motifs$bsj$motif))

        # Merge motifs
        mergedMotifs <- mergeMotifs(motifs)
        expect_is(mergedMotifs, "data.frame")

        # Retrieve motifs
        # rbp = FALSE reverse = FALSE
        motifs <-
            getMotifs(
                targets,
                width = 6,
                species = "Mmusculus",
                rbp = FALSE,
                reverse = FALSE
            )
        expect_is(motifs, "list")
        expect_equal(length(motifs), 1)

        # check nrow all data frame
        expect_equal(nrow(motifs$bsj$targets), nrow(annotatedBSJs))
        expect_equal(ncol(motifs$bsj$counts) - 1, nrow(motifs$bsj$motif))

        # Merge motifs
        mergedMotifs <- mergeMotifs(motifs)
        expect_is(mergedMotifs, "data.frame")

        # Get BSgenome object
        genome <- BSgenome::getBSgenome("BSgenome.Mmusculus.UCSC.mm10")

        # Retrieve sequences (circRNA seqs)
        targets <-
            getCircSeqs(annotatedBSJs,
                gtf,
                genome
                )

        # Retrieve motifs
        # rbp = TRUE reverse = TRUE
        motifs <-
            getMotifs(
                targets,
                width = 6,
                species = "Mmusculus",
                rbp = TRUE,
                reverse = TRUE
            )

        expect_is(motifs, "list")
        expect_equal(length(motifs), 1)

        # check nrow all data frame
        expect_equal(nrow(motifs$circ$targets), nrow(annotatedBSJs))
        expect_equal(ncol(motifs$circ$counts) - 1, nrow(motifs$circ$motif))

        # Merge motifs
        mergedMotifs <- mergeMotifs(motifs)
        expect_is(mergedMotifs, "data.frame")

        # Retrieve motifs
        # rbp = FALSE reverse = FALSE
        motifs <-
            getMotifs(
                targets,
                width = 6,
                species = "Mmusculus",
                rbp = FALSE,
                reverse = FALSE
            )

        expect_is(motifs, "list")
        expect_equal(length(motifs), 1)

        # check nrow all data frame
        expect_equal(nrow(motifs$circ$targets), nrow(annotatedBSJs))
        expect_equal(ncol(motifs$circ$counts) - 1, nrow(motifs$circ$motif))

        # Merge motifs
        mergedMotifs <- mergeMotifs(motifs)
        expect_is(mergedMotifs, "data.frame")

    }
)


test_that("plotMotifs() generates the correct data structure",
    {
        gtf <- formatGTF(pathToGTF = "gencodeVM16.gtf")

        # Create the backSplicedJunctions without retriving the missing coordinates
        backSplicedJunctions <- getBackSplicedJunctions(gtf)
        mergedBSJunctions <-
            mergeBSJunctions(backSplicedJunctions, gtf)

        # Annonate BSJs
        annotatedBSJs <-
            annotateBSJs(mergedBSJunctions[1,], gtf)

        if (!requireNamespace("BSgenome.Mmusculus.UCSC.mm10", quietly = TRUE))
            install.packages("BSgenome.Mmusculus.UCSC.mm10")

        # Get BSgenome object
        genome <- BSgenome::getBSgenome("BSgenome.Mmusculus.UCSC.mm10")

        # Retrieve sequences
        targets <-
            getSeqsFromGRs(
                annotatedBSJs,
                genome,
                lIntron = 50,
                lExon = 10,
                type = "ie"
            )

        # Retrieve motifs
        # rbp = TRUE reverse = TRUE
        motifs <-
            getMotifs(
                targets,
                width = 6,
                species = "Mmusculus",
                rbp = TRUE,
                reverse = TRUE
            )


        randomBSJunctions <-
            getRandomBSJunctions (gtf, n = nrow(annotatedBSJs), f = 10)

        # Annonate BSJs
        annotatedRBSJs <-
            annotateBSJs(randomBSJunctions, gtf)

        # Retrieve sequences
        targetsRBSJs <-
            getSeqsFromGRs(
                annotatedRBSJs,
                genome,
                lIntron = 50,
                lExon = 10,
                type = "ie"
            )

        # Retrieve motifs
        motifsRBSJs <-
            getMotifs(
                targetsRBSJs,
                width = 6,
                species = "Mmusculus",
                rbp = TRUE,
                reverse = TRUE
            )


        # Merge motifs
        mergedMotifs <- mergeMotifs(motifs)
        mergedMotifsRBSJs <- mergeMotifs(motifsRBSJs)

        # plot
        p <-
            plotMotifs(
                mergedMotifs,
                mergedMotifsRBSJs,
                log2FC = 2,
                nf1 = nrow(annotatedBSJs),
                nf2 = nrow(annotatedRBSJs)
            )
        expect_is(p, "list")

    })
