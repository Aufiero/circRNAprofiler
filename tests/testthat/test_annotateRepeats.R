setwd(file.path(getwd(), "testdata"))

context("Test that annotateRepeats() function works correctly")
test_that("annotateRepeats() generates the correct data structure", {

    gtf <- formatGTF(pathToGTF = "gencodeVM16.gtf")
    # Retrieve back-spliced junctions coordinates
    backSplicedJunctions <- getBackSplicedJunctions(gtf)
    mergedBSJunctions <- mergeBSJunctions(backSplicedJunctions, gtf)

    # Retrieve the genomic features of the circRNAs
    annotatedBSJs <- annotateBSJs(mergedBSJunctions[1, ], gtf)

    if (requireNamespace("BSgenome.Mmusculus.UCSC.mm10", quietly = TRUE)){
        # Get BSgenome object
        genome <- BSgenome::getBSgenome("BSgenome.Mmusculus.UCSC.mm10")

        # Retrieve sequences type = "ie"
        targets <-
            getSeqsFromGRs(
                annotatedBSJs,
                genome,
                lIntron = 500,
                lExon = 10,
                type = "ie"

            )

        # Retrieve overlapping repeats

        # repeats <-
        #     annotateRepeats(
        #         targets,
        #         annotationHubID  = "AH6075",
        #         complementary = TRUE
        #     )
        # 
        # expect_is(repeats, "list")
        # expect_equal(length(repeats), 2)
    }else{
        cat(
            "Missing package BSgenome.Mmusculus.UCSC.mm10.
            Use BiocManager to install it."
        )
    }





})
