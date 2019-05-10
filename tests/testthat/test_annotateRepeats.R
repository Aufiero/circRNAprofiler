setwd(paste(getwd(), "testdata" , sep = "/"))

context("Test that annotateRepeats() function works correctly")
test_that("annotateRepeats() generates the correct data structure", {

    gtf <- formatGTF(pathToGTF = "gencodeVM16.gtf")
    # Retrieve back-spliced junctions coordinates
    backSplicedJunctions <- getBackSplicedJunctions(gtf)
    mergedBSJunctions <- mergeBSJunctions(backSplicedJunctions, gtf)

    # Retrieve the genomic features of the circRNAs
    annotatedBSJs <- annotateBSJs(mergedBSJunctions, gtf)

    # Retrieve sequences type = "ie"
    targets <-
        getSeqsFromGRs(
            annotatedBSJs,
            lIntron = 500,
            lExon = 10,
            type = "ie",
            species = "Mmusculus",
            genome = "mm10"
        )

    # Retrieve overlapping repeats
    repeats <-
        annotateRepeats(
            targets,
            annotationHubID  = "AH6075",
            complementary = TRUE
        )

    expect_is(repeats, "list")
    expect_equal(length(repeats), 2)


    # Retrieve sequences type = "fi"
    targets <-
        getSeqsFromGRs(
            annotatedBSJs,
            lIntron = 500,
            lExon = 10,
            type = "fi",
            species = "Mmusculus",
            genome = "mm10"
        )

    # Retrieve overlapping repeats
    repeats <-
        annotateRepeats(
            targets,
            annotationHubID  = "AH6075",
            complementary = TRUE
        )

    expect_is(repeats, "list")
    expect_equal(length(repeats), 2)


})
