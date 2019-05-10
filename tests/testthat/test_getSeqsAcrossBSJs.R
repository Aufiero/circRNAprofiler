setwd(paste(getwd(), "testdata", sep = "/"))

context("Test that getSeqsAcrossBSJs() function works correctly")
test_that("getSeqsAcrossBSJs() retrieves the correct sequences", {
    gtf <- formatGTF(pathToGTF = "gencodeVM16.gtf")

    # Create the backSplicedJunctions data frame
    backSplicedJunctions <- getBackSplicedJunctions(gtf)
    mergedBSJunctions <- mergeBSJunctions(backSplicedJunctions, gtf)
    # Retrive the genomic features
    annotatedBSJs <- annotateBSJs(mergedBSJunctions, gtf)

    # retrieve target sequences
    targets <- getSeqsAcrossBSJs(annotatedBSJs,
        gtf,
        species = "Mmusculus",
        genome = "mm10")

    # For positive strand
    expect_equal(targets$bsj$id[11], "Arhgap5:+:chr12:52516079:52542636")
    expect_equal(targets$bsj$length[11], 22)

    # The back-spliced sequences should be the one reported below
    expect_equal(targets$bsj$seq[11],  "UGAAGACACAGAGGAAGAUGAU")

    #For negative strand
    expect_equal(targets$bsj$id[7], "Eps15l1:-:chr8:72380306:72367904")
    expect_equal(targets$bsj$length[7], 22)

    # The back-spliced sequences should be the one reported below
    expect_equal(targets$bsj$seq[7], "AGAUGUCCAAGAUCUCAUCAUU")


})
