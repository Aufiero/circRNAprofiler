setwd(paste(getwd(), "testdata", sep = "/"))

context("Test that getLength() function works correctly")
test_that("getLength() generates the correct data structure", {

    experiment <- read.table(
        "experiment.txt",
        header = TRUE,
        stringsAsFactors = FALSE,
        sep = "\t"
    )

    gtf <- formatGTF(pathToGTF = "gencodeVM16.gtf")

    # Create the backSplicedJunctions data frame
    backSplicedJunctions <- getBackSplicedJunctions(gtf)
    mergedBSJunctions <- mergeBSJunctions(backSplicedJunctions, gtf)

    # Retrive the genomic features
    annotatedBSJs <- annotateBSJs(mergedBSJunctions, gtf)

    exInLength <- getLength(annotatedBSJs)
    expect_is(exInLength, "data.frame")

    colNames <- c(
        "id",
        "lenUpIntron",
        "lenUpBSE",
        "lenDownBSE",
        "lenDownIntron",
        "meanLengthBSEs",
        "meanLengthIntrons"
    )
    expect_equal(colnames(exInLength), colNames)
    expect_equal(nrow(exInLength), nrow(annotatedBSJs))


})
