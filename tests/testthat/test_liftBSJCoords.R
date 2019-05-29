setwd(paste(getwd(), "testdata", sep = "/"))

context("Test that liftBSJCoords function works correctly")
test_that("liftBSJCoords() generates the right data structure", {
    gtf <- formatGTF(pathToGTF = "gencodeVM16.gtf")
    # Retrieve back-spliced junctions coordinates
    backSplicedJunctions <- getBackSplicedJunctions(gtf)
    mergedBSJunctions <- mergeBSJunctions(backSplicedJunctions, gtf)

    ## convert mm10 to hg38 and retrieve gwas snps
    # liftedBSJs <-
    #     liftBSJCoords(mergedBSJunctions,
    #         map = "mm10ToHg38",
    #         annotationHubID = "AH14528")
    #
    # basicColumns <- getBasicColNames()
    #
    # expect_is(liftedBSJs, "data.frame")
    #
    # expect_identical(colnames(liftedBSJs), c(basicColumns, "id_mm10"))
    #
    # expect_equal(ncol(liftedBSJs), length(basicColumns) + 1)
    #
    # # we know that our test input data has in total 12 unique identifier
    # exptectedrows <- 10
    # expect_equal(nrow(liftedBSJs), exptectedrows)

})

test_that("liftBSJCoords() generates a data frame with the right content", {
    gtf <- formatGTF(pathToGTF = "gencodeVM16.gtf")
    # Retrieve back-spliced junctions coordinates
    backSplicedJunctions <- getBackSplicedJunctions(gtf)
    mergedBSJunctions <- mergeBSJunctions(backSplicedJunctions, gtf)

    ## convert mm10 to hg38 and retrieve gwas snps
    # liftedBSJs <-
    #     liftBSJCoords(mergedBSJunctions,
    #         map = "mm10ToHg38",
    #         annotationHubID = "AH14528")
    #
    # expect_equal(liftedBSJs$chrom[10], "chr14")
    # expect_equal(liftedBSJs$startUpBSE[10], 32090502)
    # expect_equal(liftedBSJs$endDownBSE[10], 32117287)

})
