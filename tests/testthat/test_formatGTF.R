setwd(file.path(getwd(), "testdata"))
context("Test that formatGTF() function works correctly")
test_that("formatGTF() generates the correct data structure", {
    needColumns <-
        c(
            "chrom",
            "start",
            "end",
            "width",
            "strand",
            "type",
            "gene_name",
            "transcript_id",
            "exon_number"
        )

    pathToGTF <- "UCSCmm10.gtf"
    gtf <- formatGTF(pathToGTF)
    expect_equal(colnames(gtf), needColumns)
    expect_is(gtf, "data.frame")
    expect_equal(class(gtf$chrom), "character")
    expect_equal(class(gtf$strand), "character")
    expect_equal(class(gtf$exon_number), "numeric")

    pathToGTF <- "ncbi38.gtf"
    gtf <- formatGTF(pathToGTF)
    expect_equal(colnames(gtf), needColumns)
    expect_is(gtf, "data.frame")
    expect_equal(class(gtf$chrom), "character")
    expect_equal(class(gtf$strand), "character")
    expect_equal(class(gtf$exon_number), "numeric")

    pathToGTF <- "GRCm38.92.gtf"
    gtf <- formatGTF(pathToGTF)
    expect_equal(colnames(gtf), needColumns)
    expect_is(gtf, "data.frame")
    expect_equal(class(gtf$chrom), "character")
    expect_equal(class(gtf$strand), "character")
    expect_equal(class(gtf$exon_number), "numeric")

    pathToGTF <- "gencodeVM16.gtf"
    gtf <- formatGTF(pathToGTF)
    expect_equal(colnames(gtf), needColumns)
    expect_is(gtf, "data.frame")
    expect_equal(class(gtf$chrom), "character")
    expect_equal(class(gtf$strand), "character")
    expect_equal(class(gtf$exon_number), "numeric")

    
    # Check for chromosome MT

})
