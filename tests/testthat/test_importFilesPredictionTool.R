setwd(file.path(getwd(), "testdata"))

context("test that the import functions work correctly")
test_that("importMapSplice() generate the correct data structure with
    correct content ", {

    pathToGTF <- "gencodeVM16.gtf"
    gtf <- formatGTF (pathToGTF)

    experiment <- read.table(
        "experiment.txt",
        header = TRUE,
        stringsAsFactors = FALSE,
        sep = "\t"
    )
    pathToFile <-
        paste("mapsplice", experiment$fileName[1], sep = "/")

    adaptedPatientBSJunctions <- importMapSplice(pathToFile)

    expect_is(adaptedPatientBSJunctions, "data.frame")

    #check 1 step
    column <- c(.getBasicColNames(), "coverage")
    expect_identical(colnames(adaptedPatientBSJunctions), column)

    # check 2 step
    # we know that the first row of the file test circular_RNAs_001.txt  has the values reported below
    expect_identical(adaptedPatientBSJunctions$chrom[4], "chr1")
    expect_identical(adaptedPatientBSJunctions$strand[4], "-")
    expect_identical(adaptedPatientBSJunctions$gene[4], "Raph1")
})


test_that("importOther() generate the correct data structure with
    correct content ", {

    experiment <-
        read.table(
            "experiment.txt",
            header = TRUE,
            stringsAsFactors = FALSE,
            sep = "\t"
        )
    pathToFile <- paste("tool1", experiment$fileName[1], sep = "/")

    adaptedPatientBSJunctions <- importOther(pathToFile)

    expect_is(adaptedPatientBSJunctions, "data.frame")

    # TODO: consider to include a test to check the content of the needed columns

    #check 1 step
    column <- c(.getBasicColNames(), "coverage")
    expect_identical(colnames(adaptedPatientBSJunctions), column)

    # check 2 step
    # we know that the first row of the file test circular_RNAs_001.txt  has the values reported below
    expect_identical(adaptedPatientBSJunctions$chrom[1], "chr1")
    expect_identical(adaptedPatientBSJunctions$strand[1], "-")
    expect_identical(adaptedPatientBSJunctions$gene[1], "Raph1")


})
