setwd(file.path(getwd(), "testdata"))
context("Test that getBackSplicedJunctions() function works correctly")

test_that("getBackSplicedJunctions() generates the correct data structure", {

              experiment <- read.table(
                  "experiment.txt",
                  header = TRUE,
                  stringsAsFactors = FALSE,
                  sep = "\t"
              )
              gtf <- formatGTF(pathToGTF = "gencodeVM16.gtf")

              backSplicedJunctions <- getBackSplicedJunctions(gtf)

              basicColumns <- .getBasicColNames()

              expect_is(backSplicedJunctions, "data.frame")

              expect_identical(colnames(backSplicedJunctions),
                               c(basicColumns, "tool", experiment$label))

              expect_equal(ncol(backSplicedJunctions),
                           length(basicColumns) + nrow(experiment) + 1)


              # we know that our test input data has in total 6 unique identifier
              exptectedUniqieIdentifier <- 18
              expect_equal(nrow(backSplicedJunctions), exptectedUniqieIdentifier)

          })


test_that("mergeBSJunctions() generates a data frame with the correct content", {
    experiment <- read.table(
        "experiment.txt",
        header = TRUE,
        stringsAsFactors = FALSE,
        sep = "\t"
    )
    gtf <- formatGTF(pathToGTF = "gencodeVM16.gtf")

    backSplicedJunctions <- getBackSplicedJunctions(gtf)
    mergedBSJunctions <- mergeBSJunctions(backSplicedJunctions, gtf)

    # We know that the first unique identifier of our test input files is the
    # one below
    expect_identical(mergedBSJunctions$id[4], "Raph1:-:chr1:60533406:60525592")

    # For the positive strand the coordinate of the startUpBSE must be
    # less than the coordinate of the endDownBSE
    expect_lt(mergedBSJunctions$startUpBSE[2],
              mergedBSJunctions$endDownBSE[2])

    # For the negative strand the coordinate of the startUpBSE must be
    # greater than the coordinate of the endDownBSE
    expect_gt(mergedBSJunctions$startUpBSE[4],
              mergedBSJunctions$endDownBSE[1])

    expect_identical(mergedBSJunctions[mergedBSJunctions$id ==
            "Eps15l1:-:chr8:72380306:72367904", "tool"], "ms,t1")
    expect_identical(mergedBSJunctions[mergedBSJunctions$id ==
            "Arhgap5:+:chr12:52516079:52542636", "tool"], "t1")
    expect_identical(mergedBSJunctions[mergedBSJunctions$id ==
            "Ndufv1:-:chr19:4009435:4008653", "tool"], "t1")
    expect_identical(mergedBSJunctions[mergedBSJunctions$id ==
            "Raph1:-:chr1:60533406:60525592", "tool"], "ms,t1")

})
