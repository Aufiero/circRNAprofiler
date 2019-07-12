setwd(file.path(getwd(), "testdata"))

context("Test that getRandomBSJunctions() function works correctly")
test_that("getRandomBSJunctions() generate the correct data
          structure",
          {
              gtf <- formatGTF(pathToGTF = "gencodeVM16.gtf")
              # Check rows
              n <- 6 #arbitrary
              f <- 10

              randomBSJunctions <- getRandomBSJunctions (gtf, n = n, f = f)
              expect_is(randomBSJunctions, "data.frame")
              expect_equal(nrow(randomBSJunctions), n)

              # Check columns
              basicColumns <- .getBasicColNames()
              expect_equal(colnames(randomBSJunctions), basicColumns)

          })


test_that("getRandomBSJunctions() generate a data frame with the correct
    content", {

    gtf <- formatGTF(pathToGTF = "gencodeVM16.gtf")
    # Check rows
    n <- 6 #arbitrary
    f <- 10

    randomBSJunctions <- getRandomBSJunctions (gtf, n = n, f = f)

    # check the coordinates based on the strand
    if (length(unique(randomBSJunctions$strand)) == 2) {
        # For the positive strand the coordinate of the start_BSExon must be less
        # than the coordinate of the end_BSExon
        exonJunctions_pos <-
            randomBSJunctions[randomBSJunctions$strand == "+",][1, ]
        expect_lt(exonJunctions_pos$startUpBSE,
                  exonJunctions_pos$endDownBSE)

        # For the negative strand the coordinate of the start_BSExon must be
        # greater than the coordinate of the end_BSExon
        exonJunctions_neg <-
            randomBSJunctions[randomBSJunctions$strand == "-",][1, ]
        expect_gt(exonJunctions_neg$startUpBSE,
                  exonJunctions_neg$endDownBSE)

    } else if (unique(randomBSJunctions$strand) == "+") {
        exonJunctions_pos <-
            randomBSJunctions[randomBSJunctions$strand == "+",][1, ]
        expect_lt(exonJunctions_pos$startUpBSE,
                  exonJunctions_pos$endDownBSE)

    } else if (unique(randomBSJunctions$strand) == "-") {
        exonJunctions_neg <-
            randomBSJunctions[randomBSJunctions$strand == "-",][1, ]
        exonJunctions_neg <-
            randomBSJunctions[randomBSJunctions$strand == "-",][1, ]
        expect_gt(exonJunctions_neg$startUpBSE,
                  exonJunctions_neg$endDownBSE)
    }

})

