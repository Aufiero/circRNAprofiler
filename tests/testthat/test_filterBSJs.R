setwd(file.path(getwd(), "testdata"))
context("Test that filterCirc() function works correctly")
test_that("filterCirc() filters the data frame containing circRNA counts", {

              gtf <- formatGTF(pathToGTF = "gencodeVM16.gtf")
              experiment <- read.table(
                  "experiment.txt",
                  header = TRUE,
                  stringsAsFactors = FALSE,
                  sep = "\t"
              )

              backSplicedJunctions <- getBackSplicedJunctions(gtf)
              mergedBSJunctions <- mergeBSJunctions(backSplicedJunctions, gtf)

              filteredCirc <-
                  filterCirc(mergedBSJunctions,  allSamples = TRUE, min = 0)

              expect_is(filteredCirc, "data.frame")
              expect_identical(colnames(filteredCirc), colnames(mergedBSJunctions))
              expect_equal(nrow(filteredCirc), nrow(mergedBSJunctions))

              # Test allSamples = TRUE, min=3
              filteredCirc <-
                  filterCirc(mergedBSJunctions, allSamples = TRUE, min = 8)

              expect_is(filteredCirc, "data.frame")
              expect_identical(colnames(filteredCirc), colnames(mergedBSJunctions))
              # we know that only 2 circRNAs pass the filtering step (have at least 3
              # counts in all samples)
              expect_equal(nrow(filteredCirc), 2)

              # Test allSamples = TRUE, min=3
              filteredCirc <-
                  filterCirc(mergedBSJunctions, allSamples = FALSE, min = 8)

              expect_is(filteredCirc, "data.frame")
              expect_identical(colnames(filteredCirc), colnames(mergedBSJunctions))
              # we know that only 4 circRNAs pass the filtering step (have at least 3
              # counts in all samples)
              expect_equal(nrow(filteredCirc), 4)

              # Test allSamples = TRUE, min=3
              filteredCirc <-
                  filterCirc(mergedBSJunctions, allSamples = TRUE, min = 1)

              expect_is(filteredCirc, "data.frame")
              expect_identical(colnames(filteredCirc), colnames(mergedBSJunctions))
              # we know that only 4 circRNAs pass the filtering step (have at least 3
              # counts in all samples)
              expect_equal(nrow(filteredCirc), 4)


          })
