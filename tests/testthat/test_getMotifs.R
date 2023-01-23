setwd(file.path(getwd(), "testdata"))

context("Test that getMotifs() function works correctly")
test_that("getMotifs() and mergeMotifs() generate the correct data structure
    with GR seqs",
          {
              gtf <- formatGTF(pathToGTF = "gencodeVM16.gtf")

              # Create the backSplicedJunctions without retriving the missing coordinates
              backSplicedJunctions <- getBackSplicedJunctions(gtf)
              mergedBSJunctions <-
                  mergeBSJunctions(backSplicedJunctions, gtf)

              # Annonate BSJs
              annotatedBSJs <- annotateBSJs(mergedBSJunctions[9,], gtf)

              if (requireNamespace("BSgenome.Mmusculus.UCSC.mm10", quietly = TRUE)) {
                  # Get BSgenome object
                  genome <-
                      BSgenome::getBSgenome("BSgenome.Mmusculus.UCSC.mm10")

                  # Retrieve sequences
                  targets <-
                      getSeqsFromGRs(
                          annotatedBSJs,
                          genome,
                          lIntron = 101,
                          lExon = 9,
                          type = "ie"
                      )

                  # Retrieve motifs
                  # rbp = TRUE reverse = TRUE
                  motifs <-
                      getMotifs(
                          targets,
                          width = 6,
                          database = 'ATtRACT',
                          species = "Mmusculus",
                          rbp = TRUE,
                          reverse = TRUE
                      )
                  expect_is(motifs, "list")
                  expect_equal(length(motifs), 2)

                  # check nrow all data frame
                  expect_equal(nrow(motifs$upGR$targets), nrow(annotatedBSJs))
                  expect_equal(nrow(motifs$downGR$targets), nrow(annotatedBSJs))
                  expect_equal(ncol(motifs$upGR$counts) - 1, nrow(motifs$upGR$motif))

                  # Merge motifs
                  mergedMotifs <- mergeMotifs(motifs)
                  expect_is(mergedMotifs, "data.frame")

                  # Retrieve motifs
                  # rbp = FALSE reverse = FALSE
                  motifs <-
                      getMotifs(
                          targets,
                          width = 6,
                          database = 'ATtRACT',
                          species = "Mmusculus",
                          rbp = FALSE,
                          reverse = FALSE
                      )
                  expect_is(motifs, "list")
                  expect_equal(length(motifs), 2)

                  # check nrow all data frame
                  expect_equal(nrow(motifs$upGR$targets), nrow(annotatedBSJs))
                  expect_equal(nrow(motifs$downGR$targets), nrow(annotatedBSJs))
                  expect_equal(ncol(motifs$upGR$counts) - 1, nrow(motifs$upGR$motif))

                  # Merge motifs
                  mergedMotifs <- mergeMotifs(motifs)
                  expect_is(mergedMotifs, "data.frame")
              } else{
                  cat(
                      "Missing package BSgenome.Mmusculus.UCSC.mm10.
            Use BiocManager to install it."
                  )
              }



          })


test_that("getMotifs() generates a list with the correct
    content with GR seqs",
          {
              gtf <- formatGTF(pathToGTF = "gencodeVM16.gtf")

              # Create the backSplicedJunctions without retriving the missing coordinates
              backSplicedJunctions <- getBackSplicedJunctions(gtf)
              mergedBSJunctions <-
                  mergeBSJunctions(backSplicedJunctions, gtf)

              # Annonate BSJs
              annotatedBSJs <- annotateBSJs(mergedBSJunctions[4,], gtf)

              if (requireNamespace("BSgenome.Mmusculus.UCSC.mm10", quietly = TRUE)) {
                  # Get BSgenome object
                  genome <-
                      BSgenome::getBSgenome("BSgenome.Mmusculus.UCSC.mm10")

                  # Retrieve sequences
                  targets <-
                      getSeqsFromGRs(
                          annotatedBSJs,
                          genome,
                          lIntron = 101,
                          lExon = 9,
                          type = "ie"
                      )

                  # Retrieve motifs
                  motifs <-
                      getMotifs(
                          targets,
                          width = 4,
                          database = 'ATtRACT',
                          species = "Mmusculus",
                          rbp = TRUE,
                          reverse = TRUE
                      )

                  expect_equal(motifs$upGR$counts$UCUU[1], 0)
                  expect_equal(motifs$downGR$counts$UCUU[1], 1)

                  expect_equal(motifs$upGR$counts$UUCU[1], 0)
                  expect_equal(motifs$downGR$counts$UUCU[1], 2)


                  expect_equal(motifs$upGR$location$UUCU[1], "NA")
                  expect_equal(motifs$downGR$location$UUCU[1], "26,102")

                  expect_equal(motifs$upGR$location$UCUU[1], "NA")
                  expect_equal(motifs$downGR$location$UCUU[1], "103")

                  # check 2nd pattern
                  expect_equal(motifs$upGR$counts$AGAG[1], 2)
                  expect_equal(motifs$downGR$counts$AGAG[1], 2)

                  expect_equal(motifs$upGR$counts$GAGA[1], 1)
                  expect_equal(motifs$downGR$counts$GAGA[1], 2)

                  expect_equal(motifs$upGR$location$GAGA[1], "13")
                  expect_equal(motifs$downGR$location$GAGA[1], "18,56")

                  expect_equal(motifs$upGR$location$AGAG[1], "12,14")
                  expect_equal(motifs$downGR$location$AGAG[1], "17,55")
              } else{
                  cat(
                      "Missing package BSgenome.Mmusculus.UCSC.mm10.
            Use BiocManager to install it."
                  )
              }

          })




test_that(
    "getMotifs() and mergeMotifs() generates the correct data structure
    with circRNA and BSJ seqs",
    {
        gtf <- formatGTF(pathToGTF = "gencodeVM16.gtf")

        # Create the backSplicedJunctions without retriving the missing coordinates
        backSplicedJunctions <- getBackSplicedJunctions(gtf)
        mergedBSJunctions <-
            mergeBSJunctions(backSplicedJunctions, gtf)

        # Annonate BSJs
        annotatedBSJs <-
            annotateBSJs(mergedBSJunctions[9,], gtf)

        if (requireNamespace("BSgenome.Mmusculus.UCSC.mm10", quietly = TRUE)) {
            # Get BSgenome object
            genome <-
                BSgenome::getBSgenome("BSgenome.Mmusculus.UCSC.mm10")

            # Retrieve sequences (BSJs)
            targets <-
                getSeqsAcrossBSJs(annotatedBSJs,
                                  gtf,
                                  genome)

            # Retrieve motifs
            # rbp = TRUE reverse = TRUE
            motifs <-
                getMotifs(
                    targets,
                    width = 4,
                    species = "Mmusculus",
                    rbp = TRUE,
                    reverse = TRUE
                )
            expect_is(motifs, "list")
            expect_equal(length(motifs), 1)

            # check nrow all data frame
            expect_equal(nrow(motifs$bsj$targets), nrow(annotatedBSJs))
            expect_equal(ncol(motifs$bsj$counts) - 1, nrow(motifs$bsj$motif))

            # Merge motifs
            mergedMotifs <- mergeMotifs(motifs)
            expect_is(mergedMotifs, "data.frame")


            # Retrieve sequences (circRNA seqs)
            targets <-
                getCircSeqs(annotatedBSJs,
                            gtf,
                            genome)

            # Retrieve motifs
            # rbp = TRUE reverse = TRUE
            motifs <-
                getMotifs(
                    targets,
                    width = 4,
                    species = "Mmusculus",
                    rbp = TRUE,
                    reverse = TRUE
                )

            expect_is(motifs, "list")
            expect_equal(length(motifs), 1)

            # check nrow all data frame
            expect_equal(nrow(motifs$circ$targets), nrow(annotatedBSJs))
            expect_equal(ncol(motifs$circ$counts) - 1, nrow(motifs$circ$motif))

            # Merge motifs
            mergedMotifs <- mergeMotifs(motifs)
            expect_is(mergedMotifs, "data.frame")
        } else{
            cat(
                "Missing package BSgenome.Mmusculus.UCSC.mm10.
            Use BiocManager to install it."
            )
        }


})
