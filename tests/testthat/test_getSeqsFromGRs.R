setwd(file.path(getwd(), "testdata"))

context("Test that getSeqsFromGRs() function works correctly")
test_that("getSeqsFromGRs() generates the correct data structure", {
    gtf <- formatGTF(pathToGTF = "gencodeVM16.gtf")
    # Create the backSplicedJunctions without retriving the missing coordinates
    backSplicedJunctions <- getBackSplicedJunctions(gtf)
    mergedBSJunctions <-
        mergeBSJunctions(backSplicedJunctions, gtf)
    
    # Annonate BSJs
    annotatedBSJs <- annotateBSJs(mergedBSJunctions, gtf)
    
    if (requireNamespace("BSgenome.Mmusculus.UCSC.mm10", quietly = TRUE)) {
        # Get BSgenome object
        genome <-
            BSgenome::getBSgenome("BSgenome.Mmusculus.UCSC.mm10")
        
        # type = "ie"
        targets <-
            getSeqsFromGRs(
                annotatedBSJs,
                genome,
                lIntron = 100,
                lExon = 9,
                type = "ie"
            )
        
        expect_is(targets, "list")
        expect_equal(
            colnames(targets$upGR),
            c(
                "id",
                "gene",
                "transcript",
                "strand",
                "chrom",
                "startGR",
                "endGR",
                "length",
                "seq",
                "type"
            )
        )
        expect_equal(nrow(targets$upGR), nrow(annotatedBSJs))
        
        
        # type = "bse"
        targets <-
            getSeqsFromGRs(
                annotatedBSJs,
                genome,
                lIntron = 100,
                lExon = 9,
                type = "bse"
            )
        
        expect_is(targets, "list")
        expect_equal(
            colnames(targets$upGR),
            c(
                "id",
                "gene",
                "transcript",
                "strand",
                "chrom",
                "startGR",
                "endGR",
                "length",
                "seq",
                "type"
            )
        )
        expect_equal(nrow(targets$upGR), nrow(annotatedBSJs))
        
        # type = "fi"
        targets <-
            getSeqsFromGRs(
                annotatedBSJs,
                genome,
                lIntron = 100,
                lExon = 9,
                type = "fi"
            )
        
        expect_is(targets, "list")
        expect_equal(
            colnames(targets$upGR),
            c(
                "id",
                "gene",
                "transcript",
                "strand",
                "chrom",
                "startGR",
                "endGR",
                "length",
                "seq",
                "type"
            )
        )
        expect_equal(nrow(targets$upGR), nrow(annotatedBSJs))
        
    } else{
        cat(
            "Missing package BSgenome.Mmusculus.UCSC.mm10.
            Use BiocManager to install it."
        )
    }
    
})

test_that("getSeqFromGRs() retrieves the correct genomic ranges when type = ie",
          {
              gtf <- formatGTF(pathToGTF = "gencodeVM16.gtf")
              # Create the backSplicedJunctions without retriving the missing coordinates
              backSplicedJunctions <- getBackSplicedJunctions(gtf)
              mergedBSJunctions <-
                  mergeBSJunctions(backSplicedJunctions, gtf)
              
              # Annonate BSJs
              annotatedBSJs <- annotateBSJs(mergedBSJunctions, gtf)
              
              if (requireNamespace("BSgenome.Mmusculus.UCSC.mm10", quietly = TRUE)) {
                  # Get BSgenome object
                  genome <-
                      BSgenome::getBSgenome("BSgenome.Mmusculus.UCSC.mm10")
                  
                  targets <- getSeqsFromGRs(
                      annotatedBSJs,
                      genome,
                      lIntron = 100,
                      lExon = 9,
                      type = "ie"
                  )
                  
                  # Check Raph1:-:chr1:60533406:60525592
                  expect_lt(targets$upGR$startGR[4], targets$upGR$endGR[4])
                  expect_lt(targets$downGR$startGR[4], targets$downGR$endGR[4])
                  # Check Creb1:+:chr1:64550849:64576329
                  expect_lt(targets$upGR$startGR[3], targets$upGR$endGR[3])
                  expect_lt(targets$downGR$startGR[3], targets$downGR$endGR[3])
                  
                  
                  # Check
                  expect_equal(targets$upGR$startGR[4], 60533397)
                  expect_equal(targets$upGR$endGR[4], 60533506)
                  expect_equal(targets$downGR$startGR[4], 60525492)
                  expect_equal(targets$downGR$endGR[4], 60525601)
                  
                  # Check
                  expect_equal(targets$upGR$startGR[3], 64550749)
                  expect_equal(targets$upGR$endGR[3], 64550858)
                  expect_equal(targets$downGR$startGR[3], 64576320)
                  expect_equal(targets$downGR$endGR[3], 64576429)
                  
                  # Check
                  expect_equal(targets$upGR$startGR[9], 155440776)
                  expect_equal(targets$upGR$endGR[9], 155440887)
                  expect_equal(targets$downGR$startGR[9], 155437759)
                  expect_equal(targets$downGR$endGR[9], 155437869)
                  
                  # Check
                  expect_equal(targets$upGR$startGR[11], 52516077)
                  expect_equal(targets$upGR$endGR[11], 52516086)
                  expect_equal(targets$downGR$startGR[11], 52542627)
                  expect_equal(targets$downGR$endGR[11], 52542736)
                  
                  # Check
                  expect_equal(targets$upGR$startGR[1], 43704434)
                  expect_equal(targets$upGR$endGR[1], 43704543)
                  expect_equal(targets$downGR$startGR[1], 43705468)
                  expect_equal(targets$downGR$endGR[1], 43705577)
                  
                  # Check whether the right sequences are retrieved
                  expect_equal(substr(targets$upGR$seq[4], 1, 20),
                               "AUUUAACUUCAGAGAGUCAC")
                  expect_equal(substr(targets$downGR$seq[4], 1, 20),
                               "AAUUACUGAGGUAAAUAGAG")
                  
                  expect_equal(substr(targets$upGR$seq[3], 1, 20),
                               "UUAACCACCACAUAGCCUAU")
                  expect_equal(substr(targets$downGR$seq[3], 1, 20),
                               "UGAAGAACAGGUACAGAUUC")
              } else{
                  cat(
                      "Missing package BSgenome.Mmusculus.UCSC.mm10.
            Use BiocManager to install it."
                  )
              }
              
          })



test_that("getSeqsFromGRs() extracts BSJ coordinates when type = bse", {
    gtf <- formatGTF(pathToGTF = "gencodeVM16.gtf")
    # Create the backSplicedJunctions without retriving the missing coordinates
    backSplicedJunctions <- getBackSplicedJunctions(GTFFile())
    mergedBSJunctions <- mergeBSJunctions(backSplicedJunctions, gtf)
    
    # Annonate BSJs
    annotatedBSJs <- annotateBSJs(mergedBSJunctions, gtf)
    
    if (requireNamespace("BSgenome.Mmusculus.UCSC.mm10", quietly = TRUE)) {
        # Get BSgenome object
        genome <-
            BSgenome::getBSgenome("BSgenome.Mmusculus.UCSC.mm10")
        
        targets <- getSeqsFromGRs(annotatedBSJs,
                                  genome,
                                  type = "bse")
        
        expect_lt(targets$upGR$startGR[4], targets$upGR$endGR[4])
        expect_lt(targets$downGR$startGR[4], targets$downGR$endGR[4])
        expect_lt(targets$upGR$startGR[3], targets$upGR$endGR[3])
        expect_lt(targets$downGR$startGR[3], targets$downGR$endGR[3])
        
        
        #we know that for the reported row of our test data the strand is negative
        expect_equal(targets$upGR$startGR[4], 60533287)
        expect_equal(targets$upGR$endGR[4], 60533406)
        expect_equal(targets$downGR$startGR[4], 60525592)
        expect_equal(targets$downGR$endGR[4], 60526100)
        
        # we know that for the reported row in our test data the strand is positive
        expect_equal(targets$upGR$startGR[3], 64550849)
        expect_equal(targets$upGR$endGR[3], 64550970)
        expect_equal(targets$downGR$startGR[3], 64576179)
        expect_equal(targets$downGR$endGR[3], 64576329)
        
        # Check whether the right sequences are retrieved
        expect_equal(substr(targets$upGR$seq[4], 1, 20),
                     "AUGGAGCAACUAUCCGAUGA")
        expect_equal(substr(targets$downGR$seq[4], 1, 20),
                     "AAGCUCUGAAUCAGGGUGAG")
        
        expect_equal(substr(targets$upGR$seq[3], 1, 20),
                     "GUAACUAAAUGACCAUGGAA")
        expect_equal(substr(targets$downGR$seq[3], 1, 20),
                     "CUGCCUCAGGCGAUGUACAA")
    } else{
        cat(
            "Missing package BSgenome.Mmusculus.UCSC.mm10.
            Use BiocManager to install it."
        )
    }
    
    
})


test_that("getSeqsFromGRs() extracts the coordinates of the flanking introns
    whhen type = fi",
          {
              gtf <- formatGTF(pathToGTF = "gencodeVM16.gtf")
              # Create the backSplicedJunctions without retriving the missing coordinates
              backSplicedJunctions <- getBackSplicedJunctions(gtf)
              mergedBSJunctions <-
                  mergeBSJunctions(backSplicedJunctions, gtf)
              
              # Annonate BSJs
              annotatedBSJs <- annotateBSJs(mergedBSJunctions, gtf)
              
              if (requireNamespace("BSgenome.Mmusculus.UCSC.mm10", quietly = TRUE)) {
                  # Get BSgenome object
                  genome <-
                      BSgenome::getBSgenome("BSgenome.Mmusculus.UCSC.mm10")
                  
                  targets <- getSeqsFromGRs(annotatedBSJs,
                                            genome,
                                            type = "fi")
                  
                  
                  expect_lt(targets$upGR$startGR[4], targets$upGR$endGR[4])
                  expect_lt(targets$downGR$startGR[4], targets$downGR$endGR[4])
                  expect_lt(targets$upGR$startGR[3], targets$upGR$endGR[3])
                  expect_lt(targets$downGR$startGR[3], targets$downGR$endGR[3])
                  
                  # we know that for the reported row of our test data the strand is negative
                  expect_equal(targets$upGR$startGR[4], 60533407)
                  expect_equal(targets$upGR$endGR[4], 60566564)
                  expect_equal(targets$downGR$startGR[4], 60503596)
                  expect_equal(targets$downGR$endGR[4], 60525591)
                  
                  # we know that for the reported row in our test data the strand is positive
                  expect_equal(targets$upGR$startGR[3], 64533021)
                  expect_equal(targets$upGR$endGR[3], 64550848)
                  expect_equal(targets$downGR$startGR[3], 64576330)
                  expect_equal(targets$downGR$endGR[3], 64597233)
                  
                  # Check whether the right sequences are retrieved
                  expect_equal(substr(targets$upGR$seq[4], 1, 20),
                               "GUGAGUGCUCCGCCACCUGG")
                  expect_equal(substr(targets$downGR$seq[4], 1, 20),
                               "GUAAAUAGAGAAAAUUUCUC")
                  
                  expect_equal(substr(targets$upGR$seq[3], 1, 20),
                               "GUAAGAGGAGCAGGAGGAGG")
                  expect_equal(substr(targets$downGR$seq[3], 1, 20),
                               "GUACAGAUUCCAAACACUUA")
              } else{
                  cat(
                      "Missing package BSgenome.Mmusculus.UCSC.mm10.
            Use BiocManager to install it."
                  )
              }
              
              
          })
