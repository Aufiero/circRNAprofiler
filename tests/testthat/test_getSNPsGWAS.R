setwd(paste(getwd(), "testdata", sep = "/"))

context("Test that annotateSNPsGWAS() function works correctly")
test_that("annotateSNPsGWAS() generates the correct
    data structure", {
        gtf <- formatGTF(pathToGTF = "gencodeVM16.gtf")

        # Retrieve back-spliced junctions coordinates
        backSplicedJunctions <- getBackSplicedJunctions(gtf)
        mergedBSJunctions <- mergeBSJunctions(backSplicedJunctions, gtf)
        ### convert mm10 to hg38 and retrieve gwas snps
        liftedBSJs <-
            liftBSJCoords(mergedBSJunctions,
                map = "mm10ToHg38",
                annotationHubID = "AH14528")

        # Retrieve the genomic features of the circRNAs
        annotatedBSJs <- annotateBSJs(liftedBSJs, gtf)

        # Retrieve sequences
        targets <-
            getSeqsFromGRs(
                annotatedBSJs,
                lIntron = 500,
                lExon = 10,
                type = "fi",
                species = "Hsapiens",
                genome = "hg38"
            )

        #options(warn = -1)
        # Retrieve overlapping gwas
        SNPsGWAS <- annotateSNPsGWAS(targets, genome = "hg38")
        #options(warn = 0) # default

        expect_is(SNPsGWAS, "list")
        expect_equal(length(SNPsGWAS), 2)
        })
