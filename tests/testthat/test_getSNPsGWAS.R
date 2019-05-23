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

        if (!requireNamespace("BSgenome.Mmusculus.UCSC.mm10", quietly = TRUE))
            install.packages("BSgenome.Mmusculus.UCSC.mm10")

        # Get BSgenome object
        genome <- BSgenome::getBSgenome("BSgenome.Hsapiens.UCSC.hg38")

        # Retrieve sequences
        targets <-
            getSeqsFromGRs(
                annotatedBSJs,
                genome,
                lIntron = 500,
                lExon = 10,
                type = "fi"
            )

        #options(warn = -1)
        # Retrieve overlapping gwas
        SNPsGWAS <- annotateSNPsGWAS(targets, assembly = "hg38", makeCurrent = FALSE)
        #options(warn = 0) # default

        expect_is(SNPsGWAS, "list")
        expect_equal(length(SNPsGWAS), 2)
        })
