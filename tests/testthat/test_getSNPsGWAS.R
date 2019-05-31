setwd(paste(getwd(), "testdata", sep = "/"))

context("Test that liftBSJCoords() and annotateSNPsGWAS() function works correctly")
test_that("annotateSNPsGWAS() and liftBSJCoords() generates the correct
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


        basicColumns <- getBasicColNames()

        expect_is(liftedBSJs, "data.frame")

        expect_identical(colnames(liftedBSJs), c(basicColumns, "id_mm10"))

        expect_equal(ncol(liftedBSJs), length(basicColumns) + 1)

        # we know that our test input data has in total 12 unique identifier
        exptectedrows <- 10
        expect_equal(nrow(liftedBSJs), exptectedrows)

        expect_equal(liftedBSJs$chrom[10], "chr14")
        expect_equal(liftedBSJs$startUpBSE[10], 32090502)
        expect_equal(liftedBSJs$endDownBSE[10], 32117287)



        # Retrieve the genomic features of the circRNAs
        gtf <- formatGTF(pathToGTF = "gencode.v27.gtf")
        annotatedBSJs <- annotateBSJs(liftedBSJs[1,], gtf)

        if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE))
            install.packages("BSgenome.Hsapiens.UCSC.hg38")

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

        gwas <- getGWAS(assembly= "hg38")

        check <- grep("chr",as.character(seqnames(gwas[1,])))
        expect_equal(check, 1)

        #options(warn = -1)
        # Retrieve overlapping gwas
        SNPsGWAS <- annotateSNPsGWAS(targets, assembly = "hg38")
        #options(warn = 0) # default

        expect_is(SNPsGWAS, "list")
        expect_equal(length(SNPsGWAS), 2)
        })
