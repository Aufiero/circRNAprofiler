setwd(file.path(getwd(), "testdata"))

context("Test that getCircSeqs() function works correctly")
test_that("getCircSeqs() generate the correct data structure", {
    gtf <- formatGTF(pathToGTF = "gencodeVM16.gtf")
    
    # Create the backSplicedJunctions data frame
    backSplicedJunctions <- getBackSplicedJunctions(gtf)
    mergedBSJunctions <- mergeBSJunctions(backSplicedJunctions, gtf)
    # Annotate BSJs
    annotatedBSJs <- annotateBSJs(mergedBSJunctions[1, ], gtf)
    
    if (requireNamespace("BSgenome.Mmusculus.UCSC.mm10", quietly = TRUE)) {
        # Get BSgenome object
        genome <-
            BSgenome::getBSgenome("BSgenome.Mmusculus.UCSC.mm10")
        
        # retrieve target sequences
        targets <- getCircSeqs(annotatedBSJs,
                               gtf,
                               genome)
        
        expect_is(targets, "list")
        expect_equal(
            colnames(targets$circ),
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
        
    } else{
        cat(
            "Missing package BSgenome.Mmusculus.UCSC.mm10.
            Use BiocManager to install it."
        )
    }
    
})

test_that("getCircSeqs() retrieves the right sequences", {
    gtf <- formatGTF(pathToGTF = "gencodeVM16.gtf")
    
    # Create the backSplicedJunctions data frame
    backSplicedJunctions <- getBackSplicedJunctions(gtf)
    mergedBSJunctions <-
        mergeBSJunctions(backSplicedJunctions, gtf)
    # Annotate BSJs
    annotatedBSJs <- annotateBSJs(mergedBSJunctions[c(7, 11), ], gtf)
    
    if (requireNamespace("BSgenome.Mmusculus.UCSC.mm10", quietly = TRUE)) {
        # Load BSgenome object
        genome <-
            BSgenome::getBSgenome("BSgenome.Mmusculus.UCSC.mm10")
        
        # retrieve target sequences
        targets <- getCircSeqs(annotatedBSJs,
                               gtf,
                               genome)
        
        expect_equal(nrow(targets$circ), nrow(annotatedBSJs))
        
        # For positive strand
        expect_equal(targets$circ$id[2],
                     "Arhgap5:+:chr12:52516079:52542636")
        expect_equal(targets$circ$length[2], 4039)
        # To recreate the back-spliced seq the first 50 nucleotides of the upstream
        # back-spliced exon are added to the the end of the downstream
        # back-spliced exons. So if lenght of the first circRNA should be 50
        # nucleotides longer than the actual lenght
        expect_equal(nchar(targets$circ$seq[2]), 4089)
        
        # The fist 40 nucleotides of the upstream back-sliced exons ahould be equal
        # to the last 30 nuclotides of the downstream back-spliced exon.
        expect_equal(substr(targets$circ$seq[2], 1, 50),
                     substr(targets$circ$seq[2], 4040, 4089))
        
        # The back-spliced sequences should be the one reported below
        expect_equal(
            substr(targets$circ$seq[2], 4020, 4089),
            "GGAAUUUAUUGAAGACACAGAGGAAGAUGAUCCAUAUGAUCUCGAAGAAGACCUUUCUUUGGCAAUGAGG"
        )
        
        
        #For negative strand
        expect_equal(targets$circ$id[1],
                     "Eps15l1:-:chr8:72380306:72367904")
        expect_equal(targets$circ$length[1], 819)
        # To recreate the back-spliced seq the first 40 nucleotides of the upstream
        # back-spliced exon are added to the the end of the downstream
        # back-spliced exons. The lenght of the first circRNA should be 50
        # nucleotides longer than the actual lenght.
        expect_equal(nchar(targets$circ$seq[1]), 869)
        
        # The fist 30 nucleotides of the upstream back-sliced exons ahould be equal
        # to the last 50 nuclotides of the
        # downstream back-spliced exon.
        expect_equal(substr(targets$circ$seq[1], 1, 50),
                     substr(targets$circ$seq[1], 820, 869))
        
        # The back-spliced sequences should be the one reported below
        expect_equal(
            substr(targets$circ$seq[1], 800, 869),
            "ACUUCAGUCAGAUGUCCAAGAUCUCAUCAUUGAAAACCCAGAUUCAGUCUCAGGAGUCAGACUUGAAGUC"
        )
    } else{
        cat(
            "Missing package BSgenome.Mmusculus.UCSC.mm10.
            Use BiocManager to install it."
        )
    }
    
    
})
