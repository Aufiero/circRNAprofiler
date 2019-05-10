setwd(paste(getwd(), "testdata", sep = "/"))

context("Test that getCircSeqs() function works correctly")
test_that("getCircSeqs() generate the correct data structure", {

    gtf <- formatGTF(pathToGTF = "gencodeVM16.gtf")

    # Create the backSplicedJunctions data frame
    backSplicedJunctions <- getBackSplicedJunctions(gtf)
    mergedBSJunctions <- mergeBSJunctions(backSplicedJunctions, gtf)
    # Annotate BSJs
    annotatedBSJs <- annotateBSJs(mergedBSJunctions, gtf)

    # retrieve target sequences
    targets <- getCircSeqs(annotatedBSJs,
        gtf,
        species = "Mmusculus",
        genome = "mm10"
        )

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


})

test_that("getCircSeqs() retrieves the right sequences", {

    gtf <- formatGTF(pathToGTF = "gencodeVM16.gtf")

    # Create the backSplicedJunctions data frame
    backSplicedJunctions <- getBackSplicedJunctions(gtf)
    mergedBSJunctions <-
        mergeBSJunctions(backSplicedJunctions, gtf)
    # Annotate BSJs
    annotatedBSJs <- annotateBSJs(mergedBSJunctions, gtf)

    # retrieve target sequences
    targets <- getCircSeqs(annotatedBSJs,
        gtf,
        species = "Mmusculus",
        genome = "mm10"
    )

    expect_equal(nrow(targets$circ), nrow(annotatedBSJs))

    # For positive strand
    expect_equal(targets$circ$id[11], "Arhgap5:+:chr12:52516079:52542636")
    expect_equal(targets$circ$length[11], 4039)
    # To recreate the back-spliced seq the first 50 nucleotides of the upstream
    # back-spliced exon are added to the the end of the downstream
    # back-spliced exons. So if lenght of the first circRNA should be 50
    # nucleotides longer than the actual lenght
    expect_equal(nchar(targets$circ$seq[11]), 4089)

    # The fist 40 nucleotides of the upstream back-sliced exons ahould be equal
    # to the last 30 nuclotides of the downstream back-spliced exon.
    expect_equal(substr(targets$circ$seq[11], 1, 50),
        substr(targets$circ$seq[11], 4040, 4089))

    # The back-spliced sequences should be the one reported below
    expect_equal(
        substr(targets$circ$seq[11], 4020, 4089),
        "GGAAUUUAUUGAAGACACAGAGGAAGAUGAUCCAUAUGAUCUCGAAGAAGACCUUUCUUUGGCAAUGAGG"
    )


    #For negative strand
    expect_equal(targets$circ$id[7], "Eps15l1:-:chr8:72380306:72367904")
    expect_equal(targets$circ$length[7], 819)
    # To recreate the back-spliced seq the first 40 nucleotides of the upstream
    # back-spliced exon are added to the the end of the downstream
    # back-spliced exons. The lenght of the first circRNA should be 50
    # nucleotides longer than the actual lenght.
    expect_equal(nchar(targets$circ$seq[7]), 869)

    # The fist 30 nucleotides of the upstream back-sliced exons ahould be equal
    # to the last 50 nuclotides of the
    # downstream back-spliced exon.
    expect_equal(substr(targets$circ$seq[7], 1, 50),
        substr(targets$circ$seq[7], 820, 869))

    # The back-spliced sequences should be the one reported below
    expect_equal(
        substr(targets$circ$seq[7], 800, 869),
        "ACUUCAGUCAGAUGUCCAAGAUCUCAUCAUUGAAAACCCAGAUUCAGUCUCAGGAGUCAGACUUGAAGUC"
    )


})
