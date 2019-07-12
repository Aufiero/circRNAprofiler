setwd(file.path(getwd(), "testdata"))
context("Test that annotateBSJs() function works correctly")

test_that(
    "annotateBSJs() generates the correct data structure",
    {
        gtf <- formatGTF(pathToGTF = "gencodeVM16.gtf")

        # Create the backSplicedJunctions data frame
        backSplicedJunctions <- getBackSplicedJunctions(gtf)
        mergedBSJunctions <-
            mergeBSJunctions(backSplicedJunctions, gtf)

        # Retrive the genomic features
        annotatedBSJs <- annotateBSJs(mergedBSJunctions, gtf)
        expect_is(annotatedBSJs, "data.frame")

        annotatedBSJsColNames <- .getAnnotatedBSJsColNames()
        expect_equal(colnames(annotatedBSJs), annotatedBSJsColNames)

        expect_equal(nrow(annotatedBSJs), 11) # there is 1 antisense circRNA
    }
)



test_that(
    "annotateBSJs() generates a data frame with the correct content",
    {
        gtf <- formatGTF(pathToGTF = "gencodeVM16.gtf")

        # Create the backSplicedJunctions data frame
        backSplicedJunctions <- getBackSplicedJunctions(gtf)
        mergedBSJunctions <-
            mergeBSJunctions(backSplicedJunctions, gtf)

        # Retrive the genomic geatures
        annotatedBSJs <- annotateBSJs(mergedBSJunctions, gtf)

        idColStartUpIntron <-
            which(colnames(annotatedBSJs) == "startUpIntron")
        idColEndDownIntron <-
            which(colnames(annotatedBSJs) == "endDownIntron")

        # positive strand
        annotatedBSJs_pos <-
            annotatedBSJs[annotatedBSJs$strand == "+",]
        for (i in 1:nrow(annotatedBSJs_pos)) {
            # take unique coordinates so that also sigle exon circRNA can be tested
            annotatedBSJsUnique_pos <- unique(as.numeric(annotatedBSJs_pos[i, c(idColStartUpIntron:idColEndDownIntron)]))
            ck_pos <-
                !is.unsorted(annotatedBSJsUnique_pos, na.rm = TRUE)
            expect_true(ck_pos)
        }

        # negative strand
        annotatedBSJs_neg <-
            annotatedBSJs[annotatedBSJs$strand == "-",]

        for (i in 1:nrow(annotatedBSJs_neg)) {
            # take unique coordinates so that also sigle exon circRNA can be tested
            annotatedBSJsUnique_neg <-
                unique(as.numeric(annotatedBSJs_neg[i, c(idColStartUpIntron:idColEndDownIntron)]))

            ck_neg <-
                !is.unsorted(rev(annotatedBSJsUnique_neg), na.rm = TRUE)
            expect_true(ck_neg)
        }

    }
)


test_that(
    "annotateBSJs() generates the correct data structure
    when using random back-spliced junctions",
    {
        gtf <- formatGTF(pathToGTF = "gencodeVM16.gtf")

        randomBSJunctions <-
            getRandomBSJunctions (gtf, n = 6, f = 10)

        annotatedBSJs <-
            annotateBSJs(randomBSJunctions, gtf, isRandom = TRUE)
        expect_is(annotatedBSJs, "data.frame")


        annotatedBSJsColNames <- .getAnnotatedBSJsColNames()
        expect_equal(colnames(annotatedBSJs), annotatedBSJsColNames)

        expect_equal(nrow(annotatedBSJs), nrow(randomBSJunctions))
    }
)



test_that(
    "annotateBSJs() generates a data frame with the correct content when using
    random back-spliced junctions",
    {
        gtf <- formatGTF(pathToGTF = "gencodeVM16.gtf")
        randomBSJunctions <-
            getRandomBSJunctions (gtf, n = 6, f = 10)
        annotatedBSJs <-
            annotateBSJs(randomBSJunctions, gtf, isRandom = TRUE)

        idColStartUpIntron <-
            which(colnames(annotatedBSJs) == "startUpIntron")
        idColEndDownIntron <-
            which(colnames(annotatedBSJs) == "endDownIntron")

        # positive strand
        annotatedBSJs_pos <-
            annotatedBSJs[annotatedBSJs$strand == "+",]


        for (i in 1:nrow(annotatedBSJs_pos)) {
            # take unique coordinates so that also sigle exon circRNA can be tested
            annotatedBSJsUnique_pos <-
                unique(as.numeric(annotatedBSJs_pos[i, c(idColStartUpIntron:idColEndDownIntron)]))
            ck_pos <-
                !is.unsorted(annotatedBSJsUnique_pos, na.rm = TRUE)
            expect_true(ck_pos)

        }

        # negative strand
        annotatedBSJs_neg <-
            annotatedBSJs[annotatedBSJs$strand == "-",]

        for (i in 1:nrow(annotatedBSJs_neg)) {
            # take unique coordinates so that also sigle exon circRNA can be tested
            annotatedBSJsUnique_neg <-
                unique(as.numeric(annotatedBSJs_neg[i, c(idColStartUpIntron:idColEndDownIntron)]))

            ck_neg <-
                !is.unsorted(rev(annotatedBSJsUnique_neg), na.rm = TRUE)
            expect_true(ck_neg)
        }

    }
)

