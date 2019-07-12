context("Test that getRegexPattern() function works correctly")
test_that("getRegexPattern() generates the correct regex pattern", {
    iupacSeq <- "ATCGTWWRRDHS"
    regexPattern <- getRegexPattern(iupacSeq, isDNA = TRUE)
    expect_equal(regexPattern, "ATCGT[AT][AT][AG][AG][AGT][ACT][GC]")

    iupacSeq <- "AUCGUWWRRDHS"
    regexPattern <- getRegexPattern(iupacSeq, isDNA = FALSE)
    expect_equal(regexPattern, "AUCGU[AU][AU][AG][AG][AGU][ACU][GC]")

})
