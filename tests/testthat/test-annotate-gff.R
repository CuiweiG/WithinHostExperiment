# tests/testthat/test-annotate-gff.R

test_that("annotateFromGFF adds GFF_FEATURE to rowRanges", {
    ## Create a minimal GFF3
    gff_file <- tempfile(fileext = ".gff3")
    writeLines(c(
        "##gff-version 3",
        "test_segment\t.\tgene\t1\t500\t.\t+\t.\tName=geneA",
        "test_segment\t.\tgene\t501\t1000\t.\t+\t.\tName=geneB"
    ), gff_file)

    vcf <- system.file("extdata", "test_donor.vcf",
        package = "WithinHostExperiment")
    whe <- readWithinHost(vcf,
        colData = S4Vectors::DataFrame(sample_id = "donor"),
        caller = "ivar")

    whe2 <- annotateFromGFF(whe, gff_file)
    rr <- SummarizedExperiment::rowRanges(whe2)
    expect_true("GFF_FEATURE" %in% colnames(S4Vectors::mcols(rr)))
    expect_true(any(!is.na(S4Vectors::mcols(rr)$GFF_FEATURE)))

    unlink(gff_file)
})

test_that("annotateFromGFF rejects missing file", {
    gr <- GenomicRanges::GRanges("seg",
        IRanges::IRanges(1:3, width = 1))
    whe <- WithinHostExperiment(
        assays = list(altFreq = matrix(runif(3), ncol = 1)),
        rowRanges = gr,
        colData = S4Vectors::DataFrame(sample_id = "S1"))
    expect_error(annotateFromGFF(whe, "/no/such/file.gff"),
                 "not found")
})

test_that("annotateFromGFF supports CDS feature type", {
    gff_file <- tempfile(fileext = ".gff3")
    writeLines(c(
        "##gff-version 3",
        "test_segment\t.\tCDS\t100\t900\t.\t+\t.\tName=orf1"
    ), gff_file)

    vcf <- system.file("extdata", "test_donor.vcf",
        package = "WithinHostExperiment")
    whe <- readWithinHost(vcf,
        colData = S4Vectors::DataFrame(sample_id = "donor"),
        caller = "ivar")

    whe2 <- annotateFromGFF(whe, gff_file, featureType = "CDS")
    rr <- SummarizedExperiment::rowRanges(whe2)
    expect_true("GFF_FEATURE" %in% colnames(S4Vectors::mcols(rr)))

    unlink(gff_file)
})
