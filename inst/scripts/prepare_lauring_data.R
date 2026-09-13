#!/usr/bin/env Rscript
## Prepare real SARS-CoV-2 iSNV data from Lauring Lab
## Source: github.com/lauringlab/SARS-CoV-2_VOC_transmission_bottleneck
## Paper: Bendall et al. 2023 Nature Communications 14:272
##
## This script converts the merged replicate format into
## per-sample iVar-style TSV files for WithinHostExperiment.

datadir <- "inst/scripts/real_data"
vars <- read.delim(file.path(datadir, "all_variants_filtered.tsv"),
                    stringsAsFactors = FALSE, quote = "\"")
pairs <- read.csv(file.path(datadir, "Transmission_pairs.csv"),
                   stringsAsFactors = FALSE)

## Select household HH6 (MHM1143 <-> MHM1184) - has 4 shared
## variants, real replicates, genuine transmission pair
## from a published Nature Communications study.
pair_info <- pairs[pairs$pair_id == "HH6_A", ]
donor_id <- pair_info$sample[pair_info$Transmission_indiv == "A"]
recip_id <- pair_info$sample[pair_info$Transmission_indiv == "B"]

cat("Selected pair: HH6_A\n")
cat("  Donor:", donor_id, "\n")
cat("  Recipient:", recip_id, "\n")

## Also get all variants from HH17 household for richer demo
hh17_samples <- unique(pairs$sample[grepl("HH17", pairs$pair_id)])
cat("  HH17 samples:", paste(hh17_samples, collapse = ", "), "\n")

## Get all unique samples we want
all_samples <- unique(c(donor_id, recip_id, hh17_samples))

## For each sample, extract variants and write iVar-format TSV
## The Lauring data has replicate 1 (_1) and replicate 2 (_2)
## We write replicate 1 as the main sample, replicate 2 as _rep2
write_ivar_tsv <- function(sample_id, vars_df, outdir, suffix = "") {
    sv <- vars_df[vars_df$sample == sample_id, ]
    if (nrow(sv) == 0) {
        cat("  No variants for", sample_id, "\n")
        return(invisible(NULL))
    }

    ## Replicate 1
    df1 <- data.frame(
        REGION = sv$REGION,
        POS = sv$POS,
        REF = sv$REF,
        ALT = sv$ALT,
        REF_DP = sv$REF_DP_1,
        ALT_DP = sv$ALT_DP_1,
        ALT_FREQ = sv$ALT_FREQ_1,
        TOTAL_DP = sv$TOTAL_DP_1,
        PVAL = sv$PVAL_1,
        PASS = sv$PASS_1,
        GFF_FEATURE = sv$GFF_FEATURE,
        REF_AA = sv$REF_AA,
        ALT_AA = sv$ALT_AA,
        stringsAsFactors = FALSE
    )
    fn1 <- paste0(sample_id, suffix, "_rep1.tsv")
    write.table(df1, file.path(outdir, fn1),
                sep = "\t", row.names = FALSE, quote = FALSE)

    ## Replicate 2
    df2 <- data.frame(
        REGION = sv$REGION,
        POS = sv$POS,
        REF = sv$REF,
        ALT = sv$ALT,
        REF_DP = sv$REF_DP_2,
        ALT_DP = sv$ALT_DP_2,
        ALT_FREQ = sv$ALT_FREQ_2,
        TOTAL_DP = sv$TOTAL_DP_2,
        PVAL = sv$PVAL_2,
        PASS = sv$PASS_2,
        GFF_FEATURE = sv$GFF_FEATURE,
        REF_AA = sv$REF_AA,
        ALT_AA = sv$ALT_AA,
        stringsAsFactors = FALSE
    )
    fn2 <- paste0(sample_id, suffix, "_rep2.tsv")
    write.table(df2, file.path(outdir, fn2),
                sep = "\t", row.names = FALSE, quote = FALSE)

    cat("  Written:", fn1, "(", nrow(df1), "variants),",
        fn2, "(", nrow(df2), "variants)\n")
}

## Write per-sample TSV files
for (s in all_samples) {
    write_ivar_tsv(s, vars, datadir)
}

## Summary
cat("\nData source: Lauring Lab, Nature Communications 2023\n")
cat("DOI: 10.1038/s41467-023-36001-5\n")
cat("Samples written:", length(all_samples), "\n")
cat("Files per sample: 2 (replicate 1 + replicate 2)\n")
