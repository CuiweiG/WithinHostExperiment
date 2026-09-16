#!/usr/bin/env Rscript
# ============================================================
# Generate synthetic test data for WithinHostExperiment
# Simulates iSNVs on a 1000bp viral genome segment
# ============================================================

outdir <- "inst/extdata"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

set.seed(42) # set.seed allowed only in data generation scripts
bases <- c("A", "T", "G", "C")

# ---- 1. Reference genome ----
ref_seq <- paste0(sample(bases, 1000, replace = TRUE), collapse = "")
writeLines(c(">test_segment", ref_seq), file.path(outdir, "test_ref.fa"))
cat("Written: test_ref.fa\n")

# ---- 2. Helper: pick a different base ----
alt_base <- function(ref_char) {
    sample(setdiff(bases, ref_char), 1)
}

# ---- 3. Write iVar-format VCF ----
write_ivar_vcf <- function(positions, refs, alts, aafs, depths,
                           sample_name, filename) {
    header <- c(
        "##fileformat=VCFv4.2",
        "##source=iVar",
        "##reference=test_ref.fa",
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=REF_DP,Number=1,Type=Integer,Description="Ref depth">',
        '##FORMAT=<ID=ALT_DP,Number=1,Type=Integer,Description="Alt depth">',
        '##FORMAT=<ID=ALT_FREQ,Number=1,Type=Float,Description="Alt frequency">',
        '##FORMAT=<ID=ALT_QUAL,Number=1,Type=Float,Description="Alt quality">',
        paste0("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t",
               sample_name)
    )
    alt_dp <- as.integer(round(depths * aafs))
    ref_dp <- depths - alt_dp
    records <- vapply(seq_along(positions), function(i) {
        paste("test_segment", positions[i], ".", refs[i], alts[i],
              "255", "PASS",
              paste0("DP=", depths[i]),
              "GT:REF_DP:ALT_DP:ALT_FREQ:ALT_QUAL",
              paste0("0/1:", ref_dp[i], ":", alt_dp[i], ":",
                     sprintf("%.6f", aafs[i]), ":30"),
              sep = "\t")
    }, character(1))
    writeLines(c(header, records), file.path(outdir, filename))
    cat("Written:", filename, "(", length(positions), "variants )\n")
}

# ---- 4. Donor: 15 iSNVs ----
donor_pos <- c(45, 112, 198, 267, 334, 401, 456, 523, 589,
               634, 701, 768, 812, 889, 945)
donor_ref <- vapply(donor_pos, function(p) substr(ref_seq, p, p), character(1))
donor_alt <- vapply(donor_ref, alt_base, character(1))
donor_aaf <- c(0.35, 0.08, 0.22, 0.45, 0.12, 0.06, 0.18, 0.41,
               0.09, 0.33, 0.15, 0.27, 0.04, 0.38, 0.11)
donor_depth <- c(1500L, 800L, 2200L, 3100L, 1200L, 500L, 1800L, 2500L,
                 900L, 2000L, 1600L, 2800L, 150L, 3500L, 1100L)
# Site 13 (pos=812) has depth 150 -- used to test depth filter

write_ivar_vcf(donor_pos, donor_ref, donor_alt, donor_aaf, donor_depth,
               "donor", "test_donor.vcf")

# ---- 5. Recipient: 10 iSNVs (5 shared) ----
shared_idx <- c(1, 4, 8, 10, 14) # indices in donor
recip_unique_pos <- c(78, 155, 290, 478, 650)
recip_pos <- c(donor_pos[shared_idx], recip_unique_pos)
recip_ref <- vapply(recip_pos, function(p) substr(ref_seq, p, p), character(1))
recip_alt <- c(
    donor_alt[shared_idx],
    vapply(recip_ref[6:10], alt_base, character(1))
)
recip_aaf <- c(0.28, 0.39, 0.35, 0.41, 0.30,
               0.07, 0.15, 0.22, 0.10, 0.05)
recip_depth <- c(2000L, 2800L, 2100L, 1900L, 3200L,
                 600L, 1400L, 2500L, 900L, 400L)

write_ivar_vcf(recip_pos, recip_ref, recip_alt, recip_aaf, recip_depth,
               "recipient", "test_recipient.vcf")

# ---- 6. Independent: 8 iSNVs ----
indep_pos <- c(33, 127, 245, 367, 489, 611, 733, 856)
indep_ref <- vapply(indep_pos, function(p) substr(ref_seq, p, p), character(1))
indep_alt <- vapply(indep_ref, alt_base, character(1))
indep_aaf <- c(0.12, 0.25, 0.08, 0.33, 0.19, 0.41, 0.06, 0.15)
indep_depth <- c(1200L, 1800L, 900L, 2200L, 1500L, 2800L, 700L, 1600L)

write_ivar_vcf(indep_pos, indep_ref, indep_alt, indep_aaf, indep_depth,
               "independent", "test_independent.vcf")

# ---- 7. Technical replicate Rep1: donor + noise + 2 large deviations + 1 unique site ----
rep1_pos <- c(donor_pos, 500L) # 16 sites (15 + 1 rep-only)
rep1_ref <- c(donor_ref, substr(ref_seq, 500, 500))
rep1_alt <- c(donor_alt, alt_base(substr(ref_seq, 500, 500)))
rep1_aaf <- c(
    donor_aaf[1:2],
    donor_aaf[3] + 0.15, # pos=198: offset 0.15 (0.22 -> 0.37) -- discordant
    donor_aaf[4:6],
    donor_aaf[7] - 0.12, # pos=456: offset 0.12 (0.18 -> 0.06) -- discordant
    donor_aaf[8:15],
    0.04 # pos=500: only in rep1 -- rep-only
)
rep1_aaf <- pmax(rep1_aaf, 0.01) # ensure non-negative
rep1_depth <- c(donor_depth, 1000L)

write_ivar_vcf(rep1_pos, rep1_ref, rep1_alt, rep1_aaf, rep1_depth,
               "donor_rep1", "test_donor_rep1.vcf")

# ---- 8. Technical replicate Rep2: donor + different noise (pos=500 absent) ----
rep2_aaf <- donor_aaf + runif(15, -0.02, 0.02)
rep2_aaf <- pmax(pmin(rep2_aaf, 0.49), 0.02)
rep2_depth <- donor_depth + as.integer(runif(15, -100, 100))
rep2_depth <- pmax(rep2_depth, 100L)

write_ivar_vcf(donor_pos, donor_ref, donor_alt, rep2_aaf, rep2_depth,
               "donor_rep2", "test_donor_rep2.vcf")
# Note: rep2 has only 15 sites; pos=500 is absent

# ---- 9. iVar TSV format ----
write_ivar_tsv <- function(positions, refs, alts, aafs, depths, filename) {
    alt_dp <- as.integer(round(depths * aafs))
    ref_dp <- depths - alt_dp
    df <- data.frame(
        REGION = "test_segment",
        POS = positions,
        REF = refs,
        ALT = alts,
        REF_DP = ref_dp,
        ALT_DP = alt_dp,
        ALT_FREQ = aafs,
        TOTAL_DP = depths,
        PVAL = 0.001,
        PASS = "TRUE",
        stringsAsFactors = FALSE
    )
    write.table(df, file.path(outdir, filename),
                sep = "\t", row.names = FALSE, quote = FALSE)
    cat("Written:", filename, "\n")
}
write_ivar_tsv(donor_pos, donor_ref, donor_alt, donor_aaf, donor_depth,
               "test_donor_ivar.tsv")

# ---- 10. Generic CSV ----
generic_df <- data.frame(
    chrom = "test_segment",
    pos = donor_pos,
    ref = donor_ref,
    alt = donor_alt,
    aaf = donor_aaf,
    depth = donor_depth,
    stringsAsFactors = FALSE
)
write.csv(generic_df, file.path(outdir, "test_generic.csv"), row.names = FALSE)
cat("Written: test_generic.csv\n")

# ---- 11. Sample metadata ----
coldata <- data.frame(
    sample_id = c("donor", "recipient", "independent",
                  "donor_rep1", "donor_rep2"),
    host_id = c("P1", "P2", "P3", "P1", "P1"),
    replicate_group = c(NA, NA, NA, "donor_reps", "donor_reps"),
    replicate_id = c(NA, NA, NA, "rep1", "rep2"),
    role = c("donor", "recipient", "independent", "donor", "donor"),
    pair_id = c("pair_1", "pair_1", NA, "pair_1", "pair_1"),
    ct_value = c(18.5, 22.1, 15.3, 18.7, 18.4),
    stringsAsFactors = FALSE
)
write.csv(coldata, file.path(outdir, "test_coldata.csv"), row.names = FALSE)
cat("Written: test_coldata.csv\n")

# ---- 12. LoFreq-format VCF ----
write_lofreq_vcf <- function(positions, refs, alts, aafs, depths,
                             sample_name, filename) {
    header <- c(
        "##fileformat=VCFv4.2",
        "##source=lofreq call",
        "##reference=test_ref.fa",
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw Depth">',
        '##INFO=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">',
        '##INFO=<ID=SB,Number=1,Type=Integer,Description="Phred-scaled strand bias">',
        '##INFO=<ID=DP4,Number=4,Type=Integer,Description="fwd-ref,rev-ref,fwd-alt,rev-alt">',
        '##FILTER=<ID=PASS,Description="All filters passed">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    )
    alt_dp <- as.integer(round(depths * aafs))
    ref_dp <- depths - alt_dp
    fwd_ref <- as.integer(round(ref_dp * 0.5))
    rev_ref <- ref_dp - fwd_ref
    fwd_alt <- as.integer(round(alt_dp * 0.5))
    rev_alt <- alt_dp - fwd_alt

    records <- vapply(seq_along(positions), function(i) {
        info_str <- paste0(
            "DP=", depths[i],
            ";AF=", sprintf("%.6f", aafs[i]),
            ";SB=0",
            ";DP4=", fwd_ref[i], ",", rev_ref[i], ",",
            fwd_alt[i], ",", rev_alt[i]
        )
        paste("test_segment", positions[i], ".", refs[i], alts[i],
              "255", "PASS", info_str, sep = "\t")
    }, character(1))
    writeLines(c(header, records), file.path(outdir, filename))
    cat("Written:", filename, "(LoFreq format,", length(positions), "variants)\n")
}

write_lofreq_vcf(donor_pos, donor_ref, donor_alt, donor_aaf, donor_depth,
                 "donor_lofreq", "test_donor_lofreq.vcf")

# ---- 13. Freebayes-format VCF ----
write_freebayes_vcf <- function(positions, refs, alts, aafs, depths,
                                sample_name, filename) {
    header <- c(
        "##fileformat=VCFv4.2",
        "##source=freeBayes v1.3.6",
        "##reference=test_ref.fa",
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">',
        '##FORMAT=<ID=RO,Number=1,Type=Integer,Description="Ref observations">',
        '##FORMAT=<ID=AO,Number=1,Type=Integer,Description="Alt observations">',
        paste0("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t",
               sample_name)
    )
    alt_dp <- as.integer(round(depths * aafs))
    ref_dp <- depths - alt_dp
    records <- vapply(seq_along(positions), function(i) {
        paste("test_segment", positions[i], ".", refs[i], alts[i],
              "255", "PASS",
              paste0("DP=", depths[i]),
              "GT:DP:RO:AO",
              paste0("0/1:", depths[i], ":", ref_dp[i], ":", alt_dp[i]),
              sep = "\t")
    }, character(1))
    writeLines(c(header, records), file.path(outdir, filename))
    cat("Written:", filename, "(Freebayes format,", length(positions), "variants)\n")
}

write_freebayes_vcf(donor_pos, donor_ref, donor_alt, donor_aaf, donor_depth,
                    "donor_freebayes", "test_donor_freebayes.vcf")

# ---- 14. Primer BED file ----
primer_bed <- data.frame(
    chrom = "test_segment",
    start = c(40, 880),
    end   = c(65, 950),
    name  = c("primer_F1", "primer_R1"),
    stringsAsFactors = FALSE
)
write.table(primer_bed, file.path(outdir, "test_primers.bed"),
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
cat("Written: test_primers.bed\n")

cat("\n=== All test data generated ===\n")
