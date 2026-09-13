#!/usr/bin/env Rscript
## generate_readme_figures.R
##
## Regenerates the six README figures from the public data in
## inst/scripts/real_data, using the package functions throughout.
##
## Data
##   Bendall et al. (2023) Nat. Commun. 14:272 -- all_variants_filtered.tsv,
##     AvgCoverage.all, Transmission_pairs.csv
##     (github.com/lauringlab/SARS-CoV-2_VOC_transmission_bottleneck)
##   Farjo et al. (2024) J. Virol. 98:e01618-23 -- farjo_longitudinal/
##   NCBI RefSeq NC_045512.2 -- sars2_NC045512.gff3, sars2_NC045512.fasta
##
## Run from the package root of a clone:
##   Rscript inst/scripts/generate_readme_figures.R
##
## Outputs
##   man/figures/fig1..fig6 *.png            README figures (300 dpi)
##   readme_figure_output/*.pdf              vector versions (cairo_pdf)
##   readme_figure_output/readme_figure_values.csv
##                                           every number quoted in the captions
##   readme_figure_output/input_md5.csv, sessionInfo.txt

suppressPackageStartupMessages({
    library(ggplot2)
    library(patchwork)
    library(S4Vectors)
    library(SummarizedExperiment)
    library(GenomicRanges)
})
stopifnot(packageVersion("ggplot2") >= "3.5.0")
if (!file.exists("DESCRIPTION") ||
    read.dcf("DESCRIPTION", fields = "Package")[1, 1] != "WithinHostExperiment") {
    stop("Run this script from the WithinHostExperiment package root.")
}
pkgload::load_all(".", quiet = TRUE)

DATA <- file.path("inst", "scripts", "real_data")
OUT <- "readme_figure_output"
FIGDIR <- file.path("man", "figures")
dir.create(OUT, showWarnings = FALSE)
dir.create(FIGDIR, showWarnings = FALSE)

inputs <- c(
    variants = file.path(DATA, "all_variants_filtered.tsv"),
    coverage = file.path(DATA, "AvgCoverage.all"),
    pairs = file.path(DATA, "Transmission_pairs.csv"),
    gff = file.path(DATA, "sars2_NC045512.gff3"),
    fasta = file.path(DATA, "sars2_NC045512.fasta"),
    farjo = file.path(DATA, "farjo_longitudinal")
)
missing_inputs <- inputs[!file.exists(inputs)]
if (length(missing_inputs)) {
    stop("Missing input data: ", paste(missing_inputs, collapse = ", "))
}
farjo_files <- sort(list.files(inputs[["farjo"]], pattern = "ivar",
                               full.names = TRUE))
md5_files <- c(inputs[setdiff(names(inputs), "farjo")], farjo_files)
write.csv(data.frame(file = md5_files, md5 = unname(tools::md5sum(md5_files))),
          file.path(OUT, "input_md5.csv"), row.names = FALSE)

values <- data.frame(figure = character(), panel = character(),
                     quantity = character(), value = character(),
                     stringsAsFactors = FALSE)
record <- function(figure, panel, quantity, value) {
    values[nrow(values) + 1L, ] <<- list(figure, panel, quantity,
                                         format(value, digits = 6))
    invisible(value)
}

## Okabe-Ito colours
pal <- list(blue = "#0072B2", vermillion = "#D55E00", green = "#009E73",
            orange = "#E69F00", skyblue = "#56B4E9", purple = "#CC79A7",
            grey = "#999999", black = "#000000")

theme_pub <- function(bs = 10) {
    theme_classic(base_size = bs, base_family = "sans") +
        theme(axis.title = element_text(size = bs, colour = "black"),
              axis.text = element_text(size = bs - 1, colour = "black"),
              axis.line = element_line(linewidth = 0.5, colour = "black"),
              axis.ticks = element_line(linewidth = 0.35, colour = "black"),
              axis.ticks.length = unit(2.5, "pt"),
              legend.title = element_text(size = bs - 1, face = "bold"),
              legend.text = element_text(size = bs - 1),
              legend.key.size = unit(10, "pt"),
              legend.background = element_blank(),
              panel.background = element_rect(fill = "white", colour = NA),
              panel.grid = element_blank(),
              strip.background = element_blank(),
              strip.text = element_text(size = bs, face = "bold"),
              plot.background = element_rect(fill = "white", colour = NA),
              plot.margin = margin(6, 10, 6, 6),
              plot.tag = element_text(size = bs + 3, face = "bold",
                                      colour = "black"))
}
inside <- function(x, y, hjust = x, vjust = y) {
    theme(legend.position = "inside",
          legend.position.inside = c(x, y),
          legend.justification.inside = c(hjust, vjust))
}
save_figure <- function(plot, name, width, height) {
    ggsave(file.path(FIGDIR, paste0(name, ".png")), plot, width = width,
           height = height, dpi = 300, bg = "white")
    ggsave(file.path(OUT, paste0(name, ".pdf")), plot, width = width,
           height = height, device = grDevices::cairo_pdf)
    message("saved ", name)
}

GL <- 29903L
## The iVar calls are on MN908947.3, which is the same sequence as NC_045512.2
SEQNAME <- "MN908947.3"

## ---- Bendall et al. (2023) data ----
vars <- read.delim(inputs[["variants"]], stringsAsFactors = FALSE)
cov_lines <- readLines(inputs[["coverage"]])[-1L]
cov_lines <- cov_lines[nchar(trimws(cov_lines)) > 0L]
cov_fields <- strsplit(trimws(cov_lines), "\\s+")
depth_by_sample <- tapply(as.numeric(vapply(cov_fields, `[`, "", 2L)),
                          sub("_[12]$", "", vapply(cov_fields, `[`, "", 1L)),
                          mean)
median_depth <- median(depth_by_sample)
pairs <- read.csv(inputs[["pairs"]], stringsAsFactors = FALSE,
                  fileEncoding = "UTF-8-BOM")

vars$freq_diff <- abs(vars$ALT_FREQ_1 - vars$ALT_FREQ_2)
vars$concordant <- vars$freq_diff <= 0.02
vars$mean_freq <- (vars$ALT_FREQ_1 + vars$ALT_FREQ_2) / 2
concordant <- vars[vars$concordant, ]
samples <- unique(vars$sample)
n_isnv <- nrow(vars)
n_samples <- length(samples)
n_conc <- sum(vars$concordant)
n_disc <- sum(!vars$concordant)
r2 <- cor(vars$ALT_FREQ_1, vars$ALT_FREQ_2)^2
record("data", "", "iSNV calls", n_isnv)
record("data", "", "samples", n_samples)
record("data", "", "transmission pairs in metadata", length(unique(pairs$pair_id)))

sample_depth <- function(s) if (s %in% names(depth_by_sample)) depth_by_sample[[s]] else NULL
pi_sample <- function(freq, s) {
    d <- sample_depth(s)
    piISNV(freq, GL, meanDepth = if (!is.null(d) && d > 1) d else NULL)
}
div <- do.call(rbind, lapply(samples, function(s) {
    data.frame(sample = s,
               pi_naive = pi_sample(vars$ALT_FREQ_1[vars$sample == s], s) * 1e4,
               pi_qc = pi_sample(concordant$ALT_FREQ_1[concordant$sample == s], s) * 1e4,
               stringsAsFactors = FALSE)
}))
wt_pi <- wilcox.test(div$pi_naive, div$pi_qc, paired = TRUE,
                     alternative = "greater", exact = FALSE)
med_naive <- median(div$pi_naive)
med_qc <- median(div$pi_qc)

## ================================================================
## Figure 1: replicate QC
## ================================================================
record("fig1", "a", "concordant calls (|freq diff| <= 0.02)", n_conc)
record("fig1", "a", "discordant calls", n_disc)
record("fig1", "a", "replicate R2", r2)
record("fig1", "b", "median pi naive (x1e-4)", med_naive)
record("fig1", "b", "median pi QC (x1e-4)", med_qc)
record("fig1", "b", "paired one-sided Wilcoxon p (naive > QC)", wt_pi$p.value)
record("fig1", "b", "samples", n_samples)

rep_df <- data.frame(r1 = vars$ALT_FREQ_1 * 100, r2 = vars$ALT_FREQ_2 * 100,
                     status = factor(ifelse(vars$concordant, "Concordant", "Discordant"),
                                     levels = c("Concordant", "Discordant")))
p1a <- ggplot(rep_df, aes(r1, r2, fill = status)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = pal$grey, linewidth = 0.4) +
    geom_point(shape = 21, size = 2.2, stroke = 0.3, alpha = 0.85, colour = "white") +
    annotate("text", x = 97, y = 5, hjust = 1, vjust = 0,
             label = sprintf("R² = %.3f", r2), size = 3.2) +
    scale_fill_manual(name = NULL, values = c(Concordant = pal$green, Discordant = pal$vermillion),
                      labels = c(sprintf("Concordant (n = %d)", n_conc),
                                 sprintf("Discordant (n = %d)", n_disc))) +
    scale_x_continuous("Replicate 1 frequency (%)", limits = c(0, 100), breaks = seq(0, 100, 25)) +
    scale_y_continuous("Replicate 2 frequency (%)", limits = c(0, 100), breaks = seq(0, 100, 25)) +
    coord_equal() + labs(tag = "a") +
    guides(fill = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    theme_pub(10) + inside(0.35, 0.98, 0, 1) +
    theme(legend.background = element_rect(fill = alpha("white", 0.92), colour = NA))

pi_long <- data.frame(sample = rep(div$sample, each = 2L), x = rep(c(1, 2), times = n_samples),
                      pi = c(rbind(div$pi_naive, div$pi_qc)),
                      condition = rep(c("Naive", "QC"), times = n_samples))
p1b <- ggplot(pi_long, aes(x = x, y = pi)) +
    geom_line(aes(group = sample), colour = "grey60", linewidth = 0.4, alpha = 0.45) +
    geom_point(aes(colour = condition), size = 1.8, alpha = 0.85, shape = 16,
               position = position_jitter(width = 0.06, height = 0, seed = 42)) +
    annotate("point", x = c(1, 2), y = c(med_naive, med_qc), shape = 18, size = 5.5, colour = pal$black) +
    annotate("text", x = 2.45, y = max(div$pi_naive) * 0.95, hjust = 1, vjust = 1, size = 3.0,
             label = sprintf("Wilcoxon p = %.1e\nn = %d", wt_pi$p.value, n_samples)) +
    scale_colour_manual(values = c(Naive = pal$orange, QC = pal$green), guide = "none") +
    scale_x_continuous(NULL, breaks = c(1, 2), labels = c("Naive", "QC"), limits = c(0.55, 2.55)) +
    scale_y_continuous(expression(pi ~ "(" * 10^{-4} * ")"), expand = expansion(mult = c(0.03, 0.06))) +
    labs(tag = "b") + theme_pub(10)

thresholds <- seq(0, 0.20, by = 0.01)
median_pi_by_threshold <- vapply(thresholds, function(t) {
    median(vapply(samples, function(s) {
        pi_sample(vars$ALT_FREQ_1[vars$sample == s & vars$freq_diff <= t], s) * 1e4
    }, numeric(1)))
}, numeric(1))
thr_df <- data.frame(t = thresholds * 100, pi = median_pi_by_threshold)
p1c <- ggplot(thr_df, aes(t, pi)) +
    geom_ribbon(aes(ymin = 0, ymax = pi), fill = pal$blue, alpha = 0.20) +
    geom_line(colour = pal$blue, linewidth = 0.9) + geom_point(size = 1.5, colour = pal$blue) +
    geom_vline(xintercept = c(2, 5), linetype = "dashed", colour = c(pal$vermillion, pal$orange), linewidth = 0.45) +
    scale_x_continuous("|Frequency difference| threshold (%)", breaks = seq(0, 20, 5),
                       expand = expansion(mult = c(0.01, 0.03))) +
    scale_y_continuous(expression("Median" ~ pi ~ "(" * 10^{-4} * ")"), expand = expansion(mult = c(0.02, 0.10))) +
    labs(tag = "c") + theme_pub(10)
save_figure((p1a | p1b | p1c), "fig1_replicate_qc", 7.5, 3.1)

## ================================================================
## Figure 2: frequency spectrum
## ================================================================
n_sites <- length(unique(paste(vars$POS, vars$REF, vars$ALT)))
record("fig2", "", "distinct sites", n_sites)
vars$status <- factor(ifelse(vars$concordant, "Concordant", "Discordant"),
                      levels = c("Concordant", "Discordant"))
record("fig2", "", "discordant calls below 10% mean frequency", sum(!vars$concordant & vars$mean_freq < 0.10))
p2 <- ggplot(vars, aes(x = mean_freq * 100, fill = status)) +
    geom_histogram(binwidth = 2.5, colour = "white", linewidth = 0.2, alpha = 0.90, boundary = 0) +
    geom_vline(xintercept = 3, linetype = "dashed", colour = pal$black, linewidth = 0.55) +
    annotate("text", x = 4.5, y = Inf, hjust = 0, vjust = 1.6, label = "3%", size = 3.2) +
    scale_fill_manual(name = NULL, values = c(Concordant = pal$blue, Discordant = pal$vermillion),
                      labels = c(sprintf("Concordant (n = %d)", n_conc), sprintf("Discordant (n = %d)", n_disc))) +
    scale_x_continuous("Alternative allele frequency (%)", breaks = seq(0, 100, 20),
                       expand = expansion(mult = c(0.01, 0.02))) +
    scale_y_continuous("iSNVs", expand = expansion(mult = c(0, 0.10))) +
    labs(caption = sprintf("n = %d calls | %d sites | %d samples", n_isnv, n_sites, n_samples)) +
    theme_pub(10) + inside(0.98, 0.98) +
    theme(plot.caption = element_text(size = 7, colour = "grey40", margin = margin(t = 4)))
save_figure(p2, "fig2_frequency_spectrum", 4.2, 3.2)

## ================================================================
## Figure 3: depth, annotation, transmission
## ================================================================
record("fig3", "a", "median mean read depth", median_depth)
p3a <- ggplot(data.frame(d = as.numeric(depth_by_sample)), aes(x = d)) +
    geom_histogram(fill = pal$orange, colour = "white", linewidth = 0.2, binwidth = 300, alpha = 0.90, boundary = 0) +
    geom_vline(xintercept = median_depth, linetype = "dashed", colour = pal$black, linewidth = 0.5) +
    annotate("text", x = median_depth + 200, y = Inf, hjust = 0, vjust = 1.8, size = 3.0, lineheight = 1.15,
             label = sprintf("median\n%s×", format(round(median_depth), big.mark = ","))) +
    scale_x_continuous("Mean read depth", labels = scales::label_comma(), expand = expansion(mult = c(0.01, 0.05))) +
    scale_y_continuous("Samples", expand = expansion(mult = c(0, 0.15))) + labs(tag = "a") + theme_pub(10)

gff_tmp <- tempfile(fileext = ".gff3")
writeLines(gsub("NC_045512\\.2", SEQNAME, readLines(inputs[["gff"]])), gff_tmp)
fasta_tmp <- tempfile(fileext = ".fa")
fasta_lines <- readLines(inputs[["fasta"]])
fasta_lines[1L] <- paste0(">", SEQNAME)
writeLines(fasta_lines, fasta_tmp)

gr_conc <- GRanges(SEQNAME, IRanges::IRanges(concordant$POS, width = 1))
mcols(gr_conc)$ref <- concordant$REF
mcols(gr_conc)$alt <- concordant$ALT
whe_conc <- WithinHostExperiment(
    assays = list(altFreq = matrix(concordant$ALT_FREQ_1, ncol = 1, dimnames = list(NULL, "pooled"))),
    rowRanges = gr_conc, colData = DataFrame(sample_id = "pooled"))
whe_conc <- annotateFromGFF(whe_conc, gff_tmp)
gene_of_site <- mcols(rowRanges(whe_conc))$GFF_FEATURE
gene_of_site[is.na(gene_of_site)] <- "intergenic"
gene_table <- as.data.frame(table(Gene = gene_of_site), stringsAsFactors = FALSE)
gene_table <- gene_table[order(gene_table$Freq, decreasing = TRUE), ]
gene_table$Gene <- factor(gene_table$Gene, levels = rev(gene_table$Gene))
p3b <- ggplot(gene_table, aes(y = Gene, x = Freq)) +
    geom_col(fill = pal$blue, colour = "white", linewidth = 0.2, alpha = 0.90) +
    geom_text(aes(label = Freq), hjust = -0.25, size = 3.1) +
    scale_x_continuous("Concordant iSNVs", expand = expansion(mult = c(0, 0.22))) +
    labs(y = NULL, tag = "b") + theme_pub(10) + theme(axis.text.y = element_text(size = 9))

pair_df <- merge(pairs[pairs$Transmission_indiv == "A", c("sample", "pair_id")],
                 pairs[pairs$Transmission_indiv == "B", c("sample", "pair_id")],
                 by = "pair_id", suffixes = c("_donor", "_recipient"))
in_recipient <- function(recipient, pos, alt) {
    any(vars$sample == recipient & vars$POS == pos & vars$ALT == alt)
}
example_pair <- "HH46"
ep <- pair_df[pair_df$pair_id == example_pair, ]
donor_calls <- concordant[concordant$sample == ep$sample_donor[1], ]
loll <- data.frame(pos = donor_calls$POS,
                   gene = vapply(donor_calls$POS, function(p) {
                       g <- gene_of_site[concordant$POS == p][1]
                       if (is.na(g)) "intergenic" else g
                   }, character(1)),
                   donor_freq = donor_calls$ALT_FREQ_1 * 100,
                   detected = mapply(in_recipient, ep$sample_recipient[1], donor_calls$POS, donor_calls$ALT),
                   stringsAsFactors = FALSE)
loll$label <- factor(sprintf("%s (%s)", loll$pos, loll$gene), levels = rev(sprintf("%s (%s)", loll$pos, loll$gene)))
n_lost <- sum(!loll$detected)
record("fig3", "c", paste("donor iSNVs in pair", example_pair), nrow(loll))
record("fig3", "c", paste("donor iSNVs not detected in recipient, pair", example_pair), n_lost)
p3c <- ggplot(loll, aes(y = label, x = donor_freq, fill = detected)) +
    geom_col(width = 0.6, alpha = 0.90) +
    geom_text(aes(label = sprintf("%.1f%%", donor_freq)), hjust = -0.15, size = 2.8) +
    scale_fill_manual(name = NULL, values = c("TRUE" = pal$green, "FALSE" = pal$vermillion),
                      labels = c("TRUE" = "Detected", "FALSE" = "Not detected")) +
    scale_x_continuous("Donor frequency (%)", expand = expansion(mult = c(0, 0.30))) +
    annotate("text", x = max(loll$donor_freq) * 1.2, y = 0.6, hjust = 1, vjust = 0, size = 3.0, lineheight = 1.2,
             label = sprintf("Pair %s\n%d/%d not detected", example_pair, n_lost, nrow(loll))) +
    labs(y = NULL, tag = "c") + theme_pub(10) + inside(0.98, 0.02) +
    theme(axis.text.y = element_text(size = 8))

sharing <- do.call(rbind, lapply(seq_len(nrow(pair_df)), function(i) {
    dv <- concordant[concordant$sample == pair_df$sample_donor[i], ]
    if (nrow(dv) == 0L) return(NULL)
    shared <- sum(mapply(in_recipient, pair_df$sample_recipient[i], dv$POS, dv$ALT))
    data.frame(pair = pair_df$pair_id[i], donor = nrow(dv), shared = shared, lost = nrow(dv) - shared)
}))
sharing <- sharing[order(sharing$donor, decreasing = TRUE), ]
sharing$rank <- seq_len(nrow(sharing))
n_zero <- sum(sharing$shared == 0)
record("fig3", "d", "pairs with at least one concordant donor iSNV", nrow(sharing))
record("fig3", "d", "pairs sharing no donor iSNV", n_zero)
record("fig3", "d", "percentage sharing none", 100 * n_zero / nrow(sharing))
share_long <- rbind(data.frame(rank = sharing$rank, count = sharing$shared, type = "Shared"),
                    data.frame(rank = sharing$rank, count = sharing$lost, type = "Not detected"))
share_long$type <- factor(share_long$type, levels = c("Not detected", "Shared"))
p3d <- ggplot(share_long, aes(x = rank, y = count, fill = type)) +
    geom_col(width = 0.7, alpha = 0.90) +
    scale_fill_manual(name = NULL, values = c("Not detected" = pal$vermillion, Shared = pal$green)) +
    annotate("text", x = nrow(sharing) * 0.95, y = Inf, hjust = 1, vjust = 1.5, size = 3.0, lineheight = 1.2,
             label = sprintf("%d/%d (%.0f%%)\nshare none", n_zero, nrow(sharing), 100 * n_zero / nrow(sharing))) +
    scale_x_continuous("Transmission pair", expand = expansion(mult = c(0.01, 0.01))) +
    scale_y_continuous("Donor iSNVs", expand = expansion(mult = c(0, 0.12))) +
    labs(tag = "d") + theme_pub(10) + inside(0.02, 0.98)
save_figure((p3a + p3b) / (p3c + p3d), "fig3_qc_overview", 7.0, 5.8)

## ================================================================
## Figure 4: longitudinal sampling (Farjo et al. 2024)
## ================================================================
n_tp <- length(farjo_files)
whe_f <- readWithinHostTable(farjo_files, format = "ivar",
    colData = DataFrame(sample_id = paste0("t", seq_len(n_tp)), host_id = "p", timepoint = seq_len(n_tp)))
whe_f <- flagISNV(whe_f, ISNVFilter(minDepth = 100L, minFreq = 0.03, maxFreq = 0.97))
div_f <- as.data.frame(calcDiversity(whe_f, genomeLength = GL, indices = c("pi", "richness")))
div_f$timepoint <- seq_len(n_tp)
div_f$pi_e4 <- as.numeric(div_f$pi) * 1e4
div_f$richness <- as.numeric(div_f$richness)
peak <- div_f$timepoint[which.max(div_f$richness)]
record("fig4", "", "timepoints", n_tp)
record("fig4", "a", "peak QC-passed iSNVs", max(div_f$richness))
record("fig4", "a", "timepoint of peak", peak)
record("fig4", "a", "timepoint of maximum pi", div_f$timepoint[which.max(div_f$pi_e4)])
scale_factor <- max(div_f$richness) / max(div_f$pi_e4)
p4a <- ggplot(div_f, aes(x = timepoint)) +
    geom_col(aes(y = richness), fill = pal$skyblue, alpha = 0.50, width = 0.6) +
    geom_line(aes(y = pi_e4 * scale_factor), colour = pal$vermillion, linewidth = 1.0) +
    geom_point(aes(y = pi_e4 * scale_factor), colour = pal$vermillion, size = 2.5) +
    scale_y_continuous("QC-passed iSNVs (bars)",
                       sec.axis = sec_axis(~ . / scale_factor, name = expression(pi ~ "(" * 10^{-4} * ", line)")),
                       expand = expansion(mult = c(0, 0.08))) +
    scale_x_continuous("Timepoint", breaks = seq_len(n_tp)) + labs(tag = "a") + theme_pub(10) +
    theme(axis.title.y.right = element_text(colour = pal$vermillion),
          axis.text.y.right = element_text(colour = pal$vermillion))

traj <- as.data.frame(trackFrequency(whe_f, hostCol = "host_id", timeCol = "timepoint"))
freq_range <- tapply(traj$frequency, traj$variant_key, function(x) diff(range(x)))
top10 <- names(sort(freq_range, decreasing = TRUE))[seq_len(min(10L, length(freq_range)))]
p4b <- ggplot(traj[traj$variant_key %in% top10, ],
              aes(x = timepoint, y = frequency * 100, colour = variant_key, group = variant_key)) +
    geom_line(linewidth = 0.7, alpha = 0.80) + geom_point(size = 1.5, alpha = 0.90) +
    geom_hline(yintercept = 3, linetype = "dotted", colour = pal$grey, linewidth = 0.3) +
    scale_colour_manual(values = rep(unname(unlist(pal[c("blue", "vermillion", "green", "orange", "skyblue", "purple", "black", "grey")])), 2L)) +
    scale_x_continuous("Timepoint", breaks = seq_len(n_tp)) +
    scale_y_continuous("Alternative allele frequency (%)", limits = c(0, 100)) +
    labs(colour = NULL, tag = "b") + theme_pub(10) + theme(legend.position = "none")

key_tp <- unique(c(1L, peak, n_tp))
sfs_tp <- lapply(key_tp, function(d) buildSFS(whe_f, sampleIdx = d, fold = TRUE, genomeLength = GL, nBins = 8))
sfs_df <- do.call(rbind, lapply(seq_along(key_tp), function(k) {
    s <- sfs_tp[[k]]; br <- s@breaks
    data.frame(freq = (br[-length(br)] + br[-1L]) / 2 * 100, prop = s@counts / max(sum(s@counts), 1),
               tp = paste("Timepoint", key_tp[k]), stringsAsFactors = FALSE)
}))
sfs_df$tp <- factor(sfs_df$tp, levels = unique(sfs_df$tp))
p4c <- ggplot(sfs_df, aes(x = freq, y = prop, fill = tp)) +
    geom_col(position = "dodge", alpha = 0.85, colour = "white", linewidth = 0.15) +
    scale_fill_manual(name = NULL, values = c(pal$skyblue, pal$vermillion, pal$green)[seq_along(key_tp)]) +
    scale_x_continuous("Minor allele frequency (%)") +
    scale_y_continuous("Proportion", expand = expansion(mult = c(0, 0.12))) +
    labs(tag = "c") + theme_pub(10) + inside(0.98, 0.98)

whe_ft <- flagTemporalInconsistency(whe_f, minTimepoints = 2L)
temporal_class <- mcols(rowRanges(whe_ft))$temporal_class
n_persistent <- sum(temporal_class == "persistent", na.rm = TRUE)
n_transient <- sum(temporal_class == "transient", na.rm = TRUE)
record("fig4", "d", "persistent sites (>= 2 timepoints)", n_persistent)
record("fig4", "d", "transient sites (1 timepoint)", n_transient)
record("fig4", "d", "percentage transient", 100 * n_transient / (n_persistent + n_transient))
tc_df <- data.frame(class = factor(c("Persistent\n(≥ 2 timepoints)", "Transient\n(1 timepoint)"),
                                   levels = c("Persistent\n(≥ 2 timepoints)", "Transient\n(1 timepoint)")),
                    count = c(n_persistent, n_transient))
p4d <- ggplot(tc_df, aes(x = class, y = count, fill = class)) +
    geom_col(width = 0.6, alpha = 0.90) +
    geom_text(aes(label = count), vjust = -0.3, size = 3.5, fontface = "bold") +
    scale_fill_manual(values = c(pal$green, pal$vermillion), guide = "none") +
    scale_y_continuous("Variant sites", expand = expansion(mult = c(0, 0.15))) +
    labs(x = NULL, tag = "d") + theme_pub(10)
save_figure((p4a + p4b) / (p4c + p4d), "fig4_diversity_landscape", 7.0, 5.8)

## ================================================================
## Figure 5: population genetics
## ================================================================
tajima_sample <- function(freq, depth, cap = 100L) {
    S <- length(freq)
    if (S == 0L) return(NA_real_)
    n <- min(max(as.integer(round(depth)), 2L), cap)
    tajimaD(S, n, piISNV(freq, GL, meanDepth = n), GL)
}
taj <- do.call(rbind, lapply(samples, function(s) {
    d <- sample_depth(s); if (is.null(d)) d <- 100
    data.frame(sample = s,
               D_qc = tajima_sample(concordant$ALT_FREQ_1[concordant$sample == s], d),
               D_naive = tajima_sample(vars$ALT_FREQ_1[vars$sample == s], d),
               S_qc = sum(concordant$sample == s))
}))
taj_qc <- taj[taj$S_qc > 0 & is.finite(taj$D_qc), ]
record("fig5", "a", "samples with >= 1 concordant iSNV", nrow(taj_qc))
record("fig5", "a", "median Tajima D (QC, n capped at 100)", median(taj_qc$D_qc))
record("fig5", "a", "percentage of samples with D < 0", 100 * mean(taj_qc$D_qc < 0))
p5a <- ggplot(taj_qc, aes(x = D_qc)) +
    geom_histogram(bins = 15, fill = pal$blue, colour = "white", linewidth = 0.2, alpha = 0.90) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = pal$grey, linewidth = 0.5) +
    geom_vline(xintercept = median(taj_qc$D_qc), colour = pal$vermillion, linewidth = 0.6) +
    annotate("text", x = max(taj_qc$D_qc), y = Inf, vjust = 1.5, hjust = 1, size = 3.0, lineheight = 1.2,
             colour = pal$vermillion,
             label = sprintf("median = %.2f\n%.0f%% < 0", median(taj_qc$D_qc), 100 * mean(taj_qc$D_qc < 0))) +
    scale_x_continuous("Tajima's D (n capped at 100)") +
    scale_y_continuous("Samples", expand = expansion(mult = c(0, 0.12))) + labs(tag = "a") + theme_pub(10)

## b: Nei-Gojobori dN/dS, pooling the concordant iSNVs of all samples
whe_codon <- annotateCodonChange(whe_conc, gff_tmp, fasta_tmp)
dnds <- as.data.frame(dndsWithinHost(whe_codon, gff = gff_tmp, refFasta = fasta_tmp, usePassedOnly = FALSE))
record("fig5", "b", "pooled synonymous iSNVs", dnds$nS[1])
record("fig5", "b", "pooled nonsynonymous iSNVs", dnds$nN[1])
record("fig5", "b", "pooled dN/dS (Jukes-Cantor)", dnds$dNdS[1])
record("fig5", "b", "concordant iSNVs outside CDS or unannotated", nrow(concordant) - dnds$nS[1] - dnds$nN[1])
gene_dnds <- dnds[!is.na(dnds$gene), ]
for (i in seq_len(nrow(gene_dnds))) {
    record("fig5", "b", paste0("dN/dS ", gene_dnds$gene[i], " (", gene_dnds$gene_nN[i], "N/", gene_dnds$gene_nS[i], "S)"),
           gene_dnds$gene_dNdS[i])
}
gene_dnds <- gene_dnds[gene_dnds$gene_nS + gene_dnds$gene_nN >= 2L, ]
gene_dnds <- gene_dnds[order(gene_dnds$gene_nS + gene_dnds$gene_nN, decreasing = TRUE), ]
gene_dnds$gene <- factor(gene_dnds$gene, levels = rev(gene_dnds$gene))
gene_dnds$has_syn <- gene_dnds$gene_nS > 0
gene_dnds$bar <- ifelse(is.na(gene_dnds$gene_dNdS), 0, gene_dnds$gene_dNdS)
p5b <- ggplot(gene_dnds, aes(y = gene, x = bar, fill = has_syn)) +
    geom_col(width = 0.6, alpha = 0.90) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = pal$grey, linewidth = 0.5) +
    geom_text(aes(label = ifelse(has_syn, sprintf("%.2f (%dN/%dS)", gene_dNdS, gene_nN, gene_nS),
                                 sprintf("%dN, 0S", gene_nN))), hjust = -0.05, size = 2.7) +
    scale_fill_manual(values = c("TRUE" = pal$vermillion, "FALSE" = pal$grey), guide = "none") +
    scale_x_continuous("dN/dS (Nei-Gojobori, pooled)", expand = expansion(mult = c(0, 0.45))) +
    labs(y = NULL, tag = "b") + theme_pub(10) + theme(axis.text.y = element_text(size = 9))

ranked <- div[order(div$pi_naive, decreasing = TRUE), ]
ranked$rank <- seq_len(nrow(ranked))
ranked_long <- rbind(data.frame(rank = ranked$rank, pi = ranked$pi_naive, type = "Naive"),
                     data.frame(rank = ranked$rank, pi = ranked$pi_qc, type = "QC"))
p5c <- ggplot(ranked_long, aes(rank, pi, colour = type)) +
    geom_segment(data = ranked, aes(x = rank, xend = rank, y = pi_qc, yend = pi_naive),
                 colour = "grey65", linewidth = 0.5, inherit.aes = FALSE) +
    geom_point(size = 1.5, alpha = 0.90, shape = 16) +
    geom_hline(yintercept = c(med_naive, med_qc), linetype = "dotted", colour = c(pal$orange, pal$green), linewidth = 0.5) +
    scale_colour_manual(name = NULL, values = c(Naive = pal$orange, QC = pal$green),
                        guide = guide_legend(override.aes = list(size = 3))) +
    scale_x_continuous("Sample rank", expand = expansion(mult = c(0.01, 0.01))) +
    scale_y_continuous(expression(pi ~ "(" * 10^{-4} * ")"), expand = expansion(mult = c(0, 0.08))) +
    labs(tag = "c") + theme_pub(10) + inside(0.98, 0.98)

pooled_whe <- function(df, id) {
    gr <- GRanges(SEQNAME, IRanges::IRanges(df$POS, width = 1))
    mcols(gr)$ref <- df$REF; mcols(gr)$alt <- df$ALT
    WithinHostExperiment(assays = list(altFreq = matrix(df$ALT_FREQ_1, ncol = 1, dimnames = list(NULL, id)),
                                       totalDepth = matrix(as.integer(df$TOTAL_DP_1), ncol = 1, dimnames = list(NULL, id))),
                         rowRanges = gr, colData = DataFrame(sample_id = id))
}
sfs_naive <- buildSFS(pooled_whe(vars, "naive"), fold = TRUE, genomeLength = GL, nBins = 10)
sfs_qc <- buildSFS(pooled_whe(concordant, "qc"), fold = TRUE, genomeLength = GL, nBins = 10)
sfs_corrected <- correctSFSBias(sfs_qc, method = "binomial")
sfs_frame <- function(s, label) {
    br <- s@breaks
    data.frame(freq = (br[-length(br)] + br[-1L]) / 2 * 100, count = s@counts, label = label)
}
sfs5 <- rbind(sfs_frame(sfs_naive, sprintf("Naive (n = %d)", sfs_naive@nSites)),
              sfs_frame(sfs_qc, sprintf("QC (n = %d)", sfs_qc@nSites)),
              sfs_frame(sfs_corrected, sprintf("QC, bias-corrected (n = %d)", sfs_corrected@nSites)))
sfs5$label <- factor(sfs5$label, levels = unique(sfs5$label))
p5d <- ggplot(sfs5, aes(x = freq, y = count, fill = label)) +
    geom_col(position = "dodge", alpha = 0.85, width = 2.0, colour = "white", linewidth = 0.15) +
    scale_fill_manual(name = NULL, values = c(pal$orange, pal$green, pal$blue)) +
    scale_x_continuous("Minor allele frequency (%)") +
    scale_y_continuous("iSNVs", expand = expansion(mult = c(0, 0.12))) +
    labs(tag = "d") + theme_pub(10) + inside(0.98, 0.98) + theme(legend.text = element_text(size = 7.5))
save_figure((p5a + p5b) / (p5c + p5d), "fig5_evolutionary_analysis", 7.5, 5.8)

## ================================================================
## Figure 6: effect of QC on Tajima's D and the SFS
## ================================================================
taj_naive <- taj[is.finite(taj$D_naive), ]
wt_D <- wilcox.test(taj_naive$D_naive, taj_qc$D_qc, exact = FALSE)
record("fig6", "a", "median Tajima D naive", median(taj_naive$D_naive))
record("fig6", "a", "median Tajima D QC", median(taj_qc$D_qc))
record("fig6", "a", "Wilcoxon rank-sum p (naive vs QC)", wt_D$p.value)
d_long <- rbind(data.frame(D = taj_naive$D_naive, condition = "Naive"),
                data.frame(D = taj_qc$D_qc, condition = "QC"))
d_long$condition <- factor(d_long$condition, levels = c("Naive", "QC"))
p6a <- ggplot(d_long, aes(x = D, fill = condition)) +
    geom_histogram(bins = 14, colour = "white", linewidth = 0.2, alpha = 0.55, position = "identity") +
    geom_vline(xintercept = 0, linetype = "dashed", colour = pal$grey, linewidth = 0.4) +
    scale_fill_manual(name = NULL, values = c(Naive = pal$orange, QC = pal$green)) +
    annotate("text", x = max(d_long$D), y = Inf, hjust = 1, vjust = 1.5, size = 3.0,
             label = sprintf("Wilcoxon p = %.2f", wt_D$p.value)) +
    scale_x_continuous("Tajima's D") + scale_y_continuous("Samples", expand = expansion(mult = c(0, 0.15))) +
    labs(tag = "a") + theme_pub(10) + inside(0.02, 0.98)

comparison <- compareSFS(sfs_naive, sfs_qc)
record("fig6", "b", "compareSFS chi-squared", comparison$chisq_stat)
record("fig6", "b", "compareSFS p", comparison$chisq_p)
prop_df <- rbind(data.frame(freq = (sfs_naive@breaks[-length(sfs_naive@breaks)] + sfs_naive@breaks[-1L]) / 2 * 100,
                            prop = comparison$proportions1, label = "Naive"),
                 data.frame(freq = (sfs_qc@breaks[-length(sfs_qc@breaks)] + sfs_qc@breaks[-1L]) / 2 * 100,
                            prop = comparison$proportions2, label = "QC"))
prop_df$label <- factor(prop_df$label, levels = c("Naive", "QC"))
p6b <- ggplot(prop_df, aes(x = freq, y = prop, fill = label)) +
    geom_col(position = "dodge", alpha = 0.80, colour = "white", linewidth = 0.15) +
    scale_fill_manual(name = NULL, values = c(Naive = pal$orange, QC = pal$green)) +
    annotate("label", x = 40, y = max(prop_df$prop) * 0.70, size = 3.2, lineheight = 1.3,
             fill = alpha("white", 0.90), label.padding = unit(4, "pt"),
             label = sprintf("χ² = %.1f\np = %.3f", comparison$chisq_stat, comparison$chisq_p)) +
    scale_x_continuous("Minor allele frequency (%)") +
    scale_y_continuous("Proportion", expand = expansion(mult = c(0, 0.12))) +
    labs(tag = "b") + theme_pub(10) + inside(0.75, 0.85, 0, 1)
save_figure((p6a | p6b), "fig6_consensus_validation", 7.0, 3.2)

write.csv(values, file.path(OUT, "readme_figure_values.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
print(values, row.names = FALSE, right = FALSE)
message("All six figures saved. Values: ", file.path(OUT, "readme_figure_values.csv"))
