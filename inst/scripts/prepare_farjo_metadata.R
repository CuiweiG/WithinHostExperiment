## Sample metadata for the Farjo et al. (2024) participant used in Figure 4
##
## Reads three files from the study's repository at a fixed commit and writes
## inst/scripts/real_data/farjo_longitudinal/samples_432870.csv with, for each
## of the nine saliva samples: collection date, day of infection, mean genome
## coverage, and whether the study analysed the sample (mean coverage of at
## least 1000x, the low-coverage rule in the study's naivevariants.R).
##
## Source: https://github.com/BROOKELAB/SARS-CoV-2-within-host-evolution
## (MIT licence), commit fe34c8f5a8b08adeed102a0a854c5123258458eb.
## Run from the package root. Requires readxl.

commit <- "fe34c8f5a8b08adeed102a0a854c5123258458eb"
base <- paste0("https://raw.githubusercontent.com/BROOKELAB/",
               "SARS-CoV-2-within-host-evolution/", commit, "/metadata/")
participant <- "432870"
out <- file.path("inst", "scripts", "real_data", "farjo_longitudinal",
                 paste0("samples_", participant, ".csv"))

fetch <- function(name) {
    path <- file.path(tempdir(), name)
    utils::download.file(paste0(base, name), path, mode = "wb", quiet = TRUE)
    path
}

## BioSample sheet: sample, collection date (Excel serial) and source
biosample <- readxl::read_excel(fetch("SARS-CoV-2_samples.xlsx"),
                                col_names = FALSE, col_types = "text",
                                .name_repair = "minimal")
biosample <- as.data.frame(biosample)
header <- which(apply(biosample, 1, function(r) any(r %in% "*sample_name")))[1]
names(biosample) <- unlist(biosample[header, ])
biosample <- biosample[-seq_len(header), ]
depth <- utils::read.csv(fetch("naive_depth_table.csv"),
                         colClasses = "character")
depth <- depth[depth$participant == paste0("user_", participant), ]

samples <- biosample[biosample[["*sample_name"]] %in% depth$sample, ]
meta <- data.frame(
    sample = samples[["*sample_name"]],
    isolation_source = samples[["*isolation_source"]],
    collection_date = as.Date(as.numeric(samples[["*collection_date"]]),
                              origin = "1899-12-30"),
    stringsAsFactors = FALSE
)
meta$mean_coverage <- as.numeric(depth$mean_coverage[match(meta$sample,
                                                           depth$sample)])
meta <- meta[order(meta$collection_date), ]
## Day of infection as the study defines it: day 1 is the first sample
meta$day_of_infection <- as.integer(meta$collection_date -
                                    min(meta$collection_date)) + 1L
meta$analysed_by_study <- meta$mean_coverage >= 1000

## Cross-check against the study's own day and inclusion table
saliva <- as.data.frame(readxl::read_excel(fetch("all_saliva_user_info.xlsx"),
                                           col_types = "text"))
saliva <- saliva[saliva$user_id == participant, ]
listed <- match(saliva$sample_barcode, meta$sample)
stopifnot(
    nrow(meta) == 9L,
    all(meta$isolation_source == "saliva"),
    all(!is.na(listed)),
    identical(sort(saliva$sample_barcode), sort(meta$sample[meta$analysed_by_study])),
    all(as.integer(saliva$day_of_infection) == meta$day_of_infection[listed])
)

utils::write.csv(meta, out, row.names = FALSE)
message("Wrote ", out)
print(meta, row.names = FALSE)
