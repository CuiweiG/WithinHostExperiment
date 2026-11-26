# Export donor-recipient frequency table

For each transmission pair, creates a table of variant frequencies in
the donor and recipient.

## Usage

``` r
exportPairFrequencies(whe, pairId = NULL, usePassedOnly = TRUE)
```

## Arguments

- whe:

  A [`WithinHostExperiment`](WithinHostExperiment.md) with `role` and
  `pair_id` in `colData`.

- pairId:

  Character (optional). Specific pair to export. If NULL, exports all
  pairs.

- usePassedOnly:

  Logical. Use only QC-passed variants (default: TRUE).

## Value

A `DataFrame` with columns: chrom, position, ref, alt, donor_freq,
recipient_freq, shared (logical) and pair_id. A host carries the allele
at a site when its frequency there is known and above 0; `shared` is
`TRUE` where both hosts carry it. A recorded frequency of 0, such as a
site that a table lists without alternative reads, is an assessed site
without the allele: it is kept as 0 in its frequency column, and is not
shared. A frequency is `NA` where the host has no call, or its call
failed QC. Sites that neither host carries are left out.

## Examples

``` r
vcf1 <- system.file("extdata", "test_donor.vcf",
    package = "WithinHostExperiment")
vcf2 <- system.file("extdata", "test_recipient.vcf",
    package = "WithinHostExperiment")
whe <- readWithinHost(c(vcf1, vcf2),
    colData = S4Vectors::DataFrame(
        sample_id = c("donor", "recipient"),
        role = c("donor", "recipient"),
        pair_id = c("p1", "p1")),
    caller = "ivar")
exportPairFrequencies(whe)
#> DataFrame with 20 rows and 8 columns
#>            chrom  position         ref         alt donor_freq recipient_freq
#>      <character> <integer> <character> <character>  <numeric>      <numeric>
#> 1   test_segment        45           C           T       0.35           0.28
#> 2   test_segment        78           A           C         NA           0.07
#> 3   test_segment       112           A           T       0.08             NA
#> 4   test_segment       155           A           G         NA           0.15
#> 5   test_segment       198           A           G       0.22             NA
#> ...          ...       ...         ...         ...        ...            ...
#> 16  test_segment       701           T           G       0.15             NA
#> 17  test_segment       768           G           A       0.27             NA
#> 18  test_segment       812           G           T       0.04             NA
#> 19  test_segment       889           A           C       0.38            0.3
#> 20  test_segment       945           G           A       0.11             NA
#>        shared     pair_id
#>     <logical> <character>
#> 1        TRUE          p1
#> 2       FALSE          p1
#> 3       FALSE          p1
#> 4       FALSE          p1
#> 5       FALSE          p1
#> ...       ...         ...
#> 16      FALSE          p1
#> 17      FALSE          p1
#> 18      FALSE          p1
#> 19       TRUE          p1
#> 20      FALSE          p1
```
