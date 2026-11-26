# WithinHostExperiment: Container for within-host pathogen variant data

Extends
[`RangedSummarizedExperiment-class`](https://rdrr.io/pkg/SummarizedExperiment/man/RangedSummarizedExperiment-class.html)
to hold intra-host single nucleotide variant (iSNV) data with integrated
quality control audit logging.

Subsetting preserves the `qcLog` and `qcFilters` slots. These slots are
experiment-level metadata (not row/column-dependent) and are carried
forward unchanged.

## Usage

``` r
# S4 method for class 'WithinHostExperiment'
show(object)

# S4 method for class 'WithinHostExperiment,ANY,ANY,ANY'
x[i, j, ..., drop = TRUE]
```

## Arguments

- object:

  A `WithinHostExperiment` object.

- x:

  A [`WithinHostExperiment`](WithinHostExperiment.md).

- i:

  Row indices (variant sites).

- j:

  Column indices (samples).

- ...:

  Additional arguments passed to the parent method.

- drop:

  Logical (ignored for `SummarizedExperiment`).

## Value

A [`WithinHostExperiment`](WithinHostExperiment.md).

## Slots

- `qcLog`:

  `DataFrame`. Audit trail of QC filtering steps.

- `qcFilters`:

  `list`. History of [`ISNVFilter`](ISNVFilter.md) objects that have
  been applied. Validity checks each one, so a filter saved by an
  earlier version of the package, or changed with `@<-` into an invalid
  one, is reported, and the object cannot be subset or combined until
  the filter is rebuilt with [`ISNVFilter`](ISNVFilter.md) and stored
  with `qcFilters<-`.

## Standard assays

- `altCount`:

  Integer matrix. Alternative allele read counts.

- `refCount`:

  Integer matrix. Reference allele read counts.

- `totalDepth`:

  Integer matrix. Total read depth.

- `altFreq`:

  Numeric matrix. Alt allele frequencies (0-1).

- `qcPass`:

  Logical matrix. TRUE = passed QC, FALSE = flagged. Initialized to all
  TRUE; updated by [`flagISNV`](flagISNV.md). No package function writes
  `NA`, but one set by hand counts as not passed wherever an analysis
  uses only passing calls.

Assays must be ordinary R matrices. Out-of-memory backends such as
`DelayedArray` or `HDF5Array` are not supported: the validity method
refuses a `qcPass` assay that is not a base logical matrix, and
[`flagISNV`](flagISNV.md) combines assays with base R operators.

## Expected colData columns

- `sample_id`:

  Character. Unique sample identifier (required).

- `host_id`:

  Character. Host/patient identifier (optional).

- `replicate_group`:

  Character. Technical replicate group identifier (optional).

- `replicate_id`:

  Character. Replicate number within group (optional).

- `role`:

  Character. Sample role: "donor", "recipient", or "independent"
  (optional).

- `pair_id`:

  Character. Transmission pair identifier (optional).

## See also

[`readWithinHost`](readWithinHost.md) and
[`readWithinHostTable`](readWithinHostTable.md) for constructing from
files. [`flagISNV`](flagISNV.md) for quality control.
[`passedISNV`](WHE-accessors.md) for extracting QC-passed variants.
[`transmissionPairs`](WHE-accessors.md) for extracting donor-recipient
pairs.

## Examples

``` r
data(example_whe)
example_whe
#> WithinHostExperiment with 28 variant sites and 3 samples
#>   assays(5): altFreq, totalDepth, refCount, altCount, qcPass
#>   qcLog: 0 filtering step(s)
#>   transmission pairs: 1 
```
