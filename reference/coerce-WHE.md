# Convert between WithinHostExperiment, data.frame, VRanges and VCF

`as(x, "data.frame")` gives one row per site and sample with an observed
frequency; `as(x, "VRanges")` gives the same rows as a
[`VRanges`](https://rdrr.io/pkg/VariantAnnotation/man/VRanges-class.html).
In the other direction, `as(vr, "WithinHostExperiment")` builds the
object from the `refDepth`, `altDepth` and `totalDepth` of a `VRanges`,
one sample per `sampleNames` level, and
`as(vcf, "WithinHostExperiment")` from a VCF read with
[`VariantAnnotation::readVcf()`](https://rdrr.io/pkg/VariantAnnotation/man/readVcf-methods.html),
choosing the field layout from the header as
[`readWithinHost`](readWithinHost.md)`(caller = "auto")` does and taking
every sample the file holds. Coercion to and from a `VRanges` keeps the
whole `seqinfo`: every recorded sequence and its length, including
sequences without a call, since the genome length that
[`calcDiversity`](calcDiversity.md) divides by sums over all of them.
The VCF coercion keeps the lengths of the file's `##contig` lines.

## Value

A `data.frame`, a
[`VRanges`](https://rdrr.io/pkg/VariantAnnotation/man/VRanges-class.html)
or a [`WithinHostExperiment`](WithinHostExperiment.md).

## See also

[`readWithinHost`](readWithinHost.md) for VCF files.

## Examples

``` r
data(example_whe)
df <- as(example_whe, "data.frame")
head(df)
#>          chrom position ref alt sample_id altFreq totalDepth refCount altCount
#> 1 test_segment       45   C   T     donor    0.35       1500      975      525
#> 2 test_segment      112   A   T     donor    0.08        800      736       64
#> 3 test_segment      198   A   G     donor    0.22       2200     1716      484
#> 4 test_segment      267   C   A     donor    0.45       3100     1705     1395
#> 5 test_segment      334   G   T     donor    0.12       1200     1056      144
#> 6 test_segment      401   C   T     donor    0.06        500      470       30
#>   qcPass
#> 1   TRUE
#> 2   TRUE
#> 3   TRUE
#> 4   TRUE
#> 5   TRUE
#> 6   TRUE

vr <- as(example_whe, "VRanges")
vr
#> VRanges object with 33 ranges and 0 metadata columns:
#>            seqnames    ranges strand         ref              alt
#>               <Rle> <IRanges>  <Rle> <character> <characterOrRle>
#>    [1] test_segment        45      *           C                T
#>    [2] test_segment       112      *           A                T
#>    [3] test_segment       198      *           A                G
#>    [4] test_segment       267      *           C                A
#>    [5] test_segment       334      *           G                T
#>    ...          ...       ...    ...         ...              ...
#>   [29] test_segment       367      *           T                A
#>   [30] test_segment       489      *           T                G
#>   [31] test_segment       611      *           C                A
#>   [32] test_segment       733      *           A                G
#>   [33] test_segment       856      *           C                T
#>            totalDepth       refDepth       altDepth   sampleNames
#>        <integerOrRle> <integerOrRle> <integerOrRle> <factorOrRle>
#>    [1]           1500            975            525         donor
#>    [2]            800            736             64         donor
#>    [3]           2200           1716            484         donor
#>    [4]           3100           1705           1395         donor
#>    [5]           1200           1056            144         donor
#>    ...            ...            ...            ...           ...
#>   [29]           2200           1474            726   independent
#>   [30]           1500           1215            285   independent
#>   [31]           2800           1652           1148   independent
#>   [32]            700            658             42   independent
#>   [33]           1600           1360            240   independent
#>        softFilterMatrix
#>                <matrix>
#>    [1]                 
#>    [2]                 
#>    [3]                 
#>    [4]                 
#>    [5]                 
#>    ...              ...
#>   [29]                 
#>   [30]                 
#>   [31]                 
#>   [32]                 
#>   [33]                 
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
#>   hardFilters: NULL
as(vr, "WithinHostExperiment")
#> WithinHostExperiment with 28 variant sites and 3 samples
#>   assays(5): altFreq, totalDepth, altCount, refCount, qcPass
#>   qcLog: 0 filtering step(s)
#>   transmission pairs: 0 

vcf <- VariantAnnotation::readVcf(system.file("extdata",
    "test_donor.vcf", package = "WithinHostExperiment"))
as(vcf, "WithinHostExperiment")
#> WithinHostExperiment with 15 variant sites and 1 samples
#>   assays(5): altFreq, totalDepth, refCount, altCount, qcPass
#>   qcLog: 0 filtering step(s)
#>   transmission pairs: 0 
```
