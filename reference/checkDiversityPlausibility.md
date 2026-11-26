# Check diversity estimates against heuristic plausibility ranges

Compares observed diversity metrics against wide plausibility ranges for
common pathogens, as a check on units and orders of magnitude. The
ranges are heuristic: they are the package's own, set generously enough
that a correctly computed estimate on a typical data set falls inside
them, and they are not taken from the literature or from any study. Each
is listed in Details. A value inside its range is therefore no evidence
that an estimate is correct, and one outside it most often points to a
unit, depth or QC problem rather than to biology.

## Usage

``` r
checkDiversityPlausibility(
  div_df,
  pathogen = c("sars_cov_2", "influenza", "hiv", "tb", "generic")
)
```

## Arguments

- div_df:

  A `DataFrame` or `data.frame` with columns `metric` and `value`, or
  containing diversity columns such as `pi` and `shannon` (as returned
  by [`calcDiversity`](calcDiversity.md)), with an optional `sample_id`.

- pathogen:

  Character scalar. Pathogen whose plausibility range is used. One of
  `"sars_cov_2"`, `"influenza"`, `"hiv"`, `"tb"`, or `"generic"`.

## Value

A `DataFrame` with one row per sample and metric and columns:

- sample_id:

  Sample identifier, taken from `div_df` when it has that column,
  otherwise the row number; `NA` for long input.

- metric:

  Diversity metric name.

- observed:

  Observed value.

- expected_min:

  Lower bound of the plausibility range.

- expected_max:

  Upper bound of the plausibility range.

- in_range:

  Logical; whether the observed value falls within the range, and
  `FALSE` when the value is missing.

## Details

The ranges, per site for \\\pi\\ and on the scale
[`calcDiversity`](calcDiversity.md) returns for Shannon entropy:

|              |              |             |
|--------------|--------------|-------------|
| **pathogen** | **pi**       | **shannon** |
| `sars_cov_2` | 1e-6 to 1e-3 | 0 to 15     |
| `influenza`  | 1e-5 to 5e-3 | 0 to 20     |
| `hiv`        | 1e-4 to 5e-2 | 0 to 30     |
| `tb`         | 1e-9 to 1e-4 | 0 to 5      |
| `generic`    | 0 to 1       | 0 to 100    |

A wide table, as [`calcDiversity`](calcDiversity.md) returns, is checked
row by row, one row per sample and metric; a long table with `metric`
and `value` columns is checked as given. Metrics without a range are
left out.

## Examples

``` r
div <- S4Vectors::DataFrame(sample_id = c("a", "b"),
    pi = c(1e-4, 0.3), shannon = c(5.0, 7.5))
checkDiversityPlausibility(div, pathogen = "sars_cov_2")
#> DataFrame with 4 rows and 6 columns
#>     sample_id      metric  observed expected_min expected_max  in_range
#>   <character> <character> <numeric>    <numeric>    <numeric> <logical>
#> 1           a          pi    0.0001        1e-06        0.001      TRUE
#> 2           b          pi    0.3000        1e-06        0.001     FALSE
#> 3           a     shannon    5.0000        0e+00       15.000      TRUE
#> 4           b     shannon    7.5000        0e+00       15.000      TRUE
```
