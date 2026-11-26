# Fit a linear model to diversity data

Fits a linear model ([`lm`](https://rdrr.io/r/stats/lm.html)) relating a
diversity metric to one or more predictors. When `group_col` is supplied
it is included as a fixed effect (a full mixed-effects model would
require the nlme dependency).

## Usage

``` r
fitDiversityModel(div_df, formula, group_col = NULL)
```

## Arguments

- div_df:

  A `data.frame` containing the diversity metrics and predictor columns
  (e.g.\\ output of `as.data.frame(calcDiversity(...))`).

- formula:

  A formula or character string describing the model (e.g.\\
  `shannon ~ timepoint`).

- group_col:

  Character scalar (optional). Column name in `div_df` representing a
  grouping variable. Added as a fixed effect if not already present in
  `formula`.

## Value

A list with components:

- model:

  The fitted `lm` object.

- summary_df:

  A `data.frame` of coefficient estimates, standard errors, t-values,
  and p-values.

- residuals:

  Numeric vector of model residuals.

## Examples

``` r
df <- data.frame(
    shannon = c(1.2, 1.5, 0.8, 1.1, 1.9, 0.7),
    timepoint = c(1, 2, 3, 1, 2, 3),
    group = c("A", "A", "A", "B", "B", "B")
)
res <- fitDiversityModel(df, shannon ~ timepoint)
res_grouped <- fitDiversityModel(df, shannon ~ timepoint,
    group_col = "group")
```
