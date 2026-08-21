# Warn about covariates a given operation cannot say anything about

These are the cases that run to completion and return a column of
perfectly valid numbers carrying no information: the lag of a static
covariate is the covariate, the horizontal gradient of a basin-wide
index is zero. Nothing downstream can tell such a column from a real
one, and a model handed it will happily report it as uninformative
rather than as a mistake — so the warning is at the point where the
intent is still visible.

## Usage

``` r
warn_degenerate(env_dat, vars, kind)
```

## Arguments

- env_dat:

  an `sf` POINT object

- vars:

  resolved covariate names

- kind:

  `"spatial"`, `"temporal"`, or `"any"` to skip the check

## Value

invisibly `NULL`; called for the warning

## Details

A warning rather than an error because the computation is well defined
and the caller may have meant it. `vars = NULL` expanding over an
enriched object is the case this is really aimed at.
