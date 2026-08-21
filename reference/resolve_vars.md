# Check that requested covariates exist and can be operated on

Three checks, deliberately of different severities:

## Usage

``` r
resolve_vars(env_dat, vars, kind = c("any", "spatial", "temporal"))
```

## Arguments

- env_dat:

  an `sf` POINT object

- vars:

  requested covariate names, or `NULL` for all of them

- kind:

  the kind of operation the caller is about to perform, so the
  degeneracy relevant to it can be checked: `"temporal"`, `"spatial"`,
  or `"any"` (the default) for operations that are neither, such as a
  column-wise difference

## Value

the resolved covariate names

## Details

A **missing** column is an error: nothing can be done.

A **non-numeric** column is an error when named explicitly — a factor
cannot be differentiated or summed, so the request cannot be honoured —
but is skipped silently when `vars = NULL` swept it up. A caller who did
not name `CHL_source` did not mean it, and failing the whole call over a
column they never asked for would make the `NULL` default unusable on
any object that has been through
[`datamatch::fill_satellite_gaps()`](https://camilleross.org/datamatch/reference/fill_satellite_gaps.html).

A **degenerate** column — static where the operation is temporal,
spatially uniform where it is spatial — is only a warning, because the
computation is well defined and the caller may have meant it. See
[`warn_degenerate()`](https://camilleross.org/derivoce/reference/warn_degenerate.md).
