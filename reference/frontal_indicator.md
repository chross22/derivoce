# Whether each row's own cell is frontal

Shares its threshold logic with
[`distance_to_front()`](https://camilleross.org/derivoce/reference/distance_to_front.md),
so that a cell counted as frontal here is one that would have had a
distance of zero there.

## Usage

``` r
frontal_indicator(env_dat, var, threshold, quantile, scope, per)
```

## Arguments

- env_dat:

  an `sf` POINT object

- var:

  the covariate whose gradient defines a front

- threshold:

  explicit gradient cutoff, or `NULL`

- quantile:

  quantile used when `threshold` is `NULL`

- scope:

  `"record"` or `"step"`

- per:

  `"km"` or `"m"`

## Value

numeric vector of 1, 0 and `NA`
