# Decompose a single cell's series

Decompose a single cell's series

## Usage

``` r
decompose_one(y, t, month, degree)
```

## Arguments

- y:

  the covariate values for one cell

- t:

  elapsed days for those rows

- month:

  calendar month for those rows

- degree:

  polynomial degree of the trend

## Value

a list of components, plus counters for the warnings
