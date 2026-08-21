# Time steps contributing to an integration window

Time steps contributing to an integration window

## Usage

``` r
contributing_steps(steps, i, window)
```

## Arguments

- steps:

  a time-step table from
  [`time_steps()`](https://camilleross.org/derivoce/reference/time_steps.md)

- i:

  index of the step being computed

- window:

  `"year"`, `"all"`, or a number of trailing steps

## Value

integer indices into `steps`, inclusive of `i`
