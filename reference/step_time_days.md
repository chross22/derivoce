# Numeric time of each step, in days, including the hour when there is one

`HOUR` arrives with sub-daily fetches –
`accessCCMP(frequency = "6hourly")` and
`accessHYCOM(frequency = "3hourly")` – and is absent otherwise. Every
place that turns the time columns into a number goes through here, so
that hourly steps are distinct instants rather than duplicates of their
day: duplicate times would collapse steps, and in the Lagrangian
functions they would corrupt the interpolation between velocity fields.

## Usage

``` r
step_time_days(tbl)
```

## Arguments

- tbl:

  a data frame carrying the time columns

## Value

numeric days, with an hour fraction when `HOUR` is present
