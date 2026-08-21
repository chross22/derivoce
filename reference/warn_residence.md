# Report how much of the result is censored or unknown

Report how much of the result is censored or unknown

## Usage

``` r
warn_residence(released, censored, escaped, max_days)
```

## Arguments

- released:

  particles released

- censored:

  particles still inside when the window ended

- escaped:

  particles that left the velocity field without leaving the box

- max_days:

  the window

## Value

invisible `NULL`
