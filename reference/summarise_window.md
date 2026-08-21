# Apply one summary to a window, ignoring missing values

[`min()`](https://rdrr.io/r/base/Extremes.html) and
[`max()`](https://rdrr.io/r/base/Extremes.html) of an all-missing window
would be `Inf` with a warning rather than `NA`, and
[`sd()`](https://rdrr.io/r/stats/sd.html) of a single value is `NA`
already, so the empty case is handled once here.

## Usage

``` r
summarise_window(z, stat)
```

## Arguments

- z:

  the values in the window

- stat:

  the statistic to apply

## Value

a single number
