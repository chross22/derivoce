# Hobday category from how far past the threshold a value reached

The categories are multiples of the threshold's own excess over the
climatology, so a place with little variability reaches a high category
on a smaller absolute departure than a variable one.

## Usage

``` r
mhw_category(excess, margin)
```

## Arguments

- excess:

  departure from the climatology, sign-corrected

- margin:

  the threshold's departure from the climatology, sign-corrected

## Value

an integer vector, 1 to 4
