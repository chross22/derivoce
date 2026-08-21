# Warn when a Lagrangian diagnostic returns almost nothing

A field that is nearly all `NA` is the expected result of asking for a
longer trajectory than the domain can hold, and every individual `NA` is
correct. What is not acceptable is returning that silently: the caller
sees an empty column and cannot tell an unsuitable domain from a broken
function.

## Usage

``` r
warn_lagrangian(fn, na_fraction, detail, advice, threshold = 0.9)
```

## Arguments

- fn:

  name of the calling function, for the message

- na_fraction:

  fraction of the output that is `NA`

- detail:

  one or more sentences diagnosing the cause

- advice:

  what the caller can change

- threshold:

  warn at or above this fraction

## Value

invisibly `NULL`

## Details

The threshold is deliberately high. Losing a margin to the domain edge
is normal and does not need saying every time. Losing essentially
everything does.
