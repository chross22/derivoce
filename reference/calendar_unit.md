# The calendar unit a time function counts in

Every function that looks back in time counts in calendar units, never
in positions in the series, and none has a default unit. Positions are
what a gap silently corrupts, and a default unit would be right for one
cadence and wrong for another, so the caller must say which.

## Usage

``` r
calendar_unit(by, fn)
```

## Arguments

- by:

  the value the caller supplied, or `NULL` if none was given

- fn:

  name of the calling function, for the message

## Value

`"hour"`, `"day"`, `"month"` or `"year"`
