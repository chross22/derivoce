# Label runs of consecutive exceeding time steps, per cell

Works in step index rather than row order, so a cell absent from a step
breaks its run instead of being bridged by whatever row happens to come
next.

## Usage

``` r
event_runs(env_dat, exceeds, min_steps = 1L, max_gap = 0L)
```

## Arguments

- env_dat:

  an `sf` POINT object

- exceeds:

  logical, one per row

- min_steps:

  shortest run kept

- max_gap:

  non-exceeding steps that may be bridged

## Value

a list with `id` (0 outside an event) and `duration`
