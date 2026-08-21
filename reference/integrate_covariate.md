# Time-integrated covariate

Accumulates a covariate over preceding time steps at each location. A
survey does not sample the food available that instant so much as the
food that has built up since the season began, which is what the
original pipeline's `int_chl` captured for chlorophyll.

## Usage

``` r
integrate_covariate(env_dat, vars = NULL, window = "year", suffix = "_int")
```

## Arguments

- env_dat:

  an `sf` POINT object with one row per location and time step, as
  datamatch's access functions return

- vars:

  covariate columns to integrate; `NULL` does all of them

- window:

  `"year"` to accumulate from the start of each calendar year, `"all"`
  to accumulate over the whole record, or a positive integer for a
  rolling window of that many time steps

- suffix:

  suffix for the new columns

## Value

`env_dat` with an integrated column per covariate

## Details

The default `window = "year"` reproduces `int_chl` of Ross et al.
(2023), which integrated chlorophyll from January: a running sum from
January of each year, reset at the year boundary. A numeric window
instead sums over that many trailing time steps, giving a rolling total
that does not reset.

## References

Ross C, Runge J, Roberts J, Brady D, Tupper B, Record N (2023).
Estimating North Atlantic right whale prey based on Calanus finmarchicus
thresholds. *Marine Ecology Progress Series* **703**, 1-16.
[doi:10.3354/meps14204](https://doi.org/10.3354/meps14204)

## Examples

``` r
if (FALSE) { # \dontrun{
env <- integrate_covariate(env, "CHL")               # int_chl, Ross et al. 2023
env <- integrate_covariate(env, "CHL", window = 3)   # trailing 3-step total
} # }
```
