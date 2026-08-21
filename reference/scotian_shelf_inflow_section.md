# Endpoints of the named sections

Exposed so the geometry behind each index can be inspected, plotted, or
adjusted, rather than sitting as a constant inside a function body.
Neither reproduces a specific published section.

## Usage

``` r
scotian_shelf_inflow_section()

northeast_channel_section()
```

## Value

a list with `from` and `to`, each `c(longitude, latitude)`

## How these were placed

Not by eye. Both were measured against GLORYS monthly surface velocities
for January to April 2010 on two diagnostics:

- **capture fraction**, `|mean flow . n|` over `|mean flow|`, which is 1
  when the current crosses the section squarely and 0 when it runs along
  it

- **endpoint ratio**, the larger endpoint normal velocity over the peak
  along the section, which is near zero when the section spans the
  current rather than cutting through it

An earlier pair chosen from a map scored 0.65 and 0.27 on capture. The
Northeast Channel one ran nearly along the channel axis, so its
transport was the difference between two opposing halves. The current
pair score 0.99 and 0.93, with Channel endpoints at about 4% of the
peak.

The Northeast Channel line was then cross-checked against ETOPO depth,
since bathymetry is what defines that channel: along it, depth runs
roughly 120 m, 250 m, 80 m, so it starts on one bank, crosses the deep
water, and ends on the other.

Both normals point westerly, into the Gulf, which was verified
explicitly rather than inferred from the endpoint order.

`docs/methods.md` records the full derivation, and
`docs/section-placement-diagnostics.R` re-runs it on any `uo`/`vo`
extract.

## See also

[`section_transport()`](https://camilleross.org/derivoce/reference/section_transport.md)
for an arbitrary line

## Examples

``` r
scotian_shelf_inflow_section()
#> $from
#> [1] -66.15  43.54
#> 
#> $to
#> [1] -65.95  42.76
#> 
northeast_channel_section()
#> $from
#> [1] -66.19  42.72
#> 
#> $to
#> [1] -66.61  41.88
#> 
```
