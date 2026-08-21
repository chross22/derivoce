# Join event runs separated by a short interruption

Join event runs separated by a short interruption

## Usage

``` r
bridge_gaps(flag, present, max_gap)
```

## Arguments

- flag:

  logical vector in time order

- present:

  the step index of each element

- max_gap:

  longest interruption that may be bridged

## Value

the flag vector with short gaps set `TRUE`
