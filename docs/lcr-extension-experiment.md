# Can the Labrador Current retroflection index be extended from Copernicus?

**Short answer: not from monthly fields.** This records an experiment that did not
work, and why, so nobody repeats it.

## Why it was worth trying

`datamatch::attach_climate_index()` serves `LCR`, the retroflection index of
Jutras et al. (2023). It covers **1993–2014 only** and cannot be attached to
recent observations, which is a real limitation for anything modelling the last
decade.

The index is derived rather than observed, so in principle it can be recomputed.
Jutras et al. seeded virtual particles weekly across a line on the Labrador Shelf
at (53°N, 56.7°W)–(54.3°N, 52.0°W), tracked each for three years through
GLORYS12V1 with [OceanParcels](https://oceanparcels.org/), and took the
difference between the counts reaching hydrographic sections on the Labrador and
Scotian shelves. High retroflection means fewer particles continue southwest.

derivoce already advects particles: `ftle()` and `fsle()` use fourth-order
Runge-Kutta through a time-varying field, with the same m/s to degrees-per-day
conversion. Reusing that machinery to recompute the index looked cheap.

## What was done

Monthly GLORYS surface velocities, 1993–2000, over 40–58°N and 70–45°W, which
spans the release line and the Scotian Shelf. Particles seeded along the Jutras
line, advected for three years on a two-day step, counting those that reached a
Labrador Shelf box and a Scotian Shelf box.

The script is [`lcr-extension-experiment.R`](lcr-extension-experiment.R).

## What happened

Particles reached the Labrador Shelf and **never reached the Scotian Shelf**. Not
some months, none at all. They also survived a median of 112 days of the 1095
asked for, leaving the domain to the east.

Tracing where they went explains it: they travelled south but **never west of
56°W**, and the Scotian Shelf begins around 59°W.

The eight-year mean surface flow says why. There is no coherent southwestward
shelf jet in it:

| Location | u (m/s) | v (m/s) | Direction |
|---|---|---|---|
| Release line (54°W, 53.6°N) | −0.015 | −0.077 | southwest |
| Labrador shelf (55°W, 52°N) | **+0.050** | −0.286 | south**east** |
| Newfoundland (52°W, 49°N) | +0.070 | −0.029 | southeast |
| Grand Banks (50°W, 46°N) | +0.008 | −0.003 | near stagnant |
| Scotian Shelf (62°W, 44°N) | +0.033 | −0.019 | east |

The Labrador Current is a narrow shelf-edge jet, of order 100 km wide. A monthly
mean smears it, and at the surface it competes with wind-driven Ekman drift
pointing offshore. What survives averaging is a southeastward drift into the open
Atlantic rather than a current following the shelf.

The Grand Banks row is the fatal one. That is where retroflection is decided, and
in the monthly mean the flow there is essentially zero. There is nothing to steer
a particle either way, so the quantity the index measures does not exist in this
field.

## Depth does not rescue it

The surface layer is the most Ekman-contaminated, so the same comparison was run
at 40 m:

| Location | surface u | 40 m u |
|---|---|---|
| Labrador shelf | +0.050 | +0.018 |
| Grand Banks | +0.008 | −0.012 |

Better, and still not a shelf-following jet. The eastward component on the
Labrador shelf shrinks but does not reverse.

## What would be needed

**Daily fields**, as Jutras et al. used. The eddies and the jet survive daily
averaging and do not survive monthly averaging.

That is a much larger undertaking than a bigger download. For 1993–2014 over this
domain, daily GLORYS is roughly 4 GB. Worse is the computation: weekly releases
over 22 years is about 1,100 releases, each tracking of order a thousand
particles for 1,095 daily steps, with four velocity evaluations per step. That is
on the order of a billion interpolations. In R, with `terra::extract` per step,
the pilot here ran about 11 ms per step for 33 particles, which extrapolates to
roughly a day of compute, and the daily fields cannot all be held in memory.

OceanParcels exists for exactly this, which is presumably why Jutras et al. used
it.

## Conclusion

Recomputing `LCR` from monthly Copernicus fields does not work, and the failure
is physical rather than a matter of tuning. The index depends on a narrow current
and a bifurcation that monthly averaging removes.

Anyone wanting to extend the series past 2014 should plan on daily fields and a
purpose-built Lagrangian framework, not on this package's advection.

Meanwhile `AMOC` from the RAPID array runs from 2004 to the present and describes
the basin-scale circulation the Labrador Current sits within. It is not a
substitute for `LCR`, but it is available for the years `LCR` is not.

## Reference

Jutras M, Dufour CO, Mucci A, Talbot LC (2023). Large-scale control of the
retroflection of the Labrador Current. *Nature Communications* **14**, 2623.
[doi:10.1038/s41467-023-38321-y](https://doi.org/10.1038/s41467-023-38321-y)
