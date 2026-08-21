# derivoce: Derived Oceanographic Covariates for Species Distribution Models

Computes derived oceanographic covariates from gridded ocean data,
including spatial and temporal gradients, time-integrated variables,
temporal lags by calendar unit, and fluid dynamics such as eddy kinetic
energy and finite-time or finite-size Lyapunov exponents. Also derives
region-scale indices that describe an area rather than a cell: volume
transport across a section, water-mass fractions from
temperature-salinity endmember mixing, and regional property anomalies,
with named cases for Scotian Shelf inflow to the Gulf of Maine and the
Northeast Channel. Takes the output of datamatch::accessEnvDat() (an sf
point object per time step) and returns the same shape with derived
columns added, so derived covariates flow into a species distribution
model alongside the variables they came from.

## See also

Useful links:

- <https://github.com/chross22/derivoce>

- <https://camilleross.org/derivoce/>

- Report bugs at <https://github.com/chross22/derivoce/issues>

## Author

**Maintainer**: Camille Ross <camille.ross@maine.edu>
([ORCID](https://orcid.org/0000-0002-1428-2294))

Authors:

- Camille Ross <camille.ross@maine.edu>
  ([ORCID](https://orcid.org/0000-0002-1428-2294))
