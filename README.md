# ga-bear-pva
Data and code for Central Georgia black bear population viability analysis (PVA). The PVA involved forecasting with a female-based population model fit to spatial capture-recapture data.


## File descriptions

### Data files

#### [data/ga_bear_data_females.gzip](data/ga_bear_data_females.gzip)

- Compressed R data file with the encounter histories, trap coordinates, and trap operational status information.

#### [data/state-space360.tif](data/state-space360.tif)

- GeoTiff spatial layer defining the state space.


### R scripts

#### [R/sstats_maps.R](R/sstats_maps.R)

- Compute summary statistics and create maps of the detection locations, referenced by year, for each bear.


#### [R/fit.R](R/forecast_fn.R)

- Fit the model with the lowest WAIC score. The complete set of models was run on the [UGA cluster](https://gacrc.uga.edu/systems/). That code is not included in this repository. 


#### [R/forecast.R](R/forecast_fn.R)

- Forecast with the fitted model. Depends on the `forecast` function in [R/forecast_fn.R](R/forecast_fn.R)


#### [R/forecast_fn.R](R/forecast_fn.R)

- Function for forecasting population dynamics.


#### [R/testS.R](R/testS.R)

- Assess sensitivity of results to specification of the state space.


### JAGS model

#### [R/dT-s0-rPDD2-phi0-pYB-sig0.jag](R/dT-s0-rPDD2-phi0-pYB-sig0.jag)

- The most supported model. 


### [Figures](figs/)


## DOI

[![DOI](https://zenodo.org/badge/233102309.svg)](https://zenodo.org/badge/latestdoi/233102309)
