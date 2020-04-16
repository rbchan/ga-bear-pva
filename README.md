# ga-bear-pva
Data and code for Central Georgia black bear population viability analysis (PVA). The PVA involved forecasting with a female-based population model fit to spatial capture-recapture data. No data on males is included in this repository.


## File descriptions

### Data files

#### [data/ga_bear_data_females.gzip](data/ga_bear_data_females.gzip)

- Compressed R data file with the encounter histories, trap coordinates, and trap operational status information.

#### [data/state-space360.tif](data/state-space360.tif)

- GeoTiff spatial layer defining the state space.


### R scripts

#### [R/sstats_maps.R](R/sstats_maps.R)

- Compute summary statistics and create maps of the detection locations, referenced by year, for each bear.


#### [R/fit_forecast.R](R/forecast_fn.R)

- Fit the model and forecast from it. Depends on the `forecasting` function in [R/forecast_fn.R](R/forecast_fn.R)


#### [R/forecast_fn.R](R/forecast_fn.R)

- Function for forecasting population dynamics and computing posterior predictive distributions.


#### [R/testS.R](R/testS.R)

- Assess sensitivity of results to specification of the state space.


### JAGS model

#### [R/dT-s0-rPDD2-phi0-pYB-sig0.jag](R/dT-s0-rPDD2-phi0-pYB-sig0.jag)

- The most supported model. 

