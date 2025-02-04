MERRA2: Data Prep
================
<br> Updated on April 17, 2024

**Overview**: This script cleans the MERRA2 data from 1980 through 1995.

Load packages:

``` r
library(dplyr)
library(lubridate)
```

Load data:

``` r
aod_raw = read.csv("../data/tavgM_2d_aer_Nx/TOTEXTTAUall_48x24.csv")
temp_strat_raw <-
  read.csv("../data/instM_3d_asm_Np/MERRA2_3dasm_temperature_50mb_48x24.csv")
```

Clean AOD and temperature data:

``` r
aod <- 
  aod_raw |> 
  rename("aod" = "TOTEXTTAU") |>
  mutate(date = as_date(date))

temp_strat <-
  temp_strat_raw |> 
  select(-lev) |> 
  rename("temp_strat" = `T`) |>
  mutate(date = as_date(date)) |>
  mutate(day = day(date))
```

Join the cleaned data and clean:

- reorder columns
- subset desired years for analysis
- arrange data as desired
- arrange columns as desired
- assign location weights

``` r
merra2 <-
  full_join(
    aod,
    temp_strat,
    by = c("lon", "lat", "date", "year", "month", "day")
  ) |>
  filter(year >= 1980 & year <= 1995) |>
  mutate(loc_weight = sqrt(cos(lat * pi / 180))) |>
  arrange(lon, lat, date) |>
  mutate(loc_id = 1:n(), .by = date) |>
  select(lon, lat, loc_id, loc_weight, date, year, month, day, everything())
```

Compute means and standard deviations for normalized anomalies for all
years (1980-1995):

``` r
merra2_means_and_sds <-
 merra2 |>
 summarise(
   aod_mean = mean(aod, na.rm = TRUE),
   temp_strat_mean = mean(temp_strat, na.rm = TRUE),
   aod_sd = sd(aod, na.rm = TRUE),
   temp_strat_sd = sd(temp_strat, na.rm = TRUE),
   .by = c(month, lon, lat, loc_weight, loc_id)
 )
```

Add means and center variables to prepare values to be used for training
ESN on all years (1980-1995):

``` r
merra2 <-
  left_join(merra2,
            merra2_means_and_sds,
            by = c("lon", "lat", "loc_id", "loc_weight", "month")) |>
  mutate(
    aod_anom = (aod - aod_mean) / (aod_sd),
    temp_strat_anom = (temp_strat - temp_strat_mean) / (temp_strat_sd)
  ) |>
  select(
    lon:day,
    aod,
    aod_mean,
    aod_sd,
    aod_anom,
    temp_strat,
    temp_strat_mean,
    temp_strat_sd,
    temp_strat_anom
  )
```

Save data:

``` r
write.csv(x = merra2, file = "../data/merra2.csv", row.names = FALSE)
```
