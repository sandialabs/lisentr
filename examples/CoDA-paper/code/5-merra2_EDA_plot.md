MERRA2: EDA
================
<br> Updated on April 17, 2024

**Overview**: This script creates visualizations of MERRA2 data.

# Set Up

Load packages:

``` r
library(dplyr)
library(ggplot2)
library(lubridate)
library(patchwork)
library(tidyr)
```

Load data:

``` r
merra2 <-
  read.csv("../data/merra2.csv") |>
  mutate(date = as_date(date))
```

# Line Plot

Prepare data for plots:

``` r
# Prepare observed data for plotting
df_obs <- 
  merra2 |>
  summarize(
    sum_weight = sum(loc_weight),
    AOD = sum(loc_weight * aod) / sum_weight,
    Stratospheric_Temperature = sum(loc_weight * temp_strat) / sum_weight,
    .by = date
  ) |>
  pivot_longer(
    cols = AOD:Stratospheric_Temperature,
    names_to = "variable",
    values_to = "value"
  ) |>
  mutate(data = "Observed")

# Prepare anomalies for plotting
df_anom <- 
  merra2 |>
  summarize(
    sum_weight = sum(loc_weight),
    AOD = sum(loc_weight * aod_anom) / sum_weight,
    Stratospheric_Temperature = sum(loc_weight * temp_strat_anom) / sum_weight,
    .by = date
  ) |>
  pivot_longer(
    cols = AOD:Stratospheric_Temperature,
    names_to = "variable",
    values_to = "value"
  ) |>
  mutate(data = "Normalized Anomalies")

# Join data
df_wgt_aves <-
  rbind(df_obs, df_anom) |>
  mutate(variable = stringr::str_replace(variable, "_", " ")) |>
  mutate(facet = paste0(variable, ": ", data)) |>
  mutate(
    facet = forcats::fct_relevel(
      facet,
      "AOD: Observed",
      "AOD: Normalized Anomalies",
      "Stratospheric Temperature: Observed",
      "Stratospheric Temperature: Normalized Anomalies"
    )
  )
```

``` r
# Specify figure file path
fp_p_line = "../figs/merra2_weighted_global_averages.png"

# Create/load plot
if (!file.exists(fp_p_line)) {
  
  # Create plot
  p_line_plot <-
    ggplot(
      data = df_wgt_aves, 
      aes(
        x = date,
        y = value,
        group = interaction(variable, data)
      )
    ) +
    geom_vline(
      xintercept = as_date("1991-06-12"),
      linetype = 2,
      linewidth = 0.75
    ) +
    geom_vline(
      xintercept = as_date("1982-03-29"),
      linetype = 2,
      color = "gray",
      linewidth = 0.75
    ) +
    geom_line() +
    facet_wrap(facet ~ ., scales = "free") +
    theme_bw(base_size = 18) +
    labs(x = "Date", y = "") +
    theme(strip.background = element_blank(), axis.title.y = element_blank())
  
  # Save plot
  ggsave(
    filename = fp_p_line,
    plot = p_line_plot,
    dpi = 500,
    width = 14,
    height = 7
  )
  
}

# Load plot
knitr::include_graphics(fp_p_line)
```

<img src="../figs/merra2_weighted_global_averages.png" width="7000" />

# Heatmap

``` r
# Specify figure file path
fp_p_spatial = "../figs/merra2_heatmaps_1991.png"

# Create/load plot
if (!file.exists(fp_p_spatial)) {
  
  # Access world map data:
  world_map = map_data("world")

  # AOD heatmap
  p_aod <-
    ggplot(
      data = merra2 |> filter(year == 1991),
      aes(x = lon, y = lat, fill = aod_anom)
    ) +
    geom_tile() +
    geom_map(
      data = world_map,
      map = world_map,
      aes(x = long, y = lat, map_id = region),
      fill = NA,
      color = "black",
      linewidth = 0.5
    ) +
    scale_fill_distiller(palette = "RdBu") +
    facet_wrap(~ date, ncol = 3) +
    labs(
      x = "Longitude",
      y = "Latitude",
      title = "AOD Monthly Normalized Anomalies",
      fill = "Normalized Anomaly"
    ) +
    theme_bw(base_size = 16) +
    theme(
      legend.key = element_blank(),
      strip.background = element_blank(),
      legend.position = "bottom"
    )
  
  # Strat temp heatmap
  p_temp_strat <-
    ggplot(
      data = merra2 |> filter(year == 1991),
      mapping = aes(x = lon, y = lat, fill = temp_strat_anom)
    ) +
    geom_tile() +
    geom_map(
      data = world_map,
      map = world_map,
      aes(x = long, y = lat, map_id = region),
      fill = NA,
      color = "black"
    ) +
    facet_wrap( ~ date, ncol = 3) +
    scale_fill_distiller(palette = "RdBu") +
    labs(
      x = "Longitude",
      y = "Latitude",
      title = "Stratospheric Temperature Monthly Normalized Anomalies",
      fill = "Normalized Anomaly"
    ) +
    theme_bw(base_size = 16) +
    theme(
      legend.key = element_blank(),
      strip.background = element_blank(),
      legend.position = "bottom"
    )
  
  # Join plots
  p_spatial <- p_aod + p_temp_strat

  # Save plot
  ggsave(
    filename = fp_p_spatial,
    plot = p_spatial,
    dpi = 500,
    width = 20,
    height = 11
  )
  
}

# Load plot
knitr::include_graphics(fp_p_spatial)
```

<img src="../figs/merra2_heatmaps_1991.png" width="10000" />
