MERRA2: ESN Validation
================
<br> Updated on April 17, 2024

**Overview**: This script implements some validation techniques for the
MERRA-2 ESNs.

Install correct version of `listenr`:

``` r
# install.packages(
#   "../listenr-v0.3.4.tar.gz",
#   repos = NULL,
#   type = "source"
# )
```

Load packages:

``` r
library(dplyr)
library(ggplot2)
library(listenr)
library(lubridate)
library(purrr)
library(tidyr)
library(stringr)
```

Load data:

``` r
merra2 <-
  read.csv("../data/merra2.csv") |>
  mutate(date = as_date(date))
```

Functions for getting training/testing data predictions on spatial
scale:

``` r
get_train_preds <- function(esn) {
  data.frame(predict_esn(esn)$preds_ins %*% t(phi_train)) |>
    tibble::rownames_to_column(var = "date") |>
    pivot_longer(
      cols = -c(date),
      names_to = "loc_id",
      values_to = "temp_strat_anom_pred"
    ) |>
    mutate(
      loc_id = as.numeric(str_remove(loc_id, "X")),
      date = as_date(date),
      type = "Train"
    )
}

get_test_preds <- function(esn) {
  predict_esn(model = esn,
              x_new = x_test,
              t_new = t_test)$preds_new %*% t(phi_train) |>
    data.frame() |>
    tibble::rownames_to_column(var = "date") |>
    pivot_longer(
      cols = -c(date),
      names_to = "loc_id",
      values_to = "temp_strat_anom_pred"
    ) |>
    mutate(
      loc_id = as.numeric(str_remove(loc_id, "X")),
      date = as_date(str_remove(date, " \\+ 1")),
      type = "Test"
    )
}
```

Function for preparing train/test matrices:

``` r
create_mat <- function(train_test, var_name) {
  merra2_esn |>
    filter(data == {
      train_test
    }) |>
    select(date, lon, lat, {
      var_name
    }) |>
    pivot_wider(
      id_cols = c(lon, lat),
      names_from = date,
      values_from = {
        var_name
      }
    ) |>
    select(-lon,-lat) |>
    t()
}
```

Compute time series cross-validation RMSEs and predictions:

``` r
# Specify file paths
fp_preds = "../output/merra2_preds_strat.csv"
fp_rmses = "../output/merra2_rmses_strat.csv"

if (!file.exists(fp_preds) & !file.exists(fp_rmses)) {
  
  # Create lists for storing results
  preds_list <- list()
  rmses_list <- list()
  
  # Run process
  for(i in 1:13){
    
    # Split into training/testing
    merra2_train_test <-
      merra2 |>
      mutate(data = ifelse(year > 1981 + i, "test", "train")) |>
      select(data, date, month, loc_id, lon, lat, loc_weight, aod, temp_strat) |>
      arrange(lon, lat, date)
    
    # Compute monthly spatial means and sds
    merra2_means_and_sds <-
      merra2_train_test |>
      filter(data == "train") |>
      summarise(
        aod_mean = mean(aod, na.rm = TRUE),
        temp_strat_mean = mean(temp_strat, na.rm = TRUE),
        aod_sd = sd(aod, na.rm = TRUE),
        temp_strat_sd = sd(temp_strat, na.rm = TRUE),
        .by = c(month, lon, lat)
      )
    
    # Compute anomalies
    merra2_esn <-
      left_join(merra2_train_test,
                merra2_means_and_sds,
                by = c("month", "lon", "lat")) |>
      mutate(
        aod_anom = (aod - aod_mean) / (aod_sd),
        temp_strat_anom = (temp_strat - temp_strat_mean) / (temp_strat_sd)
      )
    
    # Prepare test/train matrices needed for ESN
    train_mat_aod <- create_mat('train', 'aod_anom')
    train_mat_temp_strat <- create_mat('train', 'temp_strat_anom')
    test_mat_aod <- create_mat('test', 'aod_anom')
    test_mat_temp_strat <- create_mat('test', 'temp_strat_anom')
    
    # Compute EOFs
    n_eofs = 5
    eofs_aod <-
      compute_eofs(Ztrain = train_mat_aod,
                   Ztest = test_mat_aod,
                   n_eofs = n_eofs)
    eofs_temp_strat <-
      compute_eofs(Ztrain = train_mat_temp_strat,
                   Ztest = test_mat_temp_strat,
                   n_eofs = n_eofs)
    
    # Prepare ESN inputs/outputs
    x_train = cbind(eofs_aod$train, eofs_temp_strat$train)
    y_train = eofs_temp_strat$train
    x_test = cbind(eofs_aod$test, eofs_temp_strat$test)
    
    # Extract times
    t_train = rownames(x_train)
    t_test = rownames(eofs_temp_strat$test)
    
    # Fit model
    eesn <-
      fit_Eesn(
        x = x_train,
        y = y_train,
        t = as.character(t_train),
        tau = 1,
        m = 5,
        tau_emb = 1,
        U_width = 0.1,
        U_pi = 0.1,
        W_width = 0.1,
        W_pi = 0.1,
        nu = 0.2,
        nh = 1000,
        reg_par = 5,
        add_quad = TRUE,
        n_ensembles = 25,
        internal_scaling = "joint", 
        cores = 6,
        seed = 20230
      )
    
    # Grab phi
    phi_train <- eofs_temp_strat$phi
    
    # Get train/test predictions
    preds_spatial_train <-
      map_df(.x = eesn, .f = get_train_preds, .id = "ensemble")
    preds_spatial_test <-
      map_df(.x = eesn, .f = get_test_preds, .id = "ensemble")
    
    # Join predictions
    preds_spatial <-
      rbind(preds_spatial_train, preds_spatial_test) |>
      left_join(
        merra2_esn |> select(
          loc_id,
          lon,
          lat,
          loc_weight,
          date,
          temp_strat,
          temp_strat_anom,
          temp_strat_mean,
          temp_strat_sd
        ),
        by = c("date", "loc_id")
      ) |>
      mutate(temp_strat_pred =
               (temp_strat_anom_pred * temp_strat_sd) + temp_strat_mean) |>
      select(
        ensemble,
        date,
        lat,
        lon,
        loc_weight,
        temp_strat,
        temp_strat_pred,
        temp_strat_anom,
        temp_strat_anom_pred,
        type
      ) |>
      mutate(year = 1981 + i)
    
    # Compute RMSE
    rmses_spatial <-
      preds_spatial |>
      summarize(
        rmse = sqrt(mean((temp_strat - temp_strat_anom) ^ 2)),
        rmse_wgt = sqrt(sum((
          loc_weight * (temp_strat - temp_strat_pred) ^ 2
        )) / sum(loc_weight)),
        rmse_anom = sqrt(mean((temp_strat_anom - temp_strat_anom_pred) ^ 2)),
        rmse_anom_wgt = sqrt(sum((
          loc_weight * (temp_strat_anom - temp_strat_anom_pred) ^ 2
        )) / sum(loc_weight)),
        .by = c(date, type, ensemble)
      ) |>
      mutate(year = 1981 + i)
    
    # Add results to lists
    preds_list[[i]] <- preds_spatial |> filter(ensemble == 1)
    rmses_list[[i]] <- rmses_spatial
    
  }
 
  # Join all results
  pred_df <- do.call("rbind", preds_list)
  rmse_df <- do.call("rbind", rmses_list)

  # Save results    
  write.csv(pred_df, file = fp_preds, row.names = FALSE)
  write.csv(rmse_df, file = fp_rmses, row.names = FALSE)

} else {
  
  # Load results
  pred_df = read.csv(fp_preds) |> mutate(date = as_date(date))
  rmse_df = read.csv(fp_rmses) |> mutate(date = as_date(date))

}
```

Time series cross-validation results (raw temperatures):

``` r
# Specify figure file path
fp_plot_merra_blk_rmses = "../figs/merra2_rmse.png"

# Create/load plot
if (!file.exists(fp_plot_merra_blk_rmses)) {
  
  plot_merra_blk_rmses <-
    rmse_df |>
    ggplot(aes(x = date, y = rmse_wgt, color = type)) +
    geom_line(aes(group = ensemble), alpha = 0.75, linewidth = 0.25) +
    scale_color_manual(values = c(wesanderson::wes_palettes$Zissou1[2], "black")) +
    facet_grid(year ~ .) +
    labs(
      y = "RMSE",
      x = "Date",
      title = "ESNs trained through year shown on row label",
      color = ""
    ) +
    geom_vline(xintercept = as_date("1991-06-15"), linetype = 2) +
    theme_bw(base_size = 20) +
    theme(title = element_text(size = 18)) +
    guides(
      colour = guide_legend(override.aes = list(alpha = 1, linewidth = 0.5))
    )
  
   ggsave(
     filename = fp_plot_merra_blk_rmses,
     plot = plot_merra_blk_rmses,
     dpi = 500,
     height = 17,
     width = 15
    )
   
}

# Load plot
knitr::include_graphics(fp_plot_merra_blk_rmses)
```

<img src="../figs/merra2_rmse.png" width="7500" />

Time series cross-validation results (anomalies):

``` r
# Specify figure file path
fp_plot_merra_blk_rmses_anom = "../figs/merra2_rmse_anom.png"

# Create/load plot
if (!file.exists(fp_plot_merra_blk_rmses_anom)) {
  
  plot_merra_blk_rmses_anom <-
    rmse_df |>
    ggplot(aes(x = date, y = rmse_anom_wgt, color = type)) +
    geom_line(aes(group = ensemble), alpha = 0.25, linewidth = 0.75) +
    scale_color_manual(values = c(wesanderson::wes_palettes$Zissou1[2], "black")) +
    facet_grid(year ~ .) +
    labs(
      y = "RMSE",
      x = "Date",
      title = "ESNs trained through year shown on row label",
      color = ""
    ) +
    geom_vline(xintercept = as_date("1991-06-15"), linetype = 2) +
    scale_y_continuous(breaks = seq(0, 18, 5)) + 
    theme_bw(base_size = 20) +
    theme(title = element_text(size = 18)) + 
    guides(
      colour = guide_legend(override.aes = list(alpha = 1, linewidth = 0.5))
    )
  
   ggsave(
     filename = fp_plot_merra_blk_rmses_anom,
     plot = plot_merra_blk_rmses_anom,
     dpi = 500,
     height = 17,
     width = 15
    )
   
}


# Load plot
knitr::include_graphics(fp_plot_merra_blk_rmses_anom)
```

<img src="../figs/merra2_rmse_anom.png" width="7500" />

Time series cross-validation results (anomalies, subset of years):

``` r
# Specify figure file path
fp_plot_merra_blk_rmses_anom_even_yrs = "../figs/merra2_rmse_anom_even_yrs.png"

# Create/load plot
if (!file.exists(fp_plot_merra_blk_rmses_anom_even_yrs)) {
  
  plot_merra_blk_rmses_anom_even_yrs <-
    rmse_df |>
    filter((year %% 2) == 0) |>
    ggplot(aes(x = date, y = rmse_anom_wgt, color = type)) +
    geom_line(aes(group = ensemble), alpha = 0.5, linewidth = 0.5) +
    scale_color_manual(values = c(wesanderson::wes_palettes$Zissou1[2], "black")) +
    facet_grid(year ~ .) +
    labs(
      y = "RMSE",
      x = "Date",
      title = "ESNs trained through year shown on row label",
      color = ""
    ) +
    geom_vline(xintercept = as_date("1991-06-15"), linetype = 2) +
    theme_bw(base_size = 22) +
    theme(title = element_text(size = 18)) +
    scale_y_continuous(breaks = seq(0, 18, 5)) +
    guides(
      colour = guide_legend(override.aes = list(alpha = 1, linewidth = 0.5))
    )
  
   ggsave(
     filename = fp_plot_merra_blk_rmses_anom_even_yrs,
     plot = plot_merra_blk_rmses_anom_even_yrs,
     dpi = 500,
     height = 12,
     width = 16
    )
   
}

# Load plot
knitr::include_graphics(fp_plot_merra_blk_rmses_anom_even_yrs)
```

<img src="../figs/merra2_rmse_anom_even_yrs.png" width="8000" />

Predictions versus observations:

``` r
# Specify figure file path
fp_plot_merra_blk_pred_vs_obs = "../figs/merra2_obs_pred.png"

# Create/load plot
if (!file.exists(fp_plot_merra_blk_pred_vs_obs)) {
  
  plot_merra_blk_pred_vs_obs <-
    pred_df |>
    mutate(type = factor(type, levels = c("Test", "Train"))) |>
    arrange(type) |>
    ggplot(aes(x = temp_strat, y = temp_strat_pred, color = type)) +
    geom_point(alpha = 0.1) +
    geom_abline(intercept = 0,
                slope = 1,
                color = "red") +
    scale_color_manual(
      values = c(wesanderson::wes_palettes$Zissou1[2], 'black')
    ) +
    labs(
      x = "Observed Stratospheric Temperature",
      y = 'Predicted Stratospheric Temperature',
      color = "",
      title = "ESN trained through year shown on label"
    ) +
    facet_wrap(year ~ ., ncol = 4) +
    guides(color = guide_legend(override.aes = list(alpha = 1))) +
    theme_bw(base_size = 16) +
    ylim(150, 275) +
    theme(aspect.ratio = 1)
  
  ggsave(
    filename = fp_plot_merra_blk_pred_vs_obs,
    plot = plot_merra_blk_pred_vs_obs,
    dpi = 500,
    height = 12,
    width = 12
  )
  
}

# Load plot
knitr::include_graphics(fp_plot_merra_blk_pred_vs_obs)
```

<img src="../figs/merra2_obs_pred.png" width="6000" />

Residuals versus observations (not included in paper)

``` r
# Specify figure file path
fp_plot_merra_res_vs_obs = "../figs/merra2_res_obs.png"

# Create/load plot
if (!file.exists(fp_plot_merra_res_vs_obs)) {
  
  plot_merra_res_vs_obs <-
    pred_df |>
    filter(ensemble == 1) |>
    mutate(type = factor(type, levels = c("Test", "Train"))) |>
    arrange(type) |>
    ggplot(aes(
      x = temp_strat,
      y = (temp_strat_pred - temp_strat),
      color = type
    )) +
    geom_point(alpha = 0.1) +
    geom_hline(yintercept = 0, color = "red") +
    scale_color_manual(values = c(wesanderson::wes_palettes$Zissou1[2], 'black')) +
    labs(
      x = "Observed Stratospheric Temperature",
      y = 'Predicted Stratospheric Temperature',
      color = "",
      title = "ESN trained through year shown on label"
    ) +
    facet_wrap(year ~ ., ncol = 4) +
    guides(color = guide_legend(override.aes = list(alpha = 1))) +
    theme_bw(base_size = 16) +
    theme(aspect.ratio = 0.75)
  
  ggsave(
    filename = fp_plot_merra_res_vs_obs,
    plot = plot_merra_res_vs_obs,
    dpi = 500,
    height = 12,
    width = 12
  )
  
}

# Load plot
knitr::include_graphics(fp_plot_merra_res_vs_obs)
```

<img src="../figs/merra2_res_obs.png" width="6000" />
