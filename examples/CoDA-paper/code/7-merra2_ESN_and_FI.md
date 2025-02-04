MERRA2: ESN Feature Importance Analysis
================
<br> Updated on April 17, 2024

**Overview**: This script trains ESNs on MERRA-2 data and computes
feature importance.

# Set Up

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
library(furrr)
library(future)
library(ggplot2)
library(listenr)
library(lubridate)
library(purrr)
library(stringr)
library(tidyr)
```

# Data Preparation

Load data:

``` r
merra2 <-
  read.csv("../data/merra2.csv") |>
  mutate(date = as_date(date))
```

Generate AOD and stratospheric temperature matrices for fitting ESNs
using all years (1980-1995):

``` r
train_mat_aod <-
  merra2 |>
  select(date, lon, lat, aod_anom) |>
  pivot_wider(
    id_cols = c(lon, lat),
    names_from = date,
    values_from = aod_anom
  ) |>
  select(-lon,-lat) |>
  t()

train_mat_temp_strat <-
  merra2 |>
  select(date, lon, lat, temp_strat_anom) |>
  pivot_wider(
    id_cols = c(lon, lat),
    names_from = date,
    values_from = temp_strat_anom
  ) |>
  select(-lon,-lat) |>
  t()
```

Compute EOFs:

``` r
n_eofs = 5
eofs_aod = compute_eofs(Ztrain = train_mat_aod, n_eofs = n_eofs)
eofs_temp_strat = compute_eofs(Ztrain = train_mat_temp_strat, n_eofs = n_eofs)
phi_train = eofs_temp_strat$phi
```

Specify model inputs/outputs (AOD and stratospheric temperature are used
to forecast stratospheric temperature):

``` r
x = cbind(eofs_aod$train, eofs_temp_strat$train)
y = eofs_temp_strat$train
```

Save weights for later use:

``` r
weights = merra2 |> distinct(lon, lat, loc_weight) |> pull(loc_weight)
```

# Tuned Model and Feature Importance

Specify number of cores, number of reps (for FI), and block sizes:

``` r
ncores = 3
nreps = 10
blks = 1:6
```

Function for computing different feature importance scenarios:

``` r
compute_merra2_fi <- function(eesn, nblocks, cores, fp) {
    
    # Compute (weighted) ZFI on spatial scale
    zfi_spatial_wgt <-
      compute_ens_fi_fast(
        modelList = eesn,
        type = "zfi",
        var_groups = list(1:n_eofs, (n_eofs + 1):(2 * n_eofs)),
        blockSize = nblocks,
        y_spatial = train_mat_temp_strat,
        phi = phi_train,
        weights = weights,
        seed = 202304,
        cores = cores,
      ) |>
      map_df(.f = data.frame, .id = "ensemble")
    print("ZFI weighted")
    
    # Compute (weighted) PFI on spatial scale
    pfi_spatial_wgt <-
      compute_ens_fi_fast(
        modelList = eesn,
        type = "pfi",
        var_groups = list(1:n_eofs, (n_eofs + 1):(2 * n_eofs)),
        blockSize = nblocks,
        nreps = nreps,
        y_spatial = train_mat_temp_strat,
        phi = phi_train,
        weights = weights,
        seed = 202304,
        cores = cores
      ) |>
      map_df(.f = data.frame, .id = "ensemble") |>
      summarise(
        rmses_obs = mean(rmses_obs),
        rmses_adj = mean(rmses_adj),
        fi = mean(fi), 
        .by = c(ensemble, t_adj, t_forecasted, vars_adj)
      )
    print("PFI weighted")

    # Join FI results and clean up
    fi <-
      bind_rows(
        zfi_spatial_wgt |> mutate(feat_imp = "ZFI"),
        pfi_spatial_wgt |> mutate(feat_imp = "PFI")
      ) |>
      mutate(
        variable = ifelse(
          vars_adj == paste(1:n_eofs, collapse = ","),
          "aod",
          "temp_strat"
        ),
        t_forecasted = as_date(t_forecasted),
        t_adj = as_date(t_adj)
      ) |>
      mutate(nblocks = nblocks) |>
      select(
        ensemble,
        feat_imp,
        nblocks,
        variable,
        vars_adj, 
        t_adj, 
        t_forecasted,
        everything()
      )
    
    # Print block number completed
    print(nblocks)
  
  # Save results
  write.csv(fi, file = fp, row.names = FALSE)
  
}
```

``` r
if (!file.exists("../output/merra2_FI_eesn1_blk6.csv")) {
  eesn1 <-
    fit_Eesn(
      x = x,
      y = y,
      t = as.character(rownames(x)),
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
      internal_scaling = "joint",
      n_ensembles = 25,
      cores = ncores,
      seed = 202309
    )
  for (nblks in sort(blks, decreasing = T)) {
    compute_merra2_fi(
      eesn = eesn1,
      nblocks = nblks,
      cores = ncores,
      fp = paste0("../output/merra2_FI_eesn1_blk", nblks, ".csv")
    )
  }
}
```

``` r
if (!file.exists("../output/merra2_FI_eesn2_blk6.csv")) {
  eesn2 <-
    fit_Eesn(
      x = x,
      y = y,
      t = as.character(rownames(x)),
      tau = 1,
      m = 5,
      tau_emb = 1,
      U_width = 0.1,
      U_pi = 0.1,
      W_width = 0.1, 
      W_pi = 0.1,
      nu = 0.2,
      nh = 50,
      reg_par = 0.1,
      add_quad = FALSE,
      internal_scaling = "joint",
      n_ensembles = 25,
      cores = ncores,
      seed = 202309
    )
  for (nblks in sort(blks, decreasing = T)) {
    compute_merra2_fi(
      eesn = eesn2,
      nblocks = nblks,
      cores = ncores,
      fp = paste0("../output/merra2_FI_eesn2_blk", nblks, ".csv")
    )
  }
}
```

``` r
if (!file.exists("../output/merra2_FI_eesn3_blk6.csv")) {
  eesn3 <-
    fit_Eesn(
      x = x,
      y = y,
      t = as.character(rownames(x)),
      tau = 1,
      m = 5,
      tau_emb = 1,
      U_width = 0.1,
      U_pi = 0.1,
      W_width = 0.1, 
      W_pi = 0.1,
      nu = 0.2,
      nh = 2000,
      reg_par = 50,
      add_quad = FALSE,
      internal_scaling = "joint",
      n_ensembles = 25,
      cores = ncores,
      seed = 202309
    )
  for (nblks in sort(blks, decreasing = T)) {
    compute_merra2_fi(
      eesn = eesn3,
      nblocks = nblks,
      cores = ncores,
      fp = paste0("../output/merra2_FI_eesn3_blk", nblks, ".csv")
    )
  }
}
```
