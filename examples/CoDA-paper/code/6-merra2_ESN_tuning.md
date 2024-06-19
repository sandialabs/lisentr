MERRA2: Hyperparameter Tuning
================
<br> Updated on April 17, 2024

**Overview**: This script implements a hyperparameter search and
compares RMSE variability across different ensembles sizes for ESNs
trained using data from 1980 through 1995.

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
library(latex2exp)
library(listenr)
library(lubridate)
library(purrr)
library(stringr)
library(tidyr)
```

Load data:

``` r
merra2 <-
  read.csv("../data/merra2.csv") |>
  mutate(date = as_date(date))
```

Extract weights for later use:

``` r
weights = merra2 |> distinct(lon, lat, loc_weight) |> pull(loc_weight)
```

# ESN Hyper Parameter Tuning

Assign training/testing data labels:

``` r
merra2_hs_train_test <-
  merra2 |>
  filter(year >= 1980 & year <= 1995) |>
  mutate(data = ifelse(year > 1993, "test", "train")) |>
  select(
    data,
    everything(),
    -aod_mean,
    -aod_sd,
    -aod_anom,
    -temp_strat_mean,
    -temp_strat_sd,
    -temp_strat_anom
  )
```

Compute spatial means on training data:

``` r
merra2_hs_anom_spatial_means_and_sds <-
  merra2_hs_train_test |>
  filter(data == "train") |>
  summarise(
    aod_mean = mean(aod, na.rm = TRUE),
    temp_strat_mean = mean(temp_strat, na.rm = TRUE),
    aod_sd = sd(aod, na.rm = TRUE),
    temp_strat_sd = sd(temp_strat, na.rm = TRUE),
    .by = c(month, lon, lat)
  )
```

Add means and center variables to prepare values to be used in the ESN:

``` r
merra2_hs <-
  left_join(
    merra2_hs_train_test,
    merra2_hs_anom_spatial_means_and_sds,
    by = c("month", "lon", "lat")
  ) |>
  mutate(
    aod_anom = (aod - aod_mean) / (aod_sd),
    temp_strat_anom = (temp_strat - temp_strat_mean) / (temp_strat_sd)
  )
```

Prepare training data for ESN (only AOD and stratospheric temperature):

``` r
merra2_hs_train = merra2_hs |> filter(data == "train")

train_mat_hs_aod <-
  merra2_hs_train |>
  select(date, lon, lat, aod_anom) |>
  pivot_wider(
    id_cols = c(lon, lat),
    names_from = date,
    values_from = aod_anom
  ) |>
  select(-lon,-lat) |>
  t()

train_mat_hs_temp_strat <-
  merra2_hs_train |>
  select(date, lon, lat, temp_strat_anom) |>
  pivot_wider(
    id_cols = c(lon, lat),
    names_from = date,
    values_from = temp_strat_anom
  ) |>
  select(-lon,-lat) |>
  t()
```

Prepare testing data for ESN:

``` r
merra2_hs_test = merra2_hs |> filter(data == "test")

test_mat_hs_aod <-
  merra2_hs_test |>
  select(date, lon, lat, aod_anom) |>
  pivot_wider(
    id_cols = c(lon, lat),
    names_from = date,
    values_from = aod_anom
  ) |>
  select(-lon,-lat) |>
  t()

test_mat_hs_temp_strat <-
  merra2_hs_test |>
  select(date, lon, lat, temp_strat_anom) |>
  pivot_wider(
    id_cols = c(lon, lat),
    names_from = date,
    values_from = temp_strat_anom
  ) |>
  select(-lon,-lat) |>
  t()
```

Specify values for hyper-parameter search:

``` r
tuning_params <- 
  list(
    nensm = 5,
    cores = 7,
    tau = 1,
    m = 5,
    tau_emb = 1,
    nh = c(50, 100, 500, 1000, 2000),
    U_width = c(0.1, 0.5, 1),
    W_width = c(0.1, 0.5, 1),
    U_pi = c(0.1, 0.5, 0.9),
    W_pi = c(0.1, 0.5, 0.9),
    nu = c(0, 0.2, 0.4, 0.6, 0.8, 1),
    reg_par = c(0.1, 5, 20, 50, 100, 500),
    add_quad = c(FALSE, TRUE),
    internal_scaling = "joint",
    seed = 20231002
  )
```

## 5 PCs

Compute EOFs:

``` r
eofs_hs_aod_pcs5 <-
  compute_eofs(
    Ztrain = train_mat_hs_aod,
    Ztest = test_mat_hs_aod,
    n_eofs = 5
  )
eofs_hs_temp_strat_pcs5 <-
  compute_eofs(
    Ztrain = train_mat_hs_temp_strat,
    Ztest = test_mat_hs_temp_strat,
    n_eofs = 5
  )
```

Determine percent of variability explained by first 5 PCs:

``` r
aod_first_five = sum((eofs_hs_aod_pcs5$svd_vals$d[1:5])^2)
aod_total = sum((eofs_hs_aod_pcs5$svd_vals$d)^2)
print(paste("AOD:", round((aod_first_five / aod_total) * 100, 2)))
```

    ## [1] "AOD: 85.63"

``` r
temp_first_five = sum((eofs_hs_temp_strat_pcs5$svd_vals$d[1:5])^2)
temp_total = sum((eofs_hs_temp_strat_pcs5$svd_vals$d)^2)
print(paste("Strat temp:", round((temp_first_five / temp_total) * 100, 2)))
```

    ## [1] "Strat temp: 59.86"

Specify training/testing inputs/outputs (AOD and stratospheric
temperature are used to forecast stratospheric temperature):

``` r
x_hs_train_pcs5 = cbind(eofs_hs_aod_pcs5$train, eofs_hs_temp_strat_pcs5$train)
y_hs_train_pcs5 = eofs_hs_temp_strat_pcs5$train
x_hs_test_pcs5 = cbind(eofs_hs_aod_pcs5$test, eofs_hs_temp_strat_pcs5$test)
y_hs_test_pcs5 = eofs_hs_temp_strat_pcs5$test
```

Extract the times associated with the training/testing data (note that
the times for `x` and `y` are the same):

``` r
t_hs_train_pcs5 = rownames(x_hs_train_pcs5)
t_hs_test_pcs5 = rownames(x_hs_test_pcs5)
```

Hyperparameter search (weighted RMSEs):

``` r
fp_hs_pcs5w = "../output/merra2_hp_search_pcs05_weighted.csv"
if (!file.exists(fp_hs_pcs5w)) {
  start <- Sys.time()
  hs_pcs5w <- hyperparameter_search(
    x = x_hs_train_pcs5,
    y = y_hs_train_pcs5,
    t = as.character(t_hs_train_pcs5),
    x_test = x_hs_test_pcs5,
    t_test = t_hs_test_pcs5,
    phi = eofs_hs_temp_strat_pcs5$phi,
    obs_train = train_mat_hs_temp_strat,
    obs_test = test_mat_hs_temp_strat,
    tau = tuning_params$tau,
    m = tuning_params$m,
    tau_emb = tuning_params$tau_emb,
    nh = tuning_params$nh,
    U_width = tuning_params$U_width,
    W_width = tuning_params$W_width,
    U_pi = tuning_params$U_pi,
    W_pi = tuning_params$W_pi,
    nu = tuning_params$nu,
    reg_par = tuning_params$reg_par,
    add_quad = tuning_params$add_quad,
    internal_scaling = tuning_params$internal_scaling,
    seed = tuning_params$seed,
    n_ensembles = tuning_params$nensm,
    cores = tuning_params$cores,
    weights = weights
  )
  finish <- Sys.time()
  finish - start
  write.csv(hs_pcs5w$output, file = fp_hs_pcs5w, row.names = FALSE)
}
hs_pcs5w = read.csv(fp_hs_pcs5w)
```

## 10 PCs

Compute EOFs:

``` r
eofs_hs_aod_pcs10 <-
  compute_eofs(
    Ztrain = train_mat_hs_aod,
    Ztest = test_mat_hs_aod,
    n_eofs = 10
  )
eofs_hs_temp_strat_pcs10 <- 
  compute_eofs(
    Ztrain = train_mat_hs_temp_strat,
    Ztest = test_mat_hs_temp_strat,
    n_eofs = 10
  )
```

Determine percent of variability explained by first 10 PCs:

``` r
aod_first_five = sum((eofs_hs_aod_pcs5$svd_vals$d[1:10])^2)
aod_total = sum((eofs_hs_aod_pcs5$svd_vals$d)^2)
print(paste("AOD:", round((aod_first_five / aod_total) * 100, 2)))
```

    ## [1] "AOD: 88.78"

``` r
temp_first_five = sum((eofs_hs_temp_strat_pcs5$svd_vals$d[1:10])^2)
temp_total = sum((eofs_hs_temp_strat_pcs5$svd_vals$d)^2)
print(paste("Strat temp:", round((temp_first_five / temp_total) * 100, 2)))
```

    ## [1] "Strat temp: 76.58"

Specify training/testing inputs/outputs (AOD and stratospheric
temperature are used to forecast stratospheric temperature):

``` r
x_hs_train_pcs10 = cbind(eofs_hs_aod_pcs10$train, eofs_hs_temp_strat_pcs10$train)
y_hs_train_pcs10 = eofs_hs_temp_strat_pcs10$train
x_hs_test_pcs10 = cbind(eofs_hs_aod_pcs10$test, eofs_hs_temp_strat_pcs10$test)
y_hs_test_pcs10 = eofs_hs_temp_strat_pcs10$test
```

Extract the times associated with the training/testing data (note that
the times for `x` and `y` are the same):

``` r
t_hs_train_pcs10 = rownames(x_hs_train_pcs10)
t_hs_test_pcs10 = rownames(x_hs_test_pcs10)
```

Hyperparameter search (weighted RMSEs):

``` r
fp_hs_pcs10w = "../output/merra2_hp_search_pcs10_weighted.csv"
if (!file.exists(fp_hs_pcs10w)) {
  start <- Sys.time()
  hs_pcs10w <- hyperparameter_search(
    x = x_hs_train_pcs10,
    y = y_hs_train_pcs10,
    t = as.character(t_hs_train_pcs10),
    x_test = x_hs_test_pcs10,
    t_test = t_hs_test_pcs10,
    phi = eofs_hs_temp_strat_pcs10$phi,
    obs_train = train_mat_hs_temp_strat,
    obs_test = test_mat_hs_temp_strat,
    tau = tuning_params$tau,
    m = tuning_params$m,
    tau_emb = tuning_params$tau_emb,
    nh = tuning_params$nh,
    U_width = tuning_params$U_width,
    W_width = tuning_params$W_width,
    U_pi = tuning_params$U_pi,
    W_pi = tuning_params$W_pi,
    nu = tuning_params$nu,
    reg_par = tuning_params$reg_par,
    add_quad = tuning_params$add_quad,
    internal_scaling = tuning_params$internal_scaling,
    seed = tuning_params$seed,
    n_ensembles = tuning_params$nensm,
    cores = tuning_params$cores,
    weights = weights
  )
  finish <- Sys.time()
  finish - start
  write.csv(hs_pcs10w$output, file = fp_hs_pcs10w, row.names = FALSE)
}
hs_pcs10w = read.csv(fp_hs_pcs10w)
```

## Results

Make list of params:

``` r
params <-
  c("npcs",
    "nh",
    "U_width",
    "W_width",
    "U_pi",
    "W_pi",
    "nu",
    "reg_par",
    "add_quad")
```

Join RMSEs:

``` r
hs <-
  bind_rows(hs_pcs5w |> mutate(npcs = 5),
            hs_pcs10w |> mutate(npcs = 10)) |>
  select(tau, m, tau_emb, npcs, everything()) |>
  mutate_at(
    .vars = vars("tau", "m", "tau_emb", params),
    .funs = factor
  )
```

Print hyper-parameter values with best test RMSEs:

``` r
hs |>
  select(-c(tau, m, tau_emb)) |>
  select(npcs, everything()) |>
  arrange(mRMSE_te) |>
  slice(1:10)
```

    ##    npcs   nh U_width W_width U_pi W_pi nu reg_par add_quad  mRMSE_tr
    ## 1     5 1000     0.1     0.1  0.1  0.1  0       5     TRUE 0.7180083
    ## 2     5 1000     0.1     0.5  0.1  0.1  0       5     TRUE 0.7180083
    ## 3     5 1000     0.1       1  0.1  0.1  0       5     TRUE 0.7180083
    ## 4     5 1000     0.1     0.1  0.1  0.5  0       5     TRUE 0.7180083
    ## 5     5 1000     0.1     0.5  0.1  0.5  0       5     TRUE 0.7180083
    ## 6     5 1000     0.1       1  0.1  0.5  0       5     TRUE 0.7180083
    ## 7     5 1000     0.1     0.1  0.1  0.9  0       5     TRUE 0.7180083
    ## 8     5 1000     0.1     0.5  0.1  0.9  0       5     TRUE 0.7180083
    ## 9     5 1000     0.1       1  0.1  0.9  0       5     TRUE 0.7180083
    ## 10    5 1000     0.1     0.1  0.1  0.1  0       5    FALSE 0.7295083
    ##     mRMSE_te
    ## 1  0.8228017
    ## 2  0.8228017
    ## 3  0.8228017
    ## 4  0.8228017
    ## 5  0.8228017
    ## 6  0.8228017
    ## 7  0.8228017
    ## 8  0.8228017
    ## 9  0.8228017
    ## 10 0.8236659

Print hyper-parameter values with best test RMSEs for LaTeX:

``` r
hs |>
  select(-c(tau, m, tau_emb, mRMSE_tr)) |>
  select(npcs, everything()) |>
  arrange(mRMSE_te) |>
  slice(1:25) |>
  xtable::xtable(digits = 6)
```

    ## % latex table generated in R 4.3.1 by xtable 1.8-4 package
    ## % Wed Apr 17 11:26:18 2024
    ## \begin{table}[ht]
    ## \centering
    ## \begin{tabular}{rlllllllllr}
    ##   \hline
    ##  & npcs & nh & U\_width & W\_width & U\_pi & W\_pi & nu & reg\_par & add\_quad & mRMSE\_te \\ 
    ##   \hline
    ## 1 & 5 & 1000 & 0.1 & 0.1 & 0.1 & 0.1 & 0 & 5 & TRUE & 0.822802 \\ 
    ##   2 & 5 & 1000 & 0.1 & 0.5 & 0.1 & 0.1 & 0 & 5 & TRUE & 0.822802 \\ 
    ##   3 & 5 & 1000 & 0.1 & 1 & 0.1 & 0.1 & 0 & 5 & TRUE & 0.822802 \\ 
    ##   4 & 5 & 1000 & 0.1 & 0.1 & 0.1 & 0.5 & 0 & 5 & TRUE & 0.822802 \\ 
    ##   5 & 5 & 1000 & 0.1 & 0.5 & 0.1 & 0.5 & 0 & 5 & TRUE & 0.822802 \\ 
    ##   6 & 5 & 1000 & 0.1 & 1 & 0.1 & 0.5 & 0 & 5 & TRUE & 0.822802 \\ 
    ##   7 & 5 & 1000 & 0.1 & 0.1 & 0.1 & 0.9 & 0 & 5 & TRUE & 0.822802 \\ 
    ##   8 & 5 & 1000 & 0.1 & 0.5 & 0.1 & 0.9 & 0 & 5 & TRUE & 0.822802 \\ 
    ##   9 & 5 & 1000 & 0.1 & 1 & 0.1 & 0.9 & 0 & 5 & TRUE & 0.822802 \\ 
    ##   10 & 5 & 1000 & 0.1 & 0.1 & 0.1 & 0.1 & 0 & 5 & FALSE & 0.823666 \\ 
    ##   11 & 5 & 1000 & 0.1 & 0.5 & 0.1 & 0.1 & 0 & 5 & FALSE & 0.823666 \\ 
    ##   12 & 5 & 1000 & 0.1 & 1 & 0.1 & 0.1 & 0 & 5 & FALSE & 0.823666 \\ 
    ##   13 & 5 & 1000 & 0.1 & 0.1 & 0.1 & 0.5 & 0 & 5 & FALSE & 0.823666 \\ 
    ##   14 & 5 & 1000 & 0.1 & 0.5 & 0.1 & 0.5 & 0 & 5 & FALSE & 0.823666 \\ 
    ##   15 & 5 & 1000 & 0.1 & 1 & 0.1 & 0.5 & 0 & 5 & FALSE & 0.823666 \\ 
    ##   16 & 5 & 1000 & 0.1 & 0.1 & 0.1 & 0.9 & 0 & 5 & FALSE & 0.823666 \\ 
    ##   17 & 5 & 1000 & 0.1 & 0.5 & 0.1 & 0.9 & 0 & 5 & FALSE & 0.823666 \\ 
    ##   18 & 5 & 1000 & 0.1 & 1 & 0.1 & 0.9 & 0 & 5 & FALSE & 0.823666 \\ 
    ##   19 & 5 & 1000 & 0.1 & 0.1 & 0.1 & 0.9 & 0.2 & 5 & TRUE & 0.824477 \\ 
    ##   20 & 5 & 1000 & 0.1 & 0.5 & 0.1 & 0.9 & 0.2 & 5 & TRUE & 0.824477 \\ 
    ##   21 & 5 & 1000 & 0.1 & 1 & 0.1 & 0.9 & 0.2 & 5 & TRUE & 0.824477 \\ 
    ##   22 & 5 & 1000 & 0.1 & 0.1 & 0.1 & 0.5 & 0.2 & 5 & TRUE & 0.824612 \\ 
    ##   23 & 5 & 1000 & 0.1 & 0.5 & 0.1 & 0.5 & 0.2 & 5 & TRUE & 0.824612 \\ 
    ##   24 & 5 & 1000 & 0.1 & 1 & 0.1 & 0.5 & 0.2 & 5 & TRUE & 0.824612 \\ 
    ##   25 & 5 & 1000 & 0.1 & 0.1 & 0.1 & 0.1 & 0.2 & 5 & TRUE & 0.824966 \\ 
    ##    \hline
    ## \end{tabular}
    ## \end{table}

## Supplemental Visualizations

Create a list of latex versions of the parameters:

``` r
params_tex <-
  c(
    "No. PCs",
    TeX("$n_h$"),
    TeX("$a_u$"),
    TeX("$a_w$"),
    TeX("$\\pi_u$"),
    TeX("$\\pi_w$"),
    TeX("$\\nu$"),
    TeX("$\\lambda_r$"),
    "Quadratic?"
  )
```

These visualizations lead me to select 5 PCs over 10 PCs and to
determine that unimportant hyper parameters are W_width, U_pi, and W_pi:

<img src="../figs/merra2_hyper_param_all.png" width="2400" />

These visualizations tell me that the parameters most associated with a
low testing RMSE are nh and U_width. These are followed by nu and
reg_par. It is still unclear whether add_quad is helpful or not.

<img src="../figs/merra2_hyper_param_sub.png" width="2400" />

U_with = 0.1 performs the best across all values of nh and which value
of reg_par does the best is dependent on U_width and nh:

<img src="../figs/merra2_hyper_param_uwidth.png" width="3200" />

Best reg_par does not appear to be dependent on nu and in general, a
lower nu performs better within each value of nh (going forward subset
nu = 0.2):

<img src="../figs/merra2_hyper_param_nu.png" width="3200" />

These plots highlight the best performing model (nh of 1000, reg_par of
5), but there are other cases that perform close including nh of 2000
and reg_par of 50. Even nh of 50 with reg_par of 0.1 is within a range
that might be considered not practically different. Note that this plot
also suggests that add_quad has more effect when reg_par is small and
minimal impact with reg_par increases:

<img src="../figs/merra2_hyper_param_nh_quad.png" width="3520" />

This plot considers the effect of U_pi with the best performing other
hyperparameters. We see that the value of U_pi associated with the
lowest RMSE changes based on the value of reg_par:

<img src="../figs/merra2_hyper_param_nh_upi.png" width="3520" />

# Number of Models in Ensemble

**Note: These results are not included in the manuscript.**

Generate AOD and stratospheric temperature matrices for fitting ESNs
using all years (1980-1995):

``` r
train_mat_nens_aod <-
  merra2 |>
  select(date, lon, lat, aod_anom) |>
  pivot_wider(
    id_cols = c(lon, lat),
    names_from = date,
    values_from = aod_anom
  ) |>
  select(-lon,-lat) |>
  t()

train_mat_nens_temp_strat <-
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
eofs_nens_aod = compute_eofs(Ztrain = train_mat_nens_aod, n_eofs = n_eofs)
eofs_nens_temp_strat = compute_eofs(Ztrain = train_mat_nens_temp_strat, n_eofs = n_eofs)
phi_nens_train = eofs_nens_temp_strat$phi
```

Specify model inputs/outputs (AOD and stratospheric temperature are used
to forecast stratospheric temperature):

``` r
x_nens = cbind(eofs_nens_aod$train, eofs_nens_temp_strat$train)
y_nens = eofs_nens_temp_strat$train
```

Function for training an EESN and returning the RMSE. For a specified
number of ESNs per ensemble and seed:

1.  Train ensemble of ESNs
2.  Get predictions from each ensemble
3.  For each spatial location and time, average prediction across
    ensembles
4.  Compute weighted RMSE averaged over all times and locations

``` r
get_ensemble_preds <- function(n_ens, seed, nh, reg_par, add_quad) {
  
  # Train EESN
  eesn <-
    fit_Eesn(
      x = x_nens,
      y = y_nens,
      t = as.character(row.names(x_nens)),
      tau = 1,
      m = 5,
      tau_emb = 1,
      U_width = 0.1,
      U_pi = 0.1,
      W_width = 0.1, 
      W_pi = 0.1,
      nu = 0.2,
      nh = nh,
      reg_par = reg_par, 
      add_quad = add_quad,
      internal_scaling = "joint",
      n_ensembles = n_ens,
      seed = seed, 
      cores = 1
    )
  
  # Get predictions from EESN on spatial scale and join with weights and such
  preds <-
    map_df(
      .x = eesn,
      .f = function(x) {
        predict_esn(x, phi = phi_nens_train)$preds_ins |>
          data.frame() |>
          tibble::rownames_to_column(var = "date")
      },
      .id = "ensemble"
    ) |>
    pivot_longer(
      cols = -c(ensemble, date),
      names_to = "loc_id",
      values_to = "pred_anom"
    ) |>
    summarise(temp_strat_anom_pred = mean(pred_anom),
              .by = c(date, loc_id)) |>
    mutate(date = as_date(date),
           loc_id = as.numeric(str_remove(loc_id, "X")))
  
  # Compute weighted RMSE
  wgt_rmse <-
    preds |>
    left_join(merra2 |> select(date, loc_id, loc_weight, temp_strat_anom)) |>
    summarize(wgt_rmse =
                sqrt(sum(
                  loc_weight * (temp_strat_anom_pred - temp_strat_anom) ^ 2
                ) / sum(loc_weight))) |>
    mutate(n_ens, seed) |>
    select(n_ens, seed, wgt_rmse)
  
  # Return weighted RMSE
  return(wgt_rmse)
  
}
```

Train ESN ensembles with different number of base models and extract
weighted RMSEs:

``` r
set.seed(20230921)
n_ens = c(1, seq(25, 100, 25))
n_sims = 25
seeds = sample(100000:1000000, n_sims, replace = F)
sim_cases = expand_grid(n_ens, seeds)
  
fp_nens_sims_eesn1 = "../output/merra2_nens_sim_res_eesn1.csv"
if (!file.exists(fp_nens_sims_eesn1)) {
  plan(multisession, workers = 3)
  set.seed(20230921)
  nens_sims_eesn1 <- 
    future_map2_dfr(
      .x = sim_cases$n_ens, 
      .y = sim_cases$seeds, 
      .f = get_ensemble_preds,
      nh = 1000, 
      reg_par = 5, 
      add_quad = TRUE
    )
  write.csv(nens_sims_eesn1, file = fp_nens_sims_eesn1, row.names = FALSE)
}
nens_sims_eesn1 = read.csv(fp_nens_sims_eesn1)

fp_nens_sims_eesn2 = "../output/merra2_nens_sim_res_eesn2.csv"
if (!file.exists(fp_nens_sims_eesn2)) {
  plan(multisession, workers = 3)
  set.seed(20230921)
  nens_sims_eesn2 <- 
    future_map2_dfr(
      .x = sim_cases$n_ens, 
      .y = sim_cases$seeds, 
      .f = get_ensemble_preds,
      nh = 500, 
      reg_par = 0.1, 
      add_quad = FALSE
    )
  write.csv(nens_sims_eesn2, file = fp_nens_sims_eesn2, row.names = FALSE)
}
nens_sims_eesn2 = read.csv(fp_nens_sims_eesn2)

fp_nens_sims_eesn3 = "../output/merra2_nens_sim_res_eesn3.csv"
if (!file.exists(fp_nens_sims_eesn3)) {
  plan(multisession, workers = 3)
  set.seed(20230921)
  nens_sims_eesn3 <- 
    future_map2_dfr(
      .x = sim_cases$n_ens, 
      .y = sim_cases$seeds, 
      .f = get_ensemble_preds,
      nh = 2000, 
      reg_par = 50, 
      add_quad = FALSE
    )
  write.csv(nens_sims_eesn3, file = fp_nens_sims_eesn3, row.names = FALSE)
}
nens_sims_eesn3 = read.csv(fp_nens_sims_eesn3)
```

Plot results:

<img src="../figs/merra2_nens_comparison.png" width="1920" />
