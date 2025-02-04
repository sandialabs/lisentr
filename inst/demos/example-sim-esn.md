Generation of Simulated Data
================

- [Data](#data)
- [Models](#models)
- [One Simulation](#one-simulation)
- [Multiple Simulations](#multiple-simulations)

Load R packages:

``` r
library(dplyr)
library(forcats)
library(ggplot2)
library(listenr)
library(purrr)
library(stringr)
library(tidyr)
```

# Data

Visualize covariate and response climatologies from simulated data:

![](example-sim-esn_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

Prepare training data matrices for ESN:

``` r
prep_mat_sim <- function(var_name) {
  var_name = paste0(var_name, "_stdzd")
  sim %>%
    select(easting, northing, time, all_of(var_name)) %>%
    pivot_wider(
      id_cols = c(easting, northing),
      names_from = time,
      values_from = all_of(var_name)
    ) %>%
    select(-easting, -northing) %>%
    t()
}
sim_mats <- set_names(c("y", "x")) %>% map(.f = prep_mat_sim)
```

Compute EOFs:

``` r
n_eofs <- 5
sim_eofs = map(.x = sim_mats, .f = compute_eofs, n_eofs = n_eofs)
```

Extract times:

``` r
t_sim = sort(unique(sim$time))
```

Specify model inputs/outputs:

``` r
x_sim = cbind(sim_eofs$x$train)
y_sim = sim_eofs$y$train
```

Plots of EOFs:

``` r
matplot(x_sim[, 1:5], type = 'l', main = "x", xlab = "time", ylab = "pc")
```

![](example-sim-esn_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

# Models

Fit ESNs (with all three scaling options):

``` r
esn_none <-
  fit_esn(
    x = x_sim,
    y = y_sim,
    t = as.character(t_sim),
    tau = 1,
    m = 5,
    tau_emb = 1,
    nh = 50,
    add_quad = TRUE,
    internal_scaling = "none",
    seed = 20230223
  )

esn_joint <-
  fit_esn(
    x = x_sim,
    y = y_sim,
    t = as.character(t_sim),
    tau = 1,
    m = 5,
    tau_emb = 1,
    nh = 50,
    add_quad = TRUE,
    internal_scaling = "joint",
    seed = 20230223
  )
```

# One Simulation

Simulate from the three ESNs:

``` r
set.seed(20230223)
esn_sim_none <-
  sim_esn(
    x = x_sim,
    t = t_sim,
    tau = esn_none$params_tuning$tau,
    m = esn_none$params_tuning$m,
    tau_emb = esn_none$params_tuning$tau_emb,
    nh = esn_none$params_tuning$nh,
    nu = esn_none$params_tuning$nu,
    add_quad = esn_none$add_quad,
    W = esn_none$param_est$W,
    U = esn_none$param_est$U,
    V = esn_none$param_est$V,
    sigma2 = esn_none$param_est$sigma2,
    internal_scaling = "none"
  )

set.seed(20230223)
esn_sim_joint <-
  sim_esn(
    x = x_sim,
    t = t_sim,
    tau = esn_joint$params_tuning$tau,
    m = esn_joint$params_tuning$m,
    tau_emb = esn_joint$params_tuning$tau_emb,
    nh = esn_joint$params_tuning$nh,
    nu = esn_joint$params_tuning$nu,
    add_quad = esn_joint$add_quad,
    W = esn_joint$param_est$W,
    U = esn_joint$param_est$U,
    V = esn_joint$param_est$V,
    sigma2 = esn_joint$param_est$sigma2, 
    internal_scaling = "joint",
    y_mean = esn_joint$data_train$y_train_mean,
    y_sd = esn_joint$data_train$y_train_sd
  )
```

Join observed and simulated data:

``` r
y_obs_none <- 
  data.frame(esn_none$data_train$y_train) %>%
  tibble::rownames_to_column(var = "t") %>%
  mutate(t = as.numeric(t)) %>%
  pivot_longer(cols = -t, names_to = "var") %>%
  mutate(scaling = "none")
  
y_obs_joint <- 
  data.frame(esn_joint$data_train$y_train) %>%
  tibble::rownames_to_column(var = "t") %>%
  mutate(t = as.numeric(t)) %>%
  pivot_longer(cols = -t, names_to = "var") %>%
  mutate(scaling = "joint")

y_obs <- 
  bind_rows(y_obs_none, y_obs_joint) %>%
  mutate(source = "data")

y_sim_none <- 
  data.frame(esn_sim_none$y) %>%
  mutate(t = esn_sim_none$y_times) %>%
  pivot_longer(cols = -t, names_to = "var") %>%
  mutate(scaling = "none")

y_sim_joint <- 
  data.frame(esn_sim_joint$y) %>%
  mutate(t = esn_sim_joint$y_times) %>%
  pivot_longer(cols = -t, names_to = "var") %>%
  mutate(scaling = "joint")

y_sim <- 
  bind_rows(y_sim_none, y_sim_joint) %>%
  mutate(source = "esn")

y_obs_sim <- 
  bind_rows(y_obs, y_sim) %>%
  mutate(var = str_replace(var, "X", "PC")) %>%
  select(source, scaling, t, var, value) %>%
  mutate(
    source = fct_relevel(source, "data", "esn"),
    scaling = fct_relevel(scaling, "none", "joint"))
```

Visualize EOFs and EOFs simulated from ESN:

![](example-sim-esn_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

# Multiple Simulations

Simulate multiple datasets from ESN:

``` r
nreps = 10

set.seed(20230227)
esn_sim_none_reps <-
  map(
    .x = rep(list(x_sim), nreps),
    .f = sim_esn,
    t = t_sim,
    tau = esn_none$params_tuning$tau,
    m = esn_none$params_tuning$m,
    tau_emb = esn_none$params_tuning$tau_emb,
    nh = esn_none$params_tuning$nh,
    nu = esn_none$params_tuning$nu,
    add_quad = esn_none$add_quad,
    W = esn_none$param_est$W,
    U = esn_none$param_est$U,
    V = esn_none$param_est$V,
    sigma2 = esn_none$param_est$sigma2
  )

set.seed(20230227)
esn_sim_joint_reps <-
  map(
    .x = rep(list(x_sim), nreps),
    .f = sim_esn,
    t = t_sim,
    tau = esn_joint$params_tuning$tau,
    m = esn_joint$params_tuning$m,
    tau_emb = esn_joint$params_tuning$tau_emb,
    nh = esn_joint$params_tuning$nh,
    nu = esn_joint$params_tuning$nu,
    add_quad = esn_joint$add_quad,
    W = esn_joint$param_est$W,
    U = esn_joint$param_est$U,
    V = esn_joint$param_est$V,
    sigma2 = esn_joint$param_est$sigma2, 
    internal_scaling = "joint",
    y_mean = esn_joint$data_train$y_train_mean,
    y_sd = esn_joint$data_train$y_train_sd
  )
```

Join simulated replicates:

``` r
esn_sim_none_reps_df <-
  map(.x = esn_sim_none_reps, .f = "y_times") %>%
  map_df(.f = data.frame) %>%
  rename('t' = '.x..i..') %>%
  bind_cols(map(.x = esn_sim_none_reps, .f = "y") %>%
              map_df(.f = data.frame, .id = "rep")) %>%
  pivot_longer(cols = c(-t,-rep), names_to = "var") %>%
  mutate(scaling = "none")

esn_sim_joint_reps_df <-
  map(.x = esn_sim_joint_reps, .f = "y_times") %>%
  map_df(.f = data.frame) %>%
  rename('t' = '.x..i..') %>%
  bind_cols(map(.x = esn_sim_joint_reps, .f = "y") %>%
              map_df(.f = data.frame, .id = "rep")) %>%
  pivot_longer(cols = c(-t,-rep), names_to = "var") %>%
  mutate(scaling = "joint")

y_sim_reps <- 
  bind_rows(
    esn_sim_none_reps_df, 
    esn_sim_joint_reps_df
  ) %>%
  mutate(var = str_replace(var, "X", "PC")) %>%
  select(scaling, rep, t, var, value) %>%
  mutate(scaling = fct_relevel(scaling, "none", "joint"))
```

Visualize EOFs and replicates of EOF simulated from ESN:

![](example-sim-esn_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->
