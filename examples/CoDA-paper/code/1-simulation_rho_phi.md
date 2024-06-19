Simulation: Varying Rho and Phi
================
<br> Updated on April 17, 2024

**Overview**: This code runs a simulation study for the synthetic data
when varying `rhoX`, `rhoDelta`, `phiX`, and `phiDelta`.

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
library(doParallel)
library(listenr)
library(purrr)
library(tidyr)
```

Specify values for simulation study:

``` r
# PFI reps
nreps = 10

# Number of data sets per simulation
nsims = 50

# Number of time points and spatial locations
nT = 70
nx = 10
ny = 10

# Parameter values for data simulation
x1mean = 20
x2mean = 45
xsd = 6
sigmaX = 2
sigmaDelta = 2
sigmaEpsilon = 2
phiX = c(.2,.8)
phiDelta = c(.2,.8)
rhoX = c(.2,.8)
rhoDelta = c(.2,.8)

# Number of blocks to consider when computing feature importance
nblocks = c(1,2,3)

# Generate all combinations of parameters to use for the simulation study
sim_grid <-
  expand.grid(
    nblocks = nblocks,
    sigmaX = sigmaX,
    sigmaDelta = sigmaDelta,
    sigmaEpsilon = sigmaEpsilon,
    phiX = phiX,
    phiDelta = phiDelta,
    rhoX = rhoX,
    rhoDelta = rhoDelta
  )
```

Implement simulation study for comparing effects of varying the rhos and
phis:

``` r
fp_sim_phi_rho = "../output/simulation_res_phi_rho.csv"
if (!file.exists(fp_sim_phi_rho)) {
  cl <- makeCluster(3)
  registerDoParallel(cl)
  start <- Sys.time()
  iterStart <- Sys.time()
  iterEnd <- Sys.time()
  results <-
    foreach(
      i = 1:nrow(sim_grid),
      .packages = c("listenr", "dplyr", "purrr", "tidyr")
    ) %do% {
      # Print update on process
      print(paste0(
        "Iteration ",
        i,
        " out of ",
        nrow(sim_grid),
        ". Minutes of last iteration: ",
        round(difftime(iterEnd, iterStart, units = 'mins'), 2)
      ))
      
      # Specify start time for this iteration
      iterStart <- Sys.time()
      
      # Create place to store results
      fi_res_list <- list()
      
      # Simulate synthetic datasets
      # Includes simple trend + periodic trend + time-varying covariate effect
      # + change point in beta effect + time varying hidden covariate effect
      for (j in 1:nsims) {
        dat = simulate_synthetic_data(
          nT = nT,
          nx = nx,
          ny = ny,
          x1mean = x1mean,
          x2mean = x2mean,
          xsd = xsd,
          sigmaX = sim_grid$sigmaX[i],
          sigmaDelta = sim_grid$sigmaDelta[i],
          sigmaEpsilon = sim_grid$sigmaEpsilon[i],
          phiX = sim_grid$phiX[i],
          phiDelta = sim_grid$phiDelta[i],
          rhoX = sim_grid$rhoX[i],
          rhoDelta = sim_grid$rhoDelta[i],
          seed = 893 + 82 * j + 91
        )
        
        # Add date variable (based on time)
        ta6_raw = dat |> mutate(date = time)
        
        # Create location ID associated with longitude and latitude:
        ta6_locs <-
          ta6_raw |>
          select(easting, northing, date, val) |>
          pivot_wider(
            id_cols = c(easting, northing),
            names_from = date,
            values_from = val
          ) |>
          mutate(loc_id = 1:n()) |>
          select(loc_id, easting, northing)
        
        # Add train/test labels:
        ta6 <-
          ta6_raw |>
          mutate(data = ifelse(date > 90, "test", "train")) |>
          select(data, everything())
        
        # Compute training data spatial means
        ta6_spatial_means <-
          ta6 |>
          filter(data == "train") |>
          group_by(easting, northing) |>
          summarise_at(
            .vars = vars(val, x1, x2),
            .funs = mean,
            .groups = "drop"
          ) |>
          rename_at(
            .vars = vars(val, x1, x2),
            .funs = function(var)
              paste0(var, "_mean"),
          )
        
        # Compute training data spatial SDs
        ta6_spatial_sds <-
          ta6 |>
          filter(data == "train") |>
          group_by(easting, northing) |>
          summarise_at(.vars = vars(val, x1, x2),
                       .funs = stats::sd) |>
          ungroup() |>
          rename_at(
            .vars = vars(val, x1, x2),
            .funs = function(var)
              paste0(var, "_sd"),
          )
        
        # Center and scale variables based on training data:
        ta6_esn <-
          left_join(ta6, ta6_spatial_means, by = join_by(easting, northing)) |>
          left_join(ta6_spatial_sds, by = join_by(easting, northing)) |>
          mutate(
            val_stdzd = (val - val_mean) / val_sd,
            x1_stdzd = (x1 - x1_mean) / x1_sd,
            x2_stdzd = (x2 - x2_mean) / x2_sd
          )
        
        # Extract train/test data
        ta6_train = ta6_esn |> filter(data == "train")
        
        # Prepare training data matrices for ESN
        prep_train_mat_ta6 <- function(var_name) {
          var_name = paste0(var_name, "_stdzd")
          ta6_train |>
            select(easting, northing, date, var_name) |>
            pivot_wider(
              id_cols = c(easting, northing),
              names_from = date,
              values_from = var_name
            ) |>
            select(-easting, -northing) |>
            t()
        }
        ta6_train_mats <-
          set_names(c("val", "x1", "x2")) |>
          map(.f = prep_train_mat_ta6)
        
        # Compute EOFs
        n_eofs <- 5
        ta6_eofs = map(.x = ta6_train_mats,
                       .f = compute_eofs,
                       n_eofs = n_eofs)
        
        # Extract times for training and testing data
        t_train_ta6 = sort(unique(ta6_train$date))
        
        # Specify model inputs/outputs
        x_ta6 = cbind(ta6_eofs$x1$train, ta6_eofs$x2$train)
        y_ta6 = ta6_eofs$val$train
        
        # Fit ESN
        esn <-
          fit_esn(
            x = x_ta6,
            y = y_ta6,
            t = as.character(t_train_ta6),
            tau = 1,
            m = 1,
            tau_emb = 1,
            nh = 50,
            add_quad = TRUE,
            internal_scaling = "none",
            seed = 20230411,
            reg_par = .1
          )
        
        # Compute PFI without retraining on spatial scale
        pfi_wo_retrain_spatial <-
          compute_pfi(
            model = esn,
            nreps = nreps,
            var_groups = list(1:n_eofs, (n_eofs + 1):(2 * n_eofs)),
            y_spatial = ta6_train_mats$val,
            phi = ta6_eofs$val$phi,
            blockSize = sim_grid$nblocks[i],
            seed = 20230411
          )
        
        # Clean up PFI results and join
        pfi_res <-
          bind_rows(
            pfi_wo_retrain_spatial$pfi |>
              mutate(retraining = "no retraining", spatial = "spatial scale") |>
              rename(pfi = pfi_mean),
          ) |>
          tibble::remove_rownames() |>
          rename(importance = pfi,
                 vars_adj = vars_perm,
                 t_adj = t_perm) |>
          mutate(
            feat_imp = "PFI",
            variable = ifelse(vars_adj == paste(1:n_eofs, collapse = ","), "x1", "x2"),
            t_forecasted = as.numeric(t_forecasted),
            t_adj = as.numeric(t_adj)
          ) |>
          mutate(time_diff = t_forecasted - t_adj) |>
          filter(time_diff == 1) |>
          select(
            feat_imp,
            retraining,
            spatial,
            variable,
            vars_adj:t_forecasted,
            time_diff,
            everything()
          )
        
        # Compute ZFI without retraining on spatial scale
        zfi_wo_retrain_spatial <-
          compute_zfi(
            model = esn,
            var_groups = list(1:n_eofs, (n_eofs + 1):(2 * n_eofs)),
            y_spatial = ta6_train_mats$val,
            phi = ta6_eofs$val$phi,
            blockSize = sim_grid$nblocks[i],
            seed = 20230411
          )
        
        # Clean up ZFI results and join
        zfi_res <-
          bind_rows(
            zfi_wo_retrain_spatial$zfi |>
              mutate(retraining = "no retraining", spatial = "spatial scale")
          ) |>
          rename(importance = zfi,
                 vars_adj = vars_zeroed,
                 t_adj = t_zeroed) |>
          mutate(
            feat_imp = "ZFI",
            variable = ifelse(vars_adj == paste(1:n_eofs, collapse = ","), "x1", "x2"),
            t_forecasted = as.numeric(t_forecasted),
            t_adj = as.numeric(t_adj)
          ) |>
          mutate(time_diff = t_forecasted - t_adj) |>
          filter(time_diff == 1) |>
          select(
            feat_imp,
            retraining,
            spatial,
            variable,
            vars_adj:t_forecasted,
            time_diff,
            everything()
          )
        
        # Feature Importance Results
        fi_res_list[[j]] <-
          bind_rows(pfi_res, zfi_res)#pfi_res#zfi_res#
        fi_res_list[[j]]$sim <- j
        
      }
      
      # Organize outputs
      out <- do.call("rbind", fi_res_list)
      out$sigmaX = sim_grid$sigmaX[i]
      out$sigmaDelta = sim_grid$sigmaDelta[i]
      out$sigmaEpsilon = sim_grid$sigmaEpsilon[i]
      out$phiX = sim_grid$phiX[i]
      out$phiDelta = sim_grid$phiDelta[i]
      out$rhoX = sim_grid$rhoX[i]
      out$rhoDelta = sim_grid$rhoDelta[i]
      out$nblocks = sim_grid$nblocks[i]
      out$simGroup = i
      
      iterEnd <- Sys.time()
      out
      
    }
  
  # Get end time and compute amount of time
  finish <- Sys.time()
  timeDiff <- finish - start
  timeDiff
  
  # End parallel processing
  stopCluster(cl)
  
  # Join results
  resultsdf <- do.call("rbind", results)
  
  # Save results
  write.csv(resultsdf, file = fp_sim_phi_rho)
}
```
