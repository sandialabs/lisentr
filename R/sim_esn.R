#' Simulate data from an ESN model
#'
#' Function for simulating data from an ESN model.
#'
#' @param x Matrix of training data inputs with T rows (number of times) and 
#'        P columns (number of input variables) - Note: rows must be ordered 
#'        top to bottom from earliest to latest observation)
#' @param t Vector containing times/dates associated with x matrix
#' @param tau Forecast lead time
#' @param m Embedding vector length
#' @param tau_emb Embedding lag
#' @param nh Number of hidden units
#' @param V Ridge regression parameter matrix
#' @param sigma2 Ridge regression error variance
#' @param U U matrix in hidden stage
#' @param W W matrix in hidden stage
#' @param nu Scaling parameter for W
#' @param add_quad Indicates whether to add a quadratic term to the ridge
#'        regression
#' @param internal_scaling Specifies the type of standardization that should be 
#'        applied internally to x and y. This is common in ESN literature.
#'        Options are "joint" (the mean and standard deviation from all values 
#'        in the x/y matrix are computed after appropriate times are removed
#'        and used to center and scale the values in the x/y matrix) and "none"
#'        (no internal scaling is applied). Default is "none". It is best to 
#'        only use internal scaling if the model parameters were obtained by
#'        first training a model using fit_esn'.
#' @param y_mean If internal_scaling is 'joint', mean value of y used for 
#'        scaling.
#' @param y_sd If internal_scaling is 'joint', standard deviation of y used 
#'        for scaling.
#' 
#' @importFrom stats rnorm
#' 
#' @export sim_esn
#'
#' @examples 
#' # Extract data times
#' x_sim = sim$x_stdzd
#' t_sim = sort(unique(sim$time))
#' 
#' # Prepare training data matrices for ESN:
#' library(dplyr)
#' prep_mat_sim <- function(var_name) {
#'   var_name = paste0(var_name, "_stdzd")
#'   sim %>%
#'     select(easting, northing, time, all_of(var_name)) %>%
#'      tidyr::pivot_wider(
#'       id_cols = c(easting, northing),
#'       names_from = time,
#'       values_from = all_of(var_name)
#'      ) %>%
#'    select(-easting, -northing) %>%
#'    t()
#'  }
#'  sim_mats <- purrr::set_names(c("y", "x")) %>% purrr::map(.f = prep_mat_sim)
#'  
#'  # Compute EOFs
#'  n_eofs = 5
#'  sim_eofs = purrr::map(.x = sim_mats, .f = compute_eofs, n_eofs = n_eofs)
#'  
#'  # Specify model inputs/outputs
#'  x_sim = cbind(sim_eofs$x$train)
#'  y_sim = sim_eofs$y$train
#'  
#'  # Fit ESN
#'  esn <-
#'   fit_esn(
#'     x = x_sim,
#'     y = y_sim,
#'     t = as.character(t_sim),
#'     tau = 1,
#'     m = 5,
#'     tau_emb = 1,
#'     nh = 50,
#'     add_quad = TRUE,
#'     internal_scaling = "joint",
#'     seed = 20230223
#'   )
#'  
#'  # Simulate data from ESN 
#'  esn_sim_dat <-
#'   sim_esn(
#'     x = x_sim,
#'     t = t_sim,
#'     tau = esn$params_tuning$tau,
#'     m = esn$params_tuning$m,
#'     tau_emb = esn$params_tuning$tau_emb,
#'     nh = esn$params_tuning$nh,
#'     V = esn$param_est$V,
#'     sigma2 = esn$param_est$sigma2,
#'     W = esn$param_est$W,
#'     U = esn$param_est$U,
#'     nu = esn$params_tuning$nu,
#'     add_quad = esn$add_quad
#'   )

sim_esn <- function(x, t, tau, m, tau_emb, nh, V, sigma2, W, U, nu,
                    add_quad, internal_scaling = "none", 
                    y_mean = NULL, y_sd = NULL) {
  
  # Prepare inputs to model
  data_obj <- 
    create_data_obj_x(
      x = x,
      t = t, 
      tau = tau, 
      m = m, 
      tau_emb = tau_emb,
      internal_scaling = internal_scaling
    )
  
  # Prepare U_and_W object
  spectral_radius = abs(eigen(W, only.values = TRUE)$values[1])
  W_scaled = W * nu / spectral_radius
  U_and_W <- list(
    U = U,
    W = W,
    W_scaled = W_scaled,
    spectral_radius = spectral_radius
  )
  
  # Create the h matrix (set starting values at 0)
  h_start = rep(0, nh)
  h_obj <-
    create_h(
      U_and_W = U_and_W,
      data_obj = data_obj,
      add_quad = add_quad,
      nh = nh,
      h_start = h_start
    )
  
  # Simulate ridge regression errors
  ncols_eta = dim(V)[1] 
  nrows_eta = dim(h_obj$h)[2]
  eta_vec = rnorm(n = ncols_eta * nrows_eta, mean = 0, sd = sqrt(sigma2))
  eta = matrix(eta_vec, nrow = nrows_eta, ncol = ncols_eta)
  
  # Compute simulated y values
  y = t(V %*% h_obj$h) + eta
  
  # Remove scaling (if previously applied)
  if (internal_scaling == "joint") {
    y = (y * y_sd) + y_mean
  }
  
  # Determine times associated with y
  n_input_times = length(t)
  y_index = (tau + m * tau_emb + 1):(n_input_times)
  y_times = t[y_index]
  
  # Return
  return(list(
    y = y,
    y_times = y_times,
    data_obj = data_obj,
    h_obj = h_obj,
    eta = eta
  ))
  
}
