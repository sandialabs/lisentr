#' Function for computing FI with ESN
#'
#' Computes FI using the training data specified in fit_esn
#'
#' @param model Object output from fit_esn or fit_Eesn function
#' @param type Feature importance type to be computed ('pfi' or 'zfi')
#' @param nreps Number of reps to use when computing PFI.
#' @param var_groups List of vectors where each vectors contains columns numbers
#'        indicating which x variables should be grouped when doing the
#'        zeroing process (set to NULL by default which zeros all x variables)
#' @param blockSize Total number of time points to zero out, defaults to 1
#' @param y_spatial Matrix of observed response values on the spatial scale with
#'        rows corresponding to times and columns corresponding to spatial 
#'        locations (should be the same matrix used for Ztrain in compute_eofs).
#' @param phi Matrix object phi output from compute_eofs when applied to Ztrain.
#' @param weights Vector of weights with length corresponding to the number of 
#'        columns in y (or y_spatial, if y_spatial is specified) Note: currently
#'        assumes that the weights are equivalent across time -- only varies by
#'        location.
#' @param scale_y Indicates whether y and corresponding predictions should be
#'        scaled before computing RMSEs. Scaling y values and predictions is 
#'        intended to put multiple response variables on the same scale 
#'        for comparison in the computation of FI. (Default is FALSE.)
#' @param return_adj_preds Indicates whether the permuted/zeroed predictions are
#'        returned in addition to PFI/ZFI values. Default is FALSE.
#' @param seed Random seed (NULL by default.)
#'
#' @export compute_fi
#'
#' @importFrom dplyr filter full_join group_by left_join mutate rename
#'             select starts_with summarise
#' @importFrom lubridate as_date
#' @importFrom purrr map map_depth map_df
#' @importFrom tibble remove_rownames
#' @importFrom tidyr pivot_longer pivot_wider
#' 
#' @examples 
#' # Create data
#' x = matrix(rnorm(40,0,1), ncol = 4)
#' y = matrix(rnorm(20,0,1), ncol = 2)
#' 
#' # Assign column names
#' colnames(x) = c("x11", "x12", "x21", "x22")
#' colnames(y) = c("y1", "y2")
#' 
#' # Times
#' t = paste0("t", 1:10)
#' 
#' # Fit an ESN model
#' esn <- 
#'   fit_esn(
#'     x = x, 
#'     y = y, 
#'     t = t, 
#'     tau = 2, 
#'     m = 1, 
#'     tau_emb = 1, 
#'     nh = 10, 
#'     seed = 1020349858
#'   )
#' 
#' # Compute ZFI
#' zfi <- 
#'   compute_fi(
#'     model = esn, 
#'     type = "zfi",
#'     var_groups = list(1:2, 3:4), 
#'     seed = 10203498,
#'     blockSize = 1
#'    )
#'    
#' # Compute ZFI
#' pfi <- 
#'   compute_fi(
#'     model = esn, 
#'     type = "pfi",
#'     nreps = 2,
#'     var_groups = list(1:2, 3:4), 
#'     seed = 10203498,
#'     blockSize = 1
#'    )

compute_fi <- function(model, type, nreps = NULL, var_groups = NULL, 
                       blockSize = 1, y_spatial = NULL, phi = NULL, 
                       weights = NULL, scale_y = FALSE, 
                       return_adj_preds = FALSE, seed = NULL) {
  
  #### CHECKS ####
  
  # Make code compatible with output from fit_esn too
  if(!is.null(names(model))){
    model <- list(model)
  }
  
  # Check that nreps is appropriate
  if (type == "zfi" & !is.null(nreps)) {
    print("Note: nreps must be set to NULL when 'type' is 'zfi' and has 
          been adjusted accordingly.")
  }
  if (type == "zfi") nreps = 1
  
  #### SET UP FOR FI ####
  
  # If var_groups is NULL, set var_groups to contain all x vars
  if (is.null(var_groups)) {
    var_groups = list(1:ncol(model[[1]]$data_input$x))
  }
  
  # Extract times associated with x and y
  x_train_times = model[[1]]$data_train$x_train_times
  y_train_times = model[[1]]$data_train$y_train_times
  
  # Identify first time that can be adj
  x_index_pred_start = (model[[1]]$params_tuning$m * model[[1]]$params_tuning$tau_emb) + 1
  x_index_adj_start = 
    ifelse(blockSize > x_index_pred_start, blockSize, x_index_pred_start)
  x_index_end = length(x_train_times)
  if (blockSize > x_index_pred_start) {
    print("Note: Block size is larger than first prediction, so some times will be excluded from feature importance computation.")
  }
  
  # Create data frame containing labels to attach with FIs
  y_t_pred_start = x_index_adj_start + model[[1]]$params_tuning$tau
  y_index_pred_start = 
    y_t_pred_start - 
    model[[1]]$params_tuning$tau - 
    (model[[1]]$params_tuning$m * model[[1]]$params_tuning$tau_emb)
  
  #### RMSES ON OBSERVED DATA ####
  
  # Extract observed in-sample y values or if specified, use y_spatial (with 
  # times removed according)
  if (is.null(y_spatial)) {
    y_obs = model[[1]]$data_train$y_train
  } else {
    y_obs <-
      create_data_obj_y(
        y = y_spatial,
        t = model[[1]]$data_input$t,
        tau = model[[1]]$params_tuning$tau,
        m = model[[1]]$params_tuning$m,
        tau_emb = model[[1]]$params_tuning$tau_emb,
        internal_scaling = "none"
      )$y_train
  }
  
  # Compute predictions on in-sample data (will compute on spatial 
  # scale if phi is not NULL)
  y_pred_temp <- lapply(
    1:length(model), 
    function(x) predict_esn(model = model[[x]], phi = phi)$preds_ins
  )
  y_pred = Reduce("+", y_pred_temp) / length(y_pred_temp)
  
  # If specified, scale obs and pred y
  if (scale_y) {
    y_obs = scale(y_obs)
    y_pred = scale(y_pred)
  }
  
  # Compute RMSEs on predictions from observed data
  if (is.null(weights)) {
    rmses_obs = sqrt(rowMeans((y_obs - y_pred)^2))  
  } else {
    rmses_obs = sqrt(rowSums(weights * (y_obs - y_pred)^2) / sum(weights))
  }
  
  #### RMSES ON ADJUSTED DATA ####
  
  # Generate cases for when to compute RMSE
  rmse_adj_cases <-
    expand.grid(
      index_adj = x_index_adj_start:x_index_end,
      var_group = var_groups,
      rep = 1:nreps
    )
  
  # Adjusted RMSES and predictions (if requested)
  set.seed(seed)
  rmses_adj <-
    purrr::pmap_df(
      .l = as.list(rmse_adj_cases), 
      .f = compute_adj_rmses,
      type = type,
      y_obs = y_obs,
      model = model,
      blockSize = blockSize,
      y_spatial = y_spatial, 
      phi = phi, 
      weights = weights,
      scale_y = scale_y,
      return_adj_preds = return_adj_preds
    ) %>%
    tibble::remove_rownames()
  
  #### FI COMPUTATIONS ####
  
  # Compute FI values and clean up data frame
  n_rep_groups = length(var_groups)*nreps
  fi <-
    rmses_adj %>%
    dplyr::mutate(
      t_forecasted = rep(
        y_train_times[y_index_pred_start:length(y_train_times)], 
        n_rep_groups
      ),
      rmses_obs = rep(
        rmses_obs[y_index_pred_start:length(y_train_times)], 
        n_rep_groups
      )
    ) %>%
    dplyr::mutate(fi = -.data$rmses_obs - (-.data$rmses_adj)) %>%
    dplyr::select(
      .data$t_adj,
      .data$t_forecasted,
      .data$vars_adj,
      .data$rep,
      .data$rmses_obs,
      .data$rmses_adj, 
      .data$fi,
      dplyr::everything()
    )
  
  # Adjust results based on whether ZFI or PFI
  if (type == "zfi") {
    fi <- fi %>% dplyr::select(-rep)
  }
  
  # Return fi results
  return(fi)
  
}
