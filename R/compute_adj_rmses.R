#' Function for getting "adjusted" RMSEs
#'
#' Sets the variables values for the specified time to 0 and obtains ESN
#' predictions with the "adjusted" times and computes RMSEs
#'
#' @param index_adj Integer specifying which time should be permuted or
#'        set to 0
#' @param var_group Vectors of integers indicating the numbers of columns
#'        that should be included in the permuting/zeroing process
#' @param rep Used only for keeping track of reps associated with 
#'        `compute_fi`.
#' @param type Feature importance type to be computed ('pfi' or 'zfi')
#' @param y_obs Matrix of y values as computed internally in compute_zfi
#' @param model Model output from fit_esn function
#' @param blockSize number of time points to permute or zero out
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
#'        for comparison in the computation of PFI. (Default is FALSE.)
#' @param return_adj_preds Indicates whether the permuted/zeroed predictions are
#'        returned in addition to PFI/ZFI values. Default is FALSE.
#' 
#' @importFrom magrittr %>%
#' 
#' @export compute_adj_rmses
#'
#' @importFrom dplyr .data everything select
#' @importFrom stringr str_remove

compute_adj_rmses <- function(index_adj, var_group, rep, type, y_obs, model, 
                              blockSize, y_spatial, phi, weights, scale_y, 
                              return_adj_preds) {
  
  # Determine number of ESN models in ensemble
  n_esn = length(model)
  
  # Extract observed x training values
  x = model[[1]]$data_input$x
  
  # Determine number of x variables to be permuted or set to zero
  n_xvars = length(var_group)
  
  # Permute or set specified time(s) and variables to 0
  x_adj = x
  index_adj_init = index_adj - blockSize + 1
  indices_adj = index_adj_init:index_adj
  if (type == "zfi") {
    x_adj[indices_adj,var_group] = rep(0, n_xvars)
  } else if (type == "pfi") {
    for(i in indices_adj){
      x_adj[i,var_group] <-
        sample(x = x[i,var_group], size = n_xvars, replace = F)
    }
  } else {
    stop("'type' specified incorrectly. Must be 'zfi' or 'pfi'.")
  }
  
  # Prepare data for ESN with adjusted times
  data_obj_adj <-
    create_data_obj_ood(
      model = model[[1]], 
      x_ood = x_adj, 
      t_ood = model[[1]]$data_input$t
    )
  
  # Obtain adjusted predictions
  preds_temp <- list()
  for(i in 1:n_esn){
    
    # Compute h matrix with adjusted values
    h_obj_adj <-
      create_h(
        U_and_W = model[[i]]$param_est,
        data_obj = data_obj_adj,
        add_quad = model[[i]]$add_quad,
        nh = model[[i]]$params_tuning$nh,
        h_start = rep(0, model[[i]]$params_tuning$nh)
      )
    
    # Fill in model with adjusted values
    model_adj = model[[i]]
    model_adj$data_train$x_train = data_obj_adj$x_ood
    model_adj$data_train$x_train_scaled = data_obj_adj$x_ood_scaled
    model_adj$data_train$x_oos = data_obj_adj$x_ood_oos
    model_adj$data_train$emb_vecs = data_obj_adj$emb_vecs
    model_adj$data_train$design_matrix = data_obj_adj$design_matrix
    model_adj$h = h_obj_adj
    
    # Specify t_forecast and associated index
    y_index_pred = index_adj + model[[i]]$params_tuning$tau
    
    # Determine index in h based on forecast time
    y_train_times_index = 
      y_index_pred - 
      model[[i]]$params_tuning$tau - 
      (model[[i]]$params_tuning$m * model[[i]]$params_tuning$tau_emb)
    
    # Compute predictions
    preds_temp[[i]] = t(model_adj$param_est$V %*% model_adj$h$h[,y_train_times_index])
    
  }
  preds_adj <- Reduce("+", preds_temp) / length(preds_temp)
  
  # Remove scaling (if previously applied)
  if (model[[1]]$internal_scaling == "joint") {
    y_train_mean = model_adj$data_train$y_train_mean
    y_train_sd = model_adj$data_train$y_train_sd 
    preds_adj = (preds_adj * y_train_sd) + y_train_mean
  }
  
  # Add times to predictions
  rownames(preds_adj) = model_adj$data_train$y_train_times[y_train_times_index]
  
  # Convert back to spatial scale is requested
  if (!is.null(phi)) {
    preds_adj = preds_adj %*% t(phi)
  }
  
  # Scale predictions if requested
  if (scale_y) preds_adj = scale(preds_adj)
  
  # Make sure y_obs is a matrix
  if(!is.matrix(y_obs)){
    y_obs = matrix(y_obs,ncol=1)
  }
  
  # Compute RMSEs with adjusted data
  if (is.null(weights)) {
    rmses_adj = sqrt(rowMeans((y_obs[y_train_times_index,] - preds_adj)^2))  
  } else {
    rmses_adj = 
      sqrt(
        rowSums(weights * (y_obs[y_train_times_index,] - preds_adj)^2) / 
          sum(weights)
      )
  }
  
  # Put RMSEs in data frame
  res <- 
    data.frame(
      t_adj = model[[1]]$data_train$x_train_times[index_adj], 
      vars_adj = paste(var_group, collapse = ","),
      rep = rep, 
      rmses_adj = rmses_adj
    )
  if (return_adj_preds) {
    res = cbind(res, preds_adj)
  }
  
  # Return adjusted RMSEs
  return(res)
  
}
