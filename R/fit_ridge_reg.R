#' Estimate ridge regression parameters
#'
#' This function returns the matrix of ridge regression parameters, V. 
#' 
#' @param h_obj Object output from `create_h`
#' @param add_quad Indicates whether to add a quadratic term to the linear ridge
#'        regression
#' @param reg_par Output ridge regression parameter
#' @param data_obj Output from create_embedded_data
#' @param nh Number of hidden units
#' @param internal_scaling Specifies the type of standardization that should be 
#'        applied internally to x and y. This is common in ESN literature.
#'        Options are "joint" (the mean and standard deviation from all values 
#'        in the x/y matrix are computed after appropriate times are removed
#'        and used to center and scale the values in the x/y matrix) and "none"
#'        (no internal scaling is applied). Default is "joint".
#' 
#' @export fit_ridge_reg

fit_ridge_reg <- function(h_obj, add_quad, reg_par, nh, data_obj,
                          internal_scaling) {
  
  # Set up matrix to use for ridge regression
  if (add_quad) {
    ridge_mat = reg_par * diag(2 * nh)
  } else {
    ridge_mat = reg_par * diag(nh)
  }
  
  # Select y matrix to be used for ridge regression (based on whether 
  # internal scaling is performed -- note, h is already created based on 
  # this specification)
  if (internal_scaling == "joint") {
    y_ridge = data_obj$y_train_scaled
  } else {
    y_ridge = data_obj$y_train
  }
  
  # Extract h matrix
  h = h_obj$h
  
  # Estimate ridge regression parameters
  V = t(y_ridge) %*% t(h) %*% solve(h %*% t(h) + ridge_mat)
  
  # Return estimated parameters
  return(V)
  
}