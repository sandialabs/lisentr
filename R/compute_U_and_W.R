#' Create U and W matrices
#'
#' This function creates the U and W matrices for the model hidden state
#' 
#' @param data_obj Output from create_embedded_data
#' @param nh Number of hidden units
#' @param U_width Uniform distribution sampling range for U ("width" parameter)
#' @param W_width Uniform distribution sampling range for W ("width" parameter)
#' @param U_pi Probabilities of non-zeros in U (sparseness parameter)
#' @param W_pi Probabilities of non-zeros in W (sparseness parameter)
#' @param nu Spectral-radius parameter (scaling parameter for W-weight matrix)
#'        
#' @importFrom stats rbinom runif
#' @export compute_U_and_W

compute_U_and_W <- function(data_obj, nh, U_width, W_width, U_pi, W_pi, nu) {
  
  # Create U
  n_cols_U = ncol(data_obj$design_matrix)
  U = matrix(runif(nh * n_cols_U, -U_width, U_width), nh)
  for (i in 1:n_cols_U) {
    numNonZeroU = rbinom(1, nh, U_pi)
    numZeroU = nh - numNonZeroU
    curUIndexZero = sample(1:nh, numZeroU)
    U[curUIndexZero, i] = 0
  }
  
  # Create W
  W = matrix(runif(nh * nh, -W_width, W_width), nh)
  for (i in 1:nh) {
    numNonZeroW = rbinom(1, nh, W_pi)
    numZeroW = nh - numNonZeroW
    curWIndexZero = sample(1:nh, numZeroW)
    W[curWIndexZero, i] = 0
  }
  
  # Scale W matrix
  spectral_radius = abs(eigen(W, only.values = TRUE)$values[1])
  W_scaled = W * nu / spectral_radius
  
  # Return U and W
  return(list(
    U = U,
    W = W,
    W_scaled = W_scaled,
    spectral_radius = spectral_radius
  ))
  
}