#' Create h matrix
#'
#' This function prepares the h matrices to use in the ridge regression.
#'
#' @param U_and_W Output from create_U_and_W
#' @param data_obj Output from create_embedded_data
#' @param add_quad Indicates whether to add a quadratic term to the linear ridge
#'        regression
#' @param nh Number of hidden units
#' @param h_start Values to use for initializing h (vector of length nh)
#'
#' @export create_h

create_h <- function (U_and_W, data_obj, add_quad, nh, h_start) {
  
  # Multiply U by all (in sample) embedding vectors
  U_prod_mat = U_and_W$U %*% t(data_obj$design_matrix)
  
  # Determine the number of time periods (aka no. embedding vectors)
  n_time_periods = data_obj$n_emb_vecs
  
  # Set up the H matrix structure (based on whether to include quadratic term)  
  if (add_quad) {
    h = matrix(0, 2 * nh, n_time_periods)
  } else {
    h = matrix(0, nh, n_time_periods)
  }
  
  # Set the temporary x values to the starting values
  h_temp = h_start
  
  # Fill in the h matrix
  for (t in 1:n_time_periods) {
    h_temp = tanh(U_and_W$W_scaled %*% h_temp + U_prod_mat[, t])
    if (add_quad) {
      h[, t] = c(h_temp, h_temp ^ 2)
    } else {
      h[, t] = c(h_temp)
    }
  }
  
  # Return the h matrix
  return(list(h = h, h_end = h_temp, h_start = h_start))
  
}
