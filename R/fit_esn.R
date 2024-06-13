#' Fit ESN model
#'
#' Function for training an ESN model.
#'
#' @param x Matrix of training data inputs with T rows (number of times) and 
#'        P columns (number of input variables) - Note: rows must be ordered 
#'        top to bottom from earliest to latest observation)
#' @param y Matrix of training data outputs with T rows (number of times) and 
#'        Q columns (number of response variables) - Note: times associated
#'        with the rows must match those in x
#' @param t Vector containing times/dates associated with x and y matrices
#' @param tau Forecast lead time
#' @param m Embedding vector length
#' @param tau_emb Embedding lag
#' @param nh Number of hidden units
#' @param h_start Values to use for initializing h (vector of length nh)
#'        (If NULL, will be set to a vector of 0s)
#' @param U_width Uniform distribution sampling range for U ("width" parameter)
#'        (Default is .1)
#' @param W_width Uniform distribution sampling range for W ("width" parameter)
#'        (Default is .1)
#' @param U_pi Probabilities of non-zeros in U (sparseness parameter)
#'        (Default is .1)
#' @param W_pi Probabilities of non-zeros in W (sparseness parameter)
#'        (Default is .1)
#' @param nu Scaling parameter for W-weight matrix) (Default is .35)
#' @param reg_par Output ridge regression parameter (Default is .01)
#' @param add_quad Indicates whether to add a quadratic term to the linear ridge
#'        regression (Default is TRUE)
#' @param internal_scaling Specifies the type of standardization that should be 
#'        applied internally to x and y. This is common in ESN literature.
#'        Options are "joint" (the mean and standard deviation from all values 
#'        in the x/y matrix are computed after appropriate times are removed
#'        and used to center and scale the values in the x/y matrix) and "none"
#'        (no internal scaling is applied). Default is "joint".
#' @param seed Random seed
#'
#' @returns List of lists. (See details for more information.)
#' 
#' @details
#' Objects output from `fit_esn`:
#' \itemize{
#'   \item \code{data_input}: List containing the values of `x`, `y`, and `t` input to `fit_esn`
#'   \item \code{data_train}: List containing data objects created for the training of the model (using `create_data_obj`) including the embedding vectors and design matrix
#'   \item \code{params_tuning}: List of tuning parameter values input to `fit_esn`
#'   \item \code{param_est}: List of ESN parameters computed during model training
#'   \item \code{h}: List containing the matrix of h vectors (`h`), the values used to initiate the computation of h vectors (`h_start`), and the last h vector computed (`h_end`)
#'   \item \code{add_quad}: Indicator for whether a quadratic term was included in the model
#'   \item \code{internal_scaling}: Indicator for whether internal scaling was applied
#'   \item \code{seed}: Seed used for model training
#' }
#' 
#' @export fit_esn
#'
#' @examples 
#' # Create data
#' x = matrix(c(rnorm(12,10,1), rnorm(12,0,1)), ncol = 2, byrow = FALSE)
#' y = matrix(x[,1], ncol = 1)
#' 
#' # Assign column names to data
#' colnames(x) = c("X1", "X2")
#' colnames(y) = "X1"
#' 
#' # Create a vector of times
#' t = c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
#'       "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
#' t = paste0(t, "2021")
#' 
#' # Fit an ESN model to training data
#' esn <- 
#'   fit_esn(
#'     x = x, 
#'     y = y, 
#'     t = t, 
#'     tau = 2, 
#'     m = 1, 
#'     tau_emb = 1, 
#'     nh = 50, 
#'     seed = 1020349858
#'   )

fit_esn <- function(x, y, t, tau, m, tau_emb, nh, h_start = NULL, 
                    U_width = .10, W_width = .10, U_pi = .10, W_pi = .10, 
                    nu = .35, reg_par = .01, add_quad = TRUE, 
                    internal_scaling = "joint", seed) {
  
  # Prepare data for model training
  data_obj <- 
    create_data_obj(
      x = x, 
      y = y, 
      t = t, 
      tau = tau, 
      m = m, 
      tau_emb = tau_emb,
      internal_scaling = internal_scaling
    )
  
  # Create U and W matrices
  set.seed(seed)
  U_and_W <-
    compute_U_and_W(
      data_obj = data_obj,
      nh = nh,
      U_width = U_width,
      W_width = W_width,
      U_pi = U_pi,
      W_pi = W_pi,
      nu = nu
    )
  
  # Create the h matrix (set starting values at 0)
  if (is.null(h_start)) h_start = rep(0, nh)
  h_obj <-
    create_h(
      U_and_W = U_and_W,
      data_obj = data_obj,
      add_quad = add_quad,
      nh = nh,
      h_start = h_start
    )
  
  # Estimate ridge regression parameters
  V <- 
    fit_ridge_reg(
      h_obj = h_obj,
      add_quad = add_quad, 
      reg_par = reg_par, 
      nh = nh, 
      data_obj = data_obj,
      internal_scaling = internal_scaling
    )
  
  # Estimate ridge regression error variance
  yhat = t(V %*% h_obj$h)
  if (internal_scaling == "none") {
    sigma2 = sum((data_obj$y_train - yhat)^2) / length(yhat)
  } else{
    sigma2 = sum((data_obj$y_train_scaled - yhat)^2) / length(yhat)
  }
  
  # Put all items of interest in a list to return
  out <-
    list(
      data_input = list(x = x, y = y, t = t),
      data_train = data_obj,
      params_tuning = list(
        tau = tau,
        m = m,
        tau_emb = tau_emb,
        nh = nh,
        U_width = U_width,
        W_width = W_width,
        U_pi = U_pi,
        W_pi = W_pi,
        nu = nu,
        reg_par = reg_par
      ),
      param_est = 
        list(
          U = U_and_W$U, 
          W = U_and_W$W,
          W_scaled = U_and_W$W_scaled,
          spectral_radius = U_and_W$spectral_radius,
          V = V,
          sigma2 = sigma2
        ),
      h = h_obj,
      add_quad = add_quad,
      internal_scaling = internal_scaling,
      seed = seed
    )
  
  # Return all items
  return(out)
  
}
