#' Fit Ensemble ESN model
#'
#' Function for training an ensemble ESN model.
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
#' @param n_ensembles Number of ensembles
#' @param cores Number of cores for parallelization (generally only needed
#'        for large or deep ESNs)
#'
#' @export fit_Eesn
#' 
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel 
#' @importFrom foreach foreach %dopar%
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
#'   fit_Eesn(
#'     x = x, 
#'     y = y, 
#'     t = t, 
#'     tau = 2, 
#'     m = 1, 
#'     tau_emb = 1, 
#'     nh = 50, 
#'     seed = 102034985,
#'     n_ensembles = 10,
#'     cores = 1
#'   )

fit_Eesn <- function(x, y, t, tau, m, tau_emb, nh, h_start = NULL, 
                     U_width = .10, W_width = .10, U_pi = .10, W_pi = .10, 
                     nu = .35, reg_par = .01, add_quad = TRUE, 
                     internal_scaling = "joint", seed, n_ensembles, cores = 1) {
  
  # Defining 'i' to get R check note of "no visible binding for 
  # global variable ‘i’" to go away
  i <- NULL
  
  if(cores > 1){
    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
    esnList <- foreach::foreach(i=1:n_ensembles,.packages="listenr") %dopar% {
      esn <-
        fit_esn(
          x = x,
          y = y,
          t = as.character(t),
          tau = tau,
          m = m,
          tau_emb = tau_emb,
          nh = nh,
          h_start = h_start,
          U_width = U_width,
          W_width = W_width,
          U_pi = U_pi,
          W_pi = W_pi,
          nu = nu,
          reg_par = reg_par,
          add_quad = add_quad,
          internal_scaling = internal_scaling,
          seed = seed*i
        )
      esn
    }
    parallel::stopCluster(cl)
  }
  else{
    esnList <- lapply(1:n_ensembles,function(i){
      esn <-
        fit_esn(
          x = x,
          y = y,
          t = as.character(t),
          tau = tau,
          m = m,
          tau_emb = tau_emb,
          nh = nh,
          h_start = h_start,
          U_width = U_width,
          W_width = W_width,
          U_pi = U_pi,
          W_pi = W_pi,
          nu = nu,
          reg_par = reg_par,
          add_quad = add_quad,
          internal_scaling = internal_scaling,
          seed = seed*i
        )
      return(esn)
    })
  }
  return(esnList)
}
