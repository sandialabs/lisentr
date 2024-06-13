#' Compute EOFs
#'
#' Function for computing EOFs given training and (optional) testing data. 
#' Test data principal component time series are computed using the 
#' training data EOFs.
#'
#' @param Ztrain Matrix of spatio-temporal training data with T (no. times) 
#'               rows and S (no. spatial locations) columns (rows must be
#'               ordered top to bottom from earliest to latest observation)
#' @param Ztest (Optional) Matrix of spatio-temporal testing data with T 
#'              (no. times) rows and S (no. spatial locations) columns (rows 
#'              must be ordered top to bottom from earliest to latest 
#'              observation) Set to NULL by default.
#' @param n_eofs Number of EOFs to return
#'
#' @export compute_eofs
#' 
#' @examples 
#' # Create training data
#' x = matrix(c(rnorm(12,10,1), rnorm(12,0,1)), ncol = 2, byrow = FALSE)
#' y = matrix(c(rnorm(12,5,3), rnorm(12,5,3)), ncol = 2, byrow = FALSE)
#' 
#' # Create a vector of times
#' t = c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
#'       "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
#' t = paste0(t, "2021")
#' 
#' # Create test data
#' x_new = matrix(c(rnorm(5,10,1), rnorm(5,0,1)), ncol = 2, byrow = FALSE)
#' t_new = c("Jan", "Feb", "Mar", "Apr", "May")
#' t_new = paste0(t_new, "2022")
#' 
#' # Compute EOFs
#' eofs <-
#'   compute_eofs(
#'     Ztrain = x,
#'     Ztest = x_new,
#'     n_eofs = 2
#'    )

compute_eofs <- function(Ztrain, Ztest = NULL, n_eofs) {
  
  # Perform SVD on training data
  E <- svd(Ztrain)
  
  # Extract specified number of EOF spatial basis functions from the 
  # training data:
  phi <- E$v[, 1:n_eofs]
  
  # Project training and testing (if specified) data onto basis functions to 
  # get the PC time series
  ts_train <- Ztrain %*% phi
  if (!is.null(Ztest)) {
    ts_test <- Ztest %*% phi 
  }
  
  # Return the PC time series (PC scores)
  if (is.null(Ztest)) {
    res = list(train = ts_train, phi = phi, svd_vals = E)
  } else {
    res = list(train = ts_train, test = ts_test, phi = phi, svd_vals = E)
  }
  return(res)
  
}
