#' Fit Ensembled (bagged) ESN model
#'
#' Function for training an Ensembled (bagged) ESN model.
#'
#' @param x Matrix of training data inputs with T rows (number of times) and 
#'        P columns (number of input variables). Rows must be ordered 
#'        top to bottom from earliest to latest observation.
#' @param y Matrix of training data outputs with T rows (number of times) and 
#'        Q columns (number of response variables). Times associated
#'        with the rows must match those in x.
#' @param t Vector containing times/dates associated with x and y matrices.
#' @param x_test Matrix of testing data inputs with T' rows (number of times) 
#'        and P columns (number of input variables). Rows must be ordered 
#'        top to bottom from earliest to latest observation.
#' @param t_test Vector containing times/dates associated with x_test.
#' @param phi Matrix object phi output from compute_eofs when applied to Ztrain.
#' @param obs_train Isn't this the same as x?
#' @param obs_test Isn't this the same as x_test?
#' @param tau Forecast lead time
#' @param m Embedding vector length
#' @param tau_emb Embedding lag
#' @param nh Number of hidden units
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
#' @param add_quad Indicates whether to add a quadratic term to the linear 
#'        ridge regression (Default is TRUE)
#' @param internal_scaling Specifies the type of standardization that should be 
#'        applied internally to x and y. This is common in ESN literature.
#'        Options are "joint" (the mean and standard deviation from all values 
#'        in the x/y matrix are computed after appropriate times are removed
#'        and used to center and scale the values in the x/y matrix) and "none"
#'        (no internal scaling is applied). Default is "joint".
#' @param weights Vector of weights for computing weighted RMSEs (if NULL, 
#'        RMSEs are computed without weighting)
#' @param seed Random seed
#' @param n_ensembles number of ensembles
#' @param cores number of cores for parallelization (generally only needed for
#'        large or deep ESNs)
#'
#' @export hyperparameter_search
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
#' # Assign column names to data
#' colnames(x) = c("X1", "X2")
#' colnames(y) = c("Y1", "Y2")
#' colnames(x_new) = c("X1", "X2")
#' 
#' # Compute EOFs
#' x_eofs <-
#'   compute_eofs(
#'     Ztrain = x,
#'     Ztest = x_new,
#'     n_eofs = 2
#'    )
#' y_eofs <-
#'   compute_eofs(
#'     Ztrain = y,
#'     n_eofs = 2
#'    )
#'    
#' # Implement hyperparameter search
#' rmses <- 
#'   hyperparameter_search(
#'     x = x_eofs$train,
#'     y = y_eofs$train,
#'     t = t,
#'     x_test = x_eofs$test,
#'     t_test = t_new,
#'     phi = x_eofs$phi,
#'     obs_train = x,
#'     obs_test = x_new,
#'     tau = 1,
#'     m = 1, 
#'     tau_emb = 1, 
#'     nh = 50, 
#'     U_width = 0.1,
#'     W_width = 0.1, 
#'     U_pi = 0.1, 
#'     W_pi = 0.1, 
#'     nu = c(0.5,1), 
#'     reg_par = c(1,2), 
#'     add_quad = TRUE, 
#'     internal_scaling = "joint", 
#'     seed = 1020349858
#'   )

hyperparameter_search <- function(x, y, t, x_test, t_test, phi, obs_train, 
                                  obs_test, tau, m, tau_emb, nh, U_width, 
                                  W_width, U_pi, W_pi, nu, reg_par, 
                                  add_quad, internal_scaling = "joint", 
                                  weights = NULL, seed, n_ensembles = 1, 
                                  cores = 1) {
  
  # Make sure input values are matrices
  if(!is.matrix(obs_train)){
    obs_train = matrix(obs_train,ncol=1)
  }
  if(!is.matrix(obs_test)){
    preds_test = matrix(obs_test,ncol=1)
  }
  if(!is.matrix(y)){
    y = matrix(y,ncol=1)
  }
  
  # Defining 'i' to get R check note of "no visible binding for 
  # global variable ‘i’" to go away
  i <- NULL
  
  # Stop computation, if tau_emb is not equal to 1  
  if (sum(tau_emb) != 1) 
    stop("hyperparameter_search only currently works for tau_emb=c(1)")
  
  # Create grid of all possible tuning parameter combinations
  output <- 
    data.frame(
      expand.grid(
        tau, m, tau_emb, nh, U_width, W_width, U_pi, W_pi, nu, reg_par, add_quad
      )
    )
  names(output) <- c("tau", "m", "tau_emb", "nh", "U_width", 
                     "W_width", "U_pi", "W_pi", "nu", "reg_par", "add_quad")
  
  # Compute RMSEs for all ESN tuning parameter combinations
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  rmse <- foreach(i = 1:nrow(output), .packages = c("listenr")) %dopar% {
    
    # Prepare matrices for storing RMSEs
    rmse_train <- matrix(NA, ncol = n_ensembles, 
                         nrow = nrow(x) - output$tau[i] - output$m[i])
    rmse_test <- matrix(NA, ncol = n_ensembles, nrow = nrow(x_test))
    
    # Iterate of number of ensemble members
    for (j in 1:n_ensembles) {
      
      # Train model
      esn <- fit_esn(
        x = x, 
        y = y, 
        t = as.character(t), 
        tau = output$tau[i], 
        m = output$m[i], 
        tau_emb = output$tau_emb[i], 
        nh = output$nh[i], 
        U_width = output$U_width[i], 
        W_width = output$W_width[i], 
        U_pi = output$U_pi[i], 
        W_pi = output$W_pi[i], 
        nu = output$nu[i],
        reg_par = output$reg_par[i], 
        add_quad = output$add_quad[i], 
        internal_scaling = internal_scaling, 
        seed = j * seed + 34
      )
      
      # Get predictions
      preds_train = predict_esn(model = esn)$preds_ins
      preds_test = predict_esn(model = esn, x_new = x_test, 
                               t_new = t_test)$preds_new
      
      # Convert to matrix if needed
      if(!is.matrix(preds_train)){
        preds_train = matrix(preds_train,ncol=1)
      }
      if(!is.matrix(preds_test)){
        preds_test = matrix(preds_test,ncol=1)
      }
      
      # Convert to spatial scale if requested
      if (!is.null(phi)) {
        preds_train = preds_train %*% t(phi)
        preds_test = preds_test %*% t(phi)
      }
      
      # Compute RMSES
      if (is.null(weights)) {
        rmse_train[, j] =
          sqrt(rowMeans((preds_train - obs_train[-c(1:(output$tau[i] + output$m[i])),]) ^
                          2))
        rmse_test[, j] = sqrt(rowMeans((preds_test - obs_test) ^ 2))
      } else {
        rmse_train[, j] =
          sqrt(rowSums(weights * (preds_train - obs_train[-c(1:(output$tau[i] + output$m[i])),]) ^
                         2) / sum(weights))
        rmse_test[, j] =
          sqrt(rowSums(weights * (preds_test - obs_test) ^ 2) / sum(weights))
      }
      
    }
    
    # Join RMSEs
    list(rmse_train = rmse_train, rmse_test = rmse_test)
    
  }
  stopCluster(cl)
  
  # Join all results
  output$mRMSE_tr <- sapply(1:length(rmse), function(i) mean(rmse[[i]]$rmse_train))
  output$mRMSE_te <- sapply(1:length(rmse), function(i) mean(rmse[[i]]$rmse_test))
  
  # Return results
  outList <- list(output = output, rmseFull = rmse)
  return(outList)
  
}
