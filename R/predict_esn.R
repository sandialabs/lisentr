#' Obtain predictions/forecasts with ESN model
#'
#' This function computes predictions for a model trained using `fit_esn`. Both
#' predictions on the training data and out-of-sample forecasts on new data can
#' be computed. New data must sequentially follow training data (no gaps in
#' time).
#'
#' @param model Object output from fit_esn function
#' @param x_ood Out of distribution data inputs (must be the same size matrix
#'        as training data x and occur at the same times as x)
#' @param x_new Data input matrix occurring after the times in the training
#'        data x
#' @param t_new Vector of times associated with x_new matrix
#' @param phi Matrix object phi output from compute_eofs when applied to 
#'        Ztrain. If specified, predictions will be returned on spatial scale.
#'        Default is NULL.
#' 
#' @returns List. (See details for more information.)
#' 
#' @details
#' Objects output from `predict_esn` when no new `x` values are input:
#' \itemize{
#'   \item \code{preds_ins}: Predictions computed using "in-sample" data used to train the model
#'   \item \code{preds_oos}: Forecasts computed using x data input to `fit_esn` but not used to train the model (i.e., "out-of-sample" data)
#'   \item \code{data_oos}: Data object created for the out-of-sample observations to obtain `preds_oos`
#'   \item \code{h_oos}: h matrix object created for the out-of-sample observations to obtain `preds_oos`
#' }
#' 
#' Objects output from `predict_esn` when new `x` values are input:
#' \itemize{
#'   \item \code{preds_new}: Forecasts computed using `x_new` data input to `predict_esn`
#'   \item \code{data_new}: Data object created for the new observations to obtain `preds_new`
#'   \item \code{h_new}: h matrix object created for the new observations to obtain `preds_new`
#' }
#' 
#' @export predict_esn
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
#'    )
#' 
#' # Get predictions on in-sample training data and forecasts on out-of-sample 
#' # training data
#' preds = predict_esn(model = esn)
#' 
#' # Get forecasts on new data
#' x_new = matrix(c(rnorm(5,10,1), rnorm(5,0,1)), ncol = 2, byrow = FALSE)
#' t_new = c("Jan", "Feb", "Mar", "Apr", "May")
#' t_new = paste0(t_new, "2022")
#' forecasts = predict_esn(model = esn, x_new = x_new, t_new = t_new)
#' 
#' # Get predictions on new (OOD) data
#' x_ood = matrix(c(rcauchy(12, 7, 1), rcauchy(12, 2, 1)), ncol = 2, byrow = FALSE)
#' new_pred = predict_esn(model = esn, x_ood = x_ood)

predict_esn <- function(model, x_new = NULL, t_new = NULL, x_ood = NULL, phi = NULL) {
  
  # Extract necessary values from model
  V = model$param_est$V
  h = model$h$h
  
  # Extract training data mean/sd for unscaling of predictions (if needed)
  if (model$internal_scaling == "joint") {
    y_train_scaled = model$data_train$y_train_scaled
    y_train_mean = model$data_train$y_train_mean
    y_train_sd = model$data_train$y_train_sd 
  }
  
  # Create out-of-sample data object
  data_obj_oos = create_data_obj_oos(model)
  
  # Obtain out-of-sample h matrix
  h_obj_oos <-
    create_h(
      U_and_W = model$param_est,
      data_obj = data_obj_oos,
      add_quad = model$add_quad,
      nh = model$params_tuning$nh,
      h_start = model$h$h_end
    )
  
  # TRAINING DATA
  # If no new x values are supplied, compute predictions on in-sample 
  # training data and forecasts on out-out-sample training data
  if (is.null(x_new) & is.null(t_new) & is.null(x_ood)) {
    
    # Compute scaled predictions (in-sample and out-of-sample)
    preds_ins = t(V %*% h)
    preds_oos = t(V %*% h_obj_oos$h)
    
    # Remove scaling (if previously applied)
    if (model$internal_scaling == "joint") {
      preds_ins = (preds_ins * y_train_sd) + y_train_mean
      preds_oos = (preds_oos * y_train_sd) + y_train_mean
    }
    
    # Add times to predictions
    rownames(preds_ins) = model$data_train$y_train_times
    x_oos_times = data_obj_oos$x_oos_times
    y_oos_times = paste(x_oos_times, "+", model$params_tuning$tau)
    rownames(preds_oos) = y_oos_times
    
    # Convert back to spatial scale is requested
    if (!is.null(phi)) {
      preds_ins = preds_ins %*% t(phi)
      preds_oos = preds_oos %*% t(phi)
    }
    
    # Put results in an object to be returns
    res <- 
      list(
        preds_ins = preds_ins, 
        preds_oos = preds_oos,
        data_oos = data_obj_oos,
        h_oos = h_obj_oos
      )
  
  # OUT-OF-DISTRIBUTION DATA
  # If new set of x inputs associated with the same times as the training data 
  # are supplied, compute predictions on OOD data
  } else if (!is.null(x_ood)){
    
    # Create a vector of times for OOD data (same times as training)
    t_ood <- model$data_input$t
    
    # Create new data input
    data_obj_ood <-
      create_data_obj_ood(
        x_ood = x_ood,
        t_ood = t_ood,
        model = model
      )
    
    # Obtain (out of distribution) h matrix
    h_obj_ood <-
      create_h(
        U_and_W = model$param_est,
        data_obj = data_obj_ood,
        add_quad = model$add_quad,
        nh = model$params_tuning$nh,
        h_start = model$h$h_start
      )
    
    # Compute scaled predictions
    preds_ood = t(V %*% h_obj_ood$h)
    
    # Remove scaling (if previously applied)
    if (model$internal_scaling == "joint") {
      preds_ood = (preds_ood * y_train_sd) + y_train_mean
    }
    
    # Add times to predictions
    rownames(preds_ood) = model$data_train$y_train_times
    
    # Convert back to spatial scale is requested
    if (!is.null(phi)) {
      preds_ood = preds_ood %*% t(phi)
    }
    
    # Put results in an object to be returns
    res <-
      list(
        preds_ood = preds_ood,
        data_ood = data_obj_ood,
        h_ood = h_obj_ood
      )
    
  # NEW DATA
  # If new x values at times occurring after the training data are supplied,
  # compute forecasts on new data
  } else {
    
    # Create new data object
    data_obj_new <-
      create_data_obj_new(
        x_new = x_new,
        t_new = t_new,
        model = model
      )
    
    # Obtain (out of sample) h matrix
    h_obj_new <-
      create_h(
        U_and_W = model$param_est,
        data_obj = data_obj_new,
        add_quad = model$add_quad,
        nh = model$params_tuning$nh,
        h_start = h_obj_oos$h_end
      )
    
    # Compute scaled predictions
    preds_new = t(V %*% h_obj_new$h)
    
    # Remove scaling (if previously applied)
    if (model$internal_scaling == "joint") {
      preds_new = (preds_new * y_train_sd) + y_train_mean
    }
    
    # Add times to predictions
    x_new_times = data_obj_new$x_new_times
    y_new_times = paste(x_new_times, "+", model$params_tuning$tau)
    rownames(preds_new) = y_new_times
    
    # Convert back to spatial scale is requested
    if (!is.null(phi)) {
      preds_new = preds_new %*% t(phi)
    }
    
    # Put results in an object to be returns
    res <-
      list(
        preds_new = preds_new,
        data_new = data_obj_new,
        h_new = h_obj_new
      )
    
  }
  
  # Return the predictions/forecasts
  return(res)
  
}
