#' Create embedding and response matrices
#'
#' This function creates a data object containing the matrix of embedding
#' vectors and the corresponding response variable matrix. Scaling is done
#' because Wikle et al. (2019) say "the responses and inputs are scaled by 
#' their respective standard deviations, as is common in the ESN literature".
#' 
#' @param x Training data inputs
#' @param x_ood Out of distribution data inputs (must be the same size matrix
#'        as training data x and occur at the same times as x)
#' @param x_new Data input matrix occurring after the times in the training
#'        data x
#' @param y Training data outputs
#' @param t Vector containing times/dates associated with x and y matrices
#' @param t_ood Vector of times associated with x_ood matrix
#' @param t_new Vector of times associated with x_new matrix
#' @param model ESN model output from fit_esn
#' @param tau Forecast lead time
#' @param m Embedding vector length
#' @param tau_emb Embedding lag
#' @param internal_scaling Specifies the type of standardization that should be 
#'        applied internally to x and y. This is common in ESN literature.
#'        Options are "joint" (the mean and standard deviation from all values 
#'        in the x/y matrix are computed after appropriate times are removed
#'        and used to center and scale the values in the x/y matrix) and "none"
#'        (no internal scaling is applied). Default is "joint".
#' 
#' @importFrom stats sd
#'
#' @name create_data_obj
#' @rdname create_data_obj
#' 
#' @examples 
#' # Create data
#' x = matrix(c(rnorm(12,10,1), rnorm(12,0,1)), ncol = 2, byrow = FALSE)
#' y = matrix(x[,1], ncol = 1)
#' 
#' # Assign column names
#' colnames(x) = c("X1", "X2")
#' colnames(y) = "X1"
#' 
#' # Create a vector of times
#' t = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
#'       "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
#' t = paste0(t, "2021")
#' 
#' # Tuning parameters
#' tau = 2
#' m = 1
#' tau_emb = 1
#' internal_scaling = TRUE
#' 
#' # Create data obj
#' data_obj <- 
#'   create_data_obj(
#'     x = x, 
#'     y = y, 
#'     t = t, 
#'     tau = tau, 
#'     m = m, 
#'     tau_emb = tau_emb, 
#'     internal_scaling = internal_scaling
#'  )
#'  
NULL

#' @rdname create_data_obj
#' @export
create_data_obj <- function (x, y, t, tau, m, tau_emb, internal_scaling) {
  
  # Check that x and y have the same number of rows
  if (dim(x)[1] != dim(y)[1]) 
    stop("x and y do not have the same amount of rows.")
  
  # Prepare x training data
  data_obj_x <-
    create_data_obj_x(
      x = x,
      t = t,
      tau = tau,
      m = m,
      tau_emb = tau_emb,
      internal_scaling = internal_scaling
    )
  
  # Prepare y training data
  data_obj_y <-
    create_data_obj_y(
      y = y,
      t = t,
      tau = tau,
      m = m,
      tau_emb = tau_emb,
      internal_scaling = internal_scaling
    )
  
  # Join two lists
  data_obj <- c(data_obj_y, data_obj_x)
  
  # Return values
  return(data_obj)
  
}

#' @rdname create_data_obj
#' @export
create_data_obj_x <- function (x, t, tau, m, tau_emb, internal_scaling) {
  
  #### Prepare x Training Data ####
  
  # Determine length of input data
  n_input_times = dim(x)[1]
  
  # Determine x time indices based on specified forecast lag (tau)
  x_train_index = 1:(n_input_times - tau)
  
  # Select x values to be used for training
  x_train = x[x_train_index, ]
  
  # Determine times associated with x training data
  x_train_times = t[x_train_index]
  
  # Save out of sample x values (not used for training for forecasts)
  x_oos = x[-x_train_index, ]
  
  # Determine times associated with sample x values
  x_oos_times = t[-x_train_index]
  
  # If specified, scale x values (common in the ESN literature)
  if (internal_scaling == "joint") {
    x_train_mean = mean(x_train)
    x_train_sd = sd(x_train)
    x_train_scaled = (x_train - x_train_mean) / x_train_sd
    x_emb_prep = x_train_scaled
  } else {
    x_train_scaled = NULL
    x_emb_prep = x_train
  }
  
  #### Create Embedding Vectors ####
  
  # Determine number of times in x used for training
  n_train_times = length(x_train_times)
  
  # Determine number of x variables (columns)
  n_xvars = ifelse(is.matrix(x_train), ncol(x_train), 1)
  
  # Determine number of embedding vectors -- this involves calculating and 
  # subtracting the number of observations that will be used for 
  # "initialization" (m*tau_emb) from the total times in the training data
  n_emb_vecs = n_train_times - (m * tau_emb)
  
  # Input the observations into the embedding vectors (and create a matching 
  # arraying containing the times corresponding to the embedding vectors)
  emb_vecs = array(NA, c(n_emb_vecs, m + 1, n_xvars))
  emb_vecs_times = emb_vecs
  for (i in 1:n_emb_vecs) {
    if (is.matrix(x_train)) {
      emb_vecs[i,,] = x_emb_prep[seq(i, (m * tau_emb + i), by = tau_emb),]
    } else {
      emb_vecs[i,,] = x_emb_prep[seq(i, (m * tau_emb + i), by = tau_emb)]
    }
    emb_vecs_times[i,,] = x_train_times[seq(i, (m * tau_emb + i), by = tau_emb)]
  }
  
  #### Create Design Matrix ####
  
  # Create the “design matrix”: number of rows correspond to the number of
  # embedding vectors (i.e., each row is an embedding vector) and the number
  # of columns correspond to one “intercept” column and the number of values in
  # an embedding vector (i.e., ((m+1)*n_xvars)) (Create a matching matrix
  # containing the times corresponding to the design matrix values)
  design_matrix = matrix(1, n_emb_vecs, (m + 1) * n_xvars + 1)
  design_matrix_times = matrix("Intercept", n_emb_vecs, (m + 1) * n_xvars + 1)
  for (i in 1:n_emb_vecs) {
    design_matrix[i, 2:((m + 1) * n_xvars + 1)] = as.vector(emb_vecs[i, ,])
    design_matrix_times[i, 2:((m + 1) * n_xvars + 1)] <-
      as.vector(emb_vecs_times[i, ,])
  }
  
  #### Output ####
  
  # Return values
  if (internal_scaling == "joint") {
    res <- list(
      x_train = x_train,
      x_train_times = x_train_times,
      x_train_scaled = x_train_scaled,
      x_train_mean = x_train_mean,
      x_train_sd = x_train_sd,
      x_oos = x_oos,
      x_oos_times = x_oos_times,
      emb_vecs = emb_vecs,
      emb_vecs_times = emb_vecs_times,
      n_emb_vecs = n_emb_vecs,
      design_matrix = design_matrix,
      design_matrix_times = design_matrix_times
    )
  } else {
    res <- list(
      x_train = x_train,
      x_train_times = x_train_times,
      x_train_scaled = x_train_scaled,
      x_oos = x_oos,
      x_oos_times = x_oos_times,
      emb_vecs = emb_vecs,
      emb_vecs_times = emb_vecs_times,
      n_emb_vecs = n_emb_vecs,
      design_matrix = design_matrix,
      design_matrix_times = design_matrix_times
    )
  }
  return(res)
  
}

#' @rdname create_data_obj
#' @export
create_data_obj_y <- function (y, t, tau, m, tau_emb, internal_scaling) {
  
  #### Create Response Matrix ####
  
  # Determine length of input data
  n_input_times = dim(y)[1]
  
  # Determine y time indices based on specified forecast lag (tau), embedding 
  # vector lag (tau_emb), and embedding vector length (m + 1)
  y_train_index = (tau + m * tau_emb + 1):(n_input_times)
  
  # Select y values to be used for training
  y_train = y[y_train_index,]
  
  # Determine times associated with y training data
  y_train_times = t[y_train_index]
  
  # If specified, scale y values (common in the ESN literature)
  if (internal_scaling == "joint") {
    y_train_mean = mean(y_train)
    y_train_sd = sd(y_train)
    y_train_scaled = (y_train - y_train_mean) / y_train_sd
  } else {
    y_train_scaled = NULL
  }
  
  #### Output ####
  
  # Return values
  if (internal_scaling == "joint") {
    list(
      y_train = y_train,
      y_train_times = y_train_times,
      y_train_scaled = y_train_scaled,
      y_train_mean = y_train_mean,
      y_train_sd = y_train_sd
    )
  } else {
    return(
      list(
        y_train = y_train,
        y_train_times = y_train_times,
        y_train_scaled = y_train_scaled
      )
    )
  }
  
}

#' @rdname create_data_obj
#' @export
create_data_obj_oos <- function (model) {
  
  # Extract params from model
  tau = model$params_tuning$tau
  tau_emb = model$params_tuning$tau_emb
  m = model$params_tuning$m
  
  # Extract training data x and times
  x = model$data_input$x
  t = model$data_input$t
  
  # Determine number of times and variables in input data
  n_input_times = length(t)
  n_xvars = dim(x)[2]
  
  # Determine indices of x to be included in out of sample data
  x_oos_index = (n_input_times - tau + 1):n_input_times
  
  # Determine number of out of sample embedding vectors
  n_emb_vecs_oos = length(x_oos_index)
  
  # If specified, scale x values (common in the ESN literature)
  if (model$internal_scaling == "joint") {
    x_scaled = (x - model$data_train$x_train_mean) / model$data_train$x_train_sd
    x_oos_scaled = x_scaled[x_oos_index,]
    x_emb_prep = x_scaled
  } else {
    x_oos_scaled = NULL
    x_emb_prep = x
  }
  
  # Create out of sample embedding vectors
  emb_vecs_oos = array(NA, c(n_emb_vecs_oos, m + 1, n_xvars))
  emb_vecs_oos_times = emb_vecs_oos
  for (i in 1:n_emb_vecs_oos) {
    if (is.matrix(x_emb_prep)) {
      emb_vecs_oos[i,,] <- 
        x_emb_prep[seq(x_oos_index[i] - 
                         (m * tau_emb), x_oos_index[i], by = tau_emb),]
    } else {
      emb_vecs_oos[i,,] <-
        x_emb_prep[seq(x_oos_index[i] - 
                         (m * tau_emb), x_oos_index[i], by = tau_emb)]
    }
    emb_vecs_oos_times[i,,] <- 
      t[seq(x_oos_index[i] - (m * tau_emb), x_oos_index[i], by = tau_emb)]
  }
  
  # Create an out of sample design matrix
  design_matrix_oos = matrix(1, n_emb_vecs_oos, (m + 1) * n_xvars + 1)
  design_matrix_oos_times <-
    matrix("Intercept", n_emb_vecs_oos, (m + 1) * n_xvars + 1)
  for (i in 1:n_emb_vecs_oos) {
    design_matrix_oos[i, 2:((m + 1) * n_xvars + 1)] <-
      as.vector(emb_vecs_oos[i,,])
    design_matrix_oos_times[i, 2:((m + 1) * n_xvars + 1)] <- 
      as.vector(emb_vecs_oos_times[i, ,])
  }
  
  # Return values
  return(
    list(
      x_oos = model$data_train$x_oos,
      x_oos_times = model$data_train$x_oos_times,
      x_oos_scaled = x_oos_scaled,
      emb_vecs = emb_vecs_oos,
      emb_vecs_times = emb_vecs_oos_times,
      n_emb_vecs = dim(design_matrix_oos)[1],
      design_matrix = design_matrix_oos,
      design_matrix_times = design_matrix_oos_times
    )
  )
  
}

#' @rdname create_data_obj
#' @export
create_data_obj_ood <- function (model, x_ood, t_ood) {
  
  #### Adjust OOD Data Based on Forecast Lag and Scaling Requested ####
  
  # Determine length of input data
  n_input_times = length(t_ood)
  
  # Check to make sure number of OOD times matches training data
  if (n_input_times != dim(x_ood)[1])
    stop("Number of rows in x must match the number of times in the
         model input data")
  
  # Determine x time indices based on specified forecast lag (tau)
  x_ood_index = 1:(n_input_times - model$params_tuning$tau)
  
  # Extract observations from OOD data that occur after training data
  x_ood_oos = x_ood[-x_ood_index, ]
  
  # Determine times associated with OOS OOD observations
  x_ood_oos_times = t_ood[-x_ood_index]
  
  # Subset OOD data to only include times that match the training data
  x_ood = x_ood[x_ood_index, ]
  
  # Determine times associated with x training data
  x_ood_times = t_ood[x_ood_index]
  
  # If specified, scale x values (common in the ESN literature)
  if (model$internal_scaling == "joint") {
    x_train_mean = model$data_train$x_train_mean
    x_train_sd = model$data_train$x_train_sd
    x_ood_scaled = (x_ood - x_train_mean) / x_train_sd
    x_emb_prep = x_ood_scaled
  } else {
    x_ood_scaled = NULL
    x_emb_prep = x_ood
  }
  
  #### Create Embedding Vectors ####
  
  # Determine number of times in x_ood
  n_ood_times = length(x_ood_times)
  
  # Determine number of x variables (columns)
  n_xvars = ifelse(is.matrix(x_ood), ncol(x_ood), 1)
  
  # Extract model parameters
  m = model$params_tuning$m
  tau_emb = model$params_tuning$tau_emb
  
  # Determine number of embedding vectors -- this involves calculating and
  # subtracting the number of observations that will be used for
  # "initialization" (m*tau_emb) from the total times in the training data s
  n_emb_vecs = n_ood_times - (m * tau_emb)
  
  # Input the observations into the embedding vectors (and create a matching
  # arraying containing the times corresponding to the embedding vectors)
  emb_vecs = array(NA, c(n_emb_vecs, model$params_tuning$m + 1, n_xvars))
  emb_vecs_times = emb_vecs
  for (i in 1:n_emb_vecs) {
    if (is.matrix(x_ood)) {
      emb_vecs[i,,] = x_emb_prep[seq(i, (m * tau_emb + i), by = tau_emb),]
    } else {
      emb_vecs[i,,] = x_emb_prep[seq(i, (m * tau_emb + i), by = tau_emb)]
    }
    emb_vecs_times[i,,] = x_ood_times[seq(i, (m * tau_emb + i), by = tau_emb)]
  }
  
  #### Create Design Matrix ####
  
  # Create the “design matrix”: number of rows correspond to the number of
  # embedding vectors (i.e., each row is an embedding vector) and the number
  # of columns correspond to one “intercept” column and the number of values in
  # an embedding vector (i.e., ((m+1)*n_xvars)) (Create a matching matrix
  # containing the times corresponding to the design matrix values)
  design_matrix = matrix(1, n_emb_vecs, (m + 1) * n_xvars + 1)
  design_matrix_times = matrix("Intercept", n_emb_vecs, (m + 1) * n_xvars + 1)
  for (i in 1:n_emb_vecs) {
    design_matrix[i, 2:((m + 1) * n_xvars + 1)] = as.vector(emb_vecs[i, ,])
    design_matrix_times[i, 2:((m + 1) * n_xvars + 1)] <-
      as.vector(emb_vecs_times[i, ,])
  }
  
  #### Output ####
  
  # Return values
  res <- list(
    x_ood = x_ood,
    x_ood_times = x_ood_times,
    x_ood_scaled = x_ood_scaled,
    x_ood_oos = x_ood_oos,
    x_ood_oos_times = x_ood_oos_times,
    emb_vecs = emb_vecs,
    emb_vecs_times = emb_vecs_times,
    n_emb_vecs = n_emb_vecs,
    design_matrix = design_matrix,
    design_matrix_times = design_matrix_times
  )
  return(res)
  
}

#' @rdname create_data_obj
#' @export
create_data_obj_new <- function (model, x_new, t_new) {
  
  # Extract params from model
  tau = model$params_tuning$tau
  tau_emb = model$params_tuning$tau_emb
  m = model$params_tuning$m
  
  # Extract training data x and times
  x = model$data_input$x
  t = model$data_input$t
  
  # Join training and new x values and times
  x_all = rbind(x, x_new)
  t_all = c(t, t_new)
  
  # Determine number of times and variables in input data
  n_input_times = length(t)
  n_xvars = dim(x)[2]
  
  # Determine length of new data
  n_new_times = length(t_new)
  
  # Determine indices of x and y values to be included in new data
  x_new_index = seq(n_input_times + 1, by = 1, length.out = n_new_times)
  
  # Determine number of out of sample embedding vectors
  n_emb_vecs_new = length(x_new_index)
  
  # If specified, scale x values (common in the ESN literature)
  if (model$internal_scaling == "joint") {
    x_all_scaled <-
      (x_all - model$data_train$x_train_mean) / model$data_train$x_train_sd
    x_new_scaled = x_all_scaled[x_new_index, ]
    x_emb_prep = x_all_scaled
  } else {
    x_new_scaled = NULL
    x_emb_prep = x_all
  }
  
  # Create embedding vectors based on new x values
  emb_vecs_new = array(NA, c(length(x_new_index), m + 1, n_xvars))
  emb_vecs_new_times = emb_vecs_new
  for (i in 1:length(x_new_index)) {
    if (is.matrix(x_emb_prep)) {
      emb_vecs_new[i, , ] <- 
        x_emb_prep[seq(x_new_index[i] - 
                         (m * tau_emb), x_new_index[i], by = tau_emb), ]
    } else {
      emb_vecs_new[i, , ] <-
        x_emb_prep[seq(x_new_index[i] - 
                         (m * tau_emb), x_new_index[i], by = tau_emb)]
    }
    emb_vecs_new_times[i, , ] <- 
      t_all[seq(x_new_index[i] - (m * tau_emb), x_new_index[i], by = tau_emb)]
  }
  
  # Create an out of sample design matrix
  design_matrix_new = matrix(1, length(x_new_index), (m + 1) * n_xvars + 1)
  design_matrix_new_times <-
    matrix("Intercept", length(x_new_index), (m + 1) * n_xvars + 1)
  for (i in 1:length(x_new_index)) {
    design_matrix_new[i, 2:((m + 1) * n_xvars + 1)] <- 
      as.vector(emb_vecs_new[i, , ])
    design_matrix_new_times[i, 2:((m + 1) * n_xvars + 1)] <-
      as.vector(emb_vecs_new_times[i, , ])
  }
  
  # Return values associated with embedding vector
  return(
    list(
      x_new = x_new,
      x_new_times = t_new,
      x_new_scaled = x_new_scaled,
      emb_vecs = emb_vecs_new,
      emb_vecs_times = emb_vecs_new_times,
      n_emb_vecs = n_emb_vecs_new,
      design_matrix = design_matrix_new,
      design_matrix_times = design_matrix_new_times
    )
  )
  
}
