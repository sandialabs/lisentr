
test_that("seed works", {
  
  # Create data
  x = matrix(rnorm(40,0,1), ncol = 4)
  y = matrix(rnorm(20,0,1), ncol = 2)
  colnames(x) = c("x11", "x12", "x21", "x22")
  colnames(y) = c("y1", "y2")
  t = paste0("t", 1:10)
  
  # Fit an ESN model
  esn <- 
    fit_esn(
      x = x, 
      y = y, 
      t = t, 
      tau = 2, 
      m = 1, 
      tau_emb = 1, 
      nh = 10, 
      seed = 1020349858
    )
  
  # Compute ZFI twice using same seed
  zfi1 <- 
    compute_fi(
      model = esn, 
      type = "zfi",
      var_groups = list(1:2, 3:4), 
      seed = 10208,
      blockSize = 1
    )
  zfi2 <- 
    compute_fi(
      model = esn, 
      type = "zfi",
      var_groups = list(1:2, 3:4), 
      seed = 10208,
      blockSize = 1
    )
  testthat::expect_identical(zfi1,zfi2)
  
  # Compute P FI twice using same seed
  pfi1 <- 
    compute_fi(
      model = esn, 
      type = "pfi",
      nreps = 1,
      var_groups = list(1:2, 3:4), 
      seed = 10203498,
      blockSize = 1
    )
  pfi2 <- 
    compute_fi(
      model = esn, 
      type = "pfi",
      nreps = 1,
      var_groups = list(1:2, 3:4), 
      seed = 10203498,
      blockSize = 1
    )
  testthat::expect_identical(pfi1,pfi2)
  
})
