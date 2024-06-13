#' Function for simulating simple synthetic data from a linear model with
#' spatio-temporal covariates, spatio-temporal noise, and white noise
#' 
#' Covariate X2 has real effect on response, covariate X1 has no effect on response
#'
#' @param nT number of times points
#' @param nx number of x values in grid
#' @param ny number of y values in grid
#' @param x1mean location of mode for x1
#' @param x2mean location of mode for x2
#' @param xsd standard deviation for covariates
#' @param sigmaX sqrt of sill of covariates (sqrt(covariance) at distance=0)
#' @param sigmaDelta sqrt of sill of spatio-temporal random effect (sqrt(covariance) at distance=0)
#' @param sigmaEpsilon sd of white noise
#' @param phiX range of covariates (sqrt(covariance) at distance=0)
#' @param phiDelta range of spatio-temporal random effect (sqrt(covariance) at distance=0)
#' @param rhoX AR(1) coefficient for temporal correlation of covariates
#' @param rhoDelta AR(1) coefficient for temporal correlation of spatio-temporal random effect
#' @param seed Random seed
#'
#' @export simulate_synthetic_data
#'
#' @importFrom dplyr left_join  
#' @importFrom MASS mvrnorm
#' @importFrom stats dist dnorm
#' @importFrom tidyr pivot_longer
#' 
#' @examples 
#' # Simulate data
#' dat = simulate_synthetic_data(nT = 70, nx = 10, ny = 10, x1mean = 20, 
#'                               x2mean = 45, xsd = 6, sigmaX = 2, 
#'                               sigmaDelta = 2, sigmaEpsilon = 2, 
#'                               phiX = .2, phiDelta = .2, rhoX = .2, 
#'                               rhoDelta = .2)

simulate_synthetic_data <- function(nT, nx, ny, x1mean, x2mean, xsd, sigmaX, 
                                    sigmaDelta, sigmaEpsilon, phiX, phiDelta, 
                                    rhoX, rhoDelta, seed = 1234) {
  
  # Set seed
  set.seed(seed)
  
  # Generate mean function of covariates
  mu1 <- 100 * dnorm(seq(from = 1, to = nT + 1, by = 1), x1mean, xsd)
  mu2 <- 100 * dnorm(seq(from = 1, to = nT + 1, by = 1), x2mean, xsd)
  
  # Create grid of values
  x.seq <- seq(from = 0, to = 1, length.out = nx)
  y.seq <- seq(from = 0, to = 1, length.out = ny)
  unit.grid <- expand.grid(easting = x.seq, northing = y.seq)
  nSp <- nrow(unit.grid)
  dist.mat <- as.matrix(dist(unit.grid))
  
  # Simulate X covariates
  Xt1 <- matrix(NA,nrow=nSp,ncol=nT+1)
  Xt1[,1] <- 
    MASS::mvrnorm(
      n = 1, 
      mu = rep(mu1[1], nSp),
      Sigma = sigmaX^2*exp(-(dist.mat^2)/(2*phiX^2))
    )
  Xt2 <- matrix(NA,nrow=nSp,ncol=nT+1)
  Xt2[,1] <- 
    MASS::mvrnorm(
      n = 1, 
      mu = rep(mu2[1], nSp),
      Sigma = sigmaX^2*exp(-(dist.mat^2)/(2*phiX^2))
    )
  for(i in 2:(nT+1)){
    Xt1[,i] <- 
      rhoX*Xt1[,i-1] + 
      MASS::mvrnorm(
        n = 1, 
        mu = rep(mu1[i], nSp),
        Sigma = sigmaX^2*exp(-(dist.mat^2)/(2*phiX^2))
      )
    Xt2[,i] <- 
      rhoX*Xt2[,i-1] + 
      MASS::mvrnorm(
        n = 1, 
        mu = rep(mu2[i], nSp),
        Sigma = sigmaX^2*exp(-(dist.mat^2)/(2*phiX^2))
      )
  }
  
  # Data frame with locations, times, and covariates
  Xt_df <- 
    data.frame(
      easting = rep(unit.grid$easting,nT+1),
      northing = rep(unit.grid$northing,nT+1),
      x1 = c(Xt1),
      x2 = c(Xt2),
      time = rep(1:(nT+1),each=(length(x.seq)*length(y.seq)))
    )
  
  # Simulate spatio-temporal random effect
  etaMat <- matrix(NA,nrow=nSp,ncol=nT+1)
  etaMat[,1] <-
    MASS::mvrnorm(
      n = 1,
      mu = rep(0, nSp),
      Sigma = sigmaDelta^2*exp(-(dist.mat^2)/(2*phiDelta^2))
    )
  for(i in 2:(nT+1)){
    etaMat[,i] <-
      rhoDelta*etaMat[,i-1] +
      MASS::mvrnorm(
        n = 1,
        mu = rep(0, nSp),
        Sigma = sigmaDelta^2*exp(-(dist.mat^2)/(2*phiDelta^2))
      )
  }
  
  # Put spatial random effects in data frame
  etaMat_df <-
    data.frame(
      easting = rep(unit.grid$easting, nT + 1),
      northing = rep(unit.grid$northing, nT + 1),
      delta = c(etaMat),
      time = rep(1:(nT + 1), each = (length(x.seq) * length(y.seq)))
    )
  
  # Specify relationship between x2 and response variable
  beta2 <- rep(1,nT+1)
  
  # Generate response 
  Y <- 
    t(sapply(1:(nT+1), function(x) Xt2[,x]*beta2[x])) + 
    matrix(rnorm(n = (nT+1)*nSp, sd = sigmaEpsilon), ncol = nSp) +
    t(etaMat)
  
  # Clean up data
  dat <- cbind(unit.grid, t(Y))
  names(dat)[-c(1:2)] <- as.character(1:(nT+1))
  dat <- 
    dat %>% 
    tidyr::pivot_longer(
      cols = !c(.data$easting, .data$northing),
      names_to = "time", 
      values_to = "val"
    )
  dat$time <- as.numeric(dat$time)
  dat <-
    left_join(dat, Xt_df, by = c("easting", "northing", "time")) %>%
    left_join(etaMat_df, by = c("easting", "northing", "time"))
  
  # Return data
  return(dat)
  
}
