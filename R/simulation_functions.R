#' @importFrom data.table data.table :=
#' @importFrom stats predict quantile sigma
#' @importFrom utils data flush.console
NULL

#' @title Simulate Poisson Covariate Data
#' @description Generates covariate data for Poisson simulations.
#'
#' @param n Number of observations to simulate.
#' @param index Optional index to include in the dataset.
#'
#' @return A `data.table` with columns `x1`, `x2` and an index column if not provided.
#' @keywords internal
#' @noRd
sim_poisson_x <- function(n, index = NULL){
  data <- data.table(x1 = x1, x2 = x2)
  if( is.null(index) ) {
    data[,fit_data_row:=.I]
    index <- "fit_data_row"
  } else {
    data <- cbind(index,data)
  }
  return(data)
}

#' @title Get Poisson \eqn{\lambda(x_{i})}
#' @description Computes the lambda parameter for Poisson regression using linear predictors.
#'
#' @param x A `data.table` with covariates `x1`, `x2`.
#' @param betas Optional vector of regression coefficients. Defaults to `c(1,3,0)`.
#'
#' @return A numeric vector of Poisson means (lambda values).
#' @keywords internal
#' @noRd
get_lambda_poisson <- function(x, betas = NULL){
  x <- x[,.SD,.SDcols = c("x1","x2")]
  x <- cbind(1,x)
  if(is.null(betas)){
    betas <- c(1,3,0)
  }
  lambda <- exp(tcrossprod(as.matrix(x),matrix(betas,nrow=1)))
  return(c(lambda))
}


#' Generate Data from a Two-Component Discrete Mixture with Dynamic Tail
#'
#' Simulates a count outcome based on a discrete mixture model. The lower tail is modeled via a discrete gamma distribution with covariate-dependent scale, and the upper tail is complemented by a discrete Generalized Pareto distribution (DGP), adjusted to match a given exceedance probability threshold.
#'
#' @param n Integer. Number of observations to simulate.
#' @param phi Numeric. Exceedance probability threshold (e.g., 0.05 defines the top 5% tail).
#' @param xi Numeric. Shape parameter of the discrete Generalized Pareto (DGP) distribution.
#' @param scale_dgp Numeric. Scale parameter for the DGP distribution.
#' @param shape_dgamma Numeric. Shape parameter for the discrete gamma distribution.
#' @param beta_scale_dgamma Numeric vector. Coefficients for modeling the scale parameter of the gamma distribution as a function of covariates.
#' @param probs Numeric vector. Quantile levels at which to compute the true quantiles of the generated distribution.
#' @param index Optional. If not NULL, an index variable will be added for tracking observations.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{data}{A \code{data.table} with covariates and generated outcome \code{y}.}
#'   \item{pmf_dggp}{A matrix containing the probability mass function (PMF) values for the mixture.}
#'   \item{true_quantiles}{A matrix of true quantiles for each observation at the levels specified in \code{probs}.}
#' }
#'
#' @details This function dynamically determines a cut-off threshold from the discrete gamma distribution such that the upper tail (defined by \code{phi}) is filled using a DGP distribution with adjusted mass. Covariate effects are modeled through the scale parameter using a linear predictor on the log scale.
#'
#' @examples
#' out <- mixture_dgamma_constant_tail()
#' head(out$data)
#'
#' @export
mixture_dgamma_constant_tail <- function(n = 100,
                                        phi = 0.05,
                                        xi = 0.3,
                                        scale_dgp = 15,
                                        shape_dgamma = 5,
                                        beta_scale_dgamma = c(1,1),
                                        probs = c(0.001,0.01,0.05,0.25,0.5,0.75,0.8,0.9,0.95,0.99,0.999,0.9999),
                                        index = NULL){

  scale_x <- function(x,a,b){x*(b-a)+a}

  x1 <- scale_x(x = stats::rbeta(n = n,shape1 = 5/3,shape2 = 5/3),
                a = -5/3,
                b = 5/3)

  x <- data.table(x1 = x1)
  x <- cbind(1,x)
  betas_scale_dgamma <- c(1,1)

  # Using the scale as a function of x;
  scale_parameter_dgamma <- exp(tcrossprod(as.matrix(x),matrix(betas_scale_dgamma,nrow=1)))

  # Using a fixed scale
  # scale_parameter_dgamma <- rep(scale_dgamma,n)

  sample <- extraDistr::rdgamma(n = n,shape = 1,scale = scale_parameter_dgamma)

  # Getting the 95% quantile
  dgamma_cdf <- do.call(rbind,lapply(scale_parameter_dgamma,function(x){extraDistr::pdgamma(q = 0:1000,shape = shape_dgamma,scale = x)}))
  dgamma_pmf <- do.call(rbind,lapply(scale_parameter_dgamma,function(x){extraDistr::ddgamma(x = 0:1000,shape = shape_dgamma,scale = x)}))

  index_quantile <- apply(dgamma_cdf,1, function(y){which.min(abs(y[(y-(1-phi))<0]-(1-phi)))})
  update_phi <- numeric()
  for(i in 1:nrow(dgamma_cdf)){
    update_phi[i] <- (1-dgamma_cdf[i,index_quantile[i]])
  }

  # Getting the probabilities for the mixture
  xx_max <- ncol(dgamma_cdf)-1
  cdf_dggp <- pmf_dggp <- matrix(0,nrow = nrow(x),ncol = ncol(dgamma_cdf))
  colnames(cdf_dggp) <- colnames(pmf_dggp) <- as.character(0:xx_max)
  true_quantiles <- matrix(0, nrow = nrow(x), ncol = length(probs))

  y <- numeric(nrow(x))

  for(j in 1:nrow(x)){
    pmf_dggp[j,1:(index_quantile[j])] <- dgamma_pmf[j,1:(index_quantile[j])]
    pmf_dggp[j,(index_quantile[j]+1):(xx_max+1)] <- update_phi[j]*ddgp(x = 0:(xx_max-index_quantile[j]),shape = xi,scale = scale_dgp)
    cdf_dggp[j,] <- cumsum(pmf_dggp[j,])
    true_quantiles[j, ] <- sapply(probs, function(p){which((cdf_dggp[j,]-p)>0)[1]-2})
    y[j] <- sample(x=0:xx_max,size = 1, prob = pmf_dggp[j,],replace = TRUE)
  }

  true_quantiles <- ifelse(true_quantiles<0,0,true_quantiles)
  colnames(true_quantiles) <- probs


  # To check the plot do
  plot(table(y), xlab = "Observed values", ylab = "Frequency")
  return(list( data = data.table(x[,.SD,.SDcols = c("x1")], y = y),
               pmf_dggp = pmf_dggp,
               true_quantiles = true_quantiles))
}


#' Generate Data from a Two-Component Discrete Mixture with Fixed Cut-off
#'
#' Simulates data from a discrete mixture where the lower tail is from a discrete gamma distribution and the upper tail is filled by a discrete Generalized Pareto (DGP) beyond a fixed quantile threshold.
#'
#' @param n Integer. Number of observations to simulate.
#' @param phi Numeric. Tail probability beyond which the DGP distribution is used.
#' @param xi Numeric. Shape parameter for the DGP distribution.
#' @param scale_dgp Numeric. Scale parameter for the DGP distribution.
#' @param shape_dgamma Numeric. Shape parameter for the discrete gamma distribution.
#' @param scale_dgamma Numeric. Fixed scale parameter for the gamma distribution.
#' @param probs Numeric vector. Quantile levels for which to compute the distribution quantiles.
#' @param index Optional. If not NULL, an index variable will be added to the output.
#'
#' @return A list with:
#' \describe{
#'   \item{data}{A \code{data.table} with covariates and generated outcome \code{y}.}
#'   \item{pmf_dggp}{Named vector containing the PMF of the full mixture.}
#'   \item{true_quantiles}{Vector of quantiles corresponding to \code{probs}.}
#' }
#'
#' @details This function defines a fixed upper-tail threshold using the gamma quantile and then reweights the DGP to fill the tail such that the mixture retains a valid PMF.
#'
#' @examples
#' out <- mixture_dgamma_constant()
#' head(out$data)
#'
#' @export

mixture_dgamma_constant <- function(n = 100,
                                    phi = 0.05,
                                    xi = 0.3,
                                    scale_dgp = 15,
                                    shape_dgamma = 5,
                                    scale_dgamma = 4,
                                    probs = c(0.001,0.01,0.05,0.25,0.5,0.75,0.8,0.9,
                                              0.95,0.99,0.999,0.9999,0.99999),
                                    index = NULL){

  scale_x <- function(x,a,b){x*(b-a)+a}

  x1 <- scale_x(x = stats::rbeta(n = n,shape1 = 5/3,shape2 = 5/3),
                a = -5/3,
                b = 5/3)

  # x1 <- matrix(1,nrow=n,ncol=1)

  x <- data.table(x1 = x1)
  colnames(x) <- "x1"
  betas_scale_dgamma <- c(1)
  betas_scale_dgp <- c(1)


  # Using the scale as a function of x;
  scale_parameter_dgamma <- exp(tcrossprod(as.matrix(x),matrix(betas_scale_dgamma,nrow=1)))

  # Using a fixed scale

  sample <- extraDistr::rdgamma(n = n,shape = 1,scale = scale_parameter_dgamma)
  extraDistr::pdgamma(q = 0:20,shape = 1,scale = 5)
  # Getting the 95% quantile
  u <- quantile(extraDistr::rdgamma(n = 10000000,shape = shape_dgamma,scale = scale_dgamma),prob = (1-phi))
  update_phi <- extraDistr::pdgamma(q = u-1,shape = shape_dgamma,scale = scale_dgamma,lower.tail = FALSE)

  # Getting the probabilities for the mixture
  if(xi<0){
    xx_max <- -scale_dgp/xi
  } else {
     xx_max <- 1000
  }

  pmf_dggp <- numeric(xx_max+1)
  names(pmf_dggp) <- 0:xx_max
  pmf_dggp[as.character(0:(u-1))] <- extraDistr::ddgamma(x = 0:(u-1),shape = shape_dgamma,scale = scale_dgamma)
  pmf_dggp[as.character(u:xx_max)] <-  pmf_dggp[as.character(u:xx_max)] + update_phi*ddgp(x = 0:(length(u:xx_max)-1),shape = xi,scale = scale_dgp)

  cdf_dggp <- cumsum(pmf_dggp)
  true_quantiles <- sapply(probs, function(p){as.numeric(names(which((cdf_dggp-p)>0)[1]))-1})
  true_quantiles[true_quantiles<0] <- 0

  y <- sample(x=0:xx_max,size = n, prob = pmf_dggp,replace = TRUE)


  if(!is.null(index)){
    x[,fit_data_row:=1:.N]
  }

  # To check the plot do
  plot(table(y), xlab = "Observed values", ylab = "Frequency")
  return(list( data = data.table(x, y = y),
               pmf_dggp = pmf_dggp,
               true_quantiles = true_quantiles))
}

#' Simulate Data from a Conditional Discrete Mixture with Gamma and Generalized Pareto
#'
#' Simulates a continuous-like discrete outcome by generating values from a discrete gamma distribution and applying DGP shifts to values exceeding a covariate-dependent threshold.
#'
#' @param n Integer. Number of observations.
#' @param alpha_right Numeric. Right-tail probability used to define the exceedance threshold per observation.
#' @param xi Numeric. Shape parameter of the DGP distribution.
#'
#' @return A list with:
#' \describe{
#'   \item{data}{A \code{data.table} with covariates \code{x1}, \code{x2}, and simulated outcome \code{y}.}
#'   \item{exceedances}{Indices of the observations exceeding the gamma quantile threshold.}
#'   \item{scale_parameter_dgp}{Vector of scale parameters used in the DGP component.}
#' }
#'
#' @details For observations with values exceeding a certain quantile of the gamma distribution, an extra count from the DGP is added. This structure allows modeling heavy-tailed count data with heteroscedastic features.
#'
#' @examples
#' out <- mixture_dgamma_dgp()
#' summary(out$data$y)
#'
#' @export

mixture_dgamma_dgp <- function(n = 100,
                               alpha_right = 0.8,
                               xi = 0.3){


    scale_x <- function(x,a,b){x*(b-a)+a}

    x1 <- scale_x(x = stats::rbeta(n = n,shape1 = 5/3,shape2 = 5/3),
                  a = -5/3,
                  b = 5/3)
    x2 <- scale_x(x = stats::rbeta(n = n,shape1 = 5/3,shape2 = 5/3),
                  a = -5/3,
                  b = 5/3)

    x <- data.table(x1 = x1, x2 = x2)
    x <- cbind(1,x)
    betas_scale_dgamma <- c(1,0.5,0)
    betas_scale_dgp <- c(1,0,2)

    # Using the scale as a function of x;
    scale_parameter_dgamma <- exp(tcrossprod(as.matrix(x),matrix(betas_scale_dgamma,nrow=1)))

    # Using a fixed scale
    scale_parameter_dgamma <- rep(1,n)

    sample <- extraDistr::rdgamma(n = n,shape = 1,scale = scale_parameter_dgamma)
    q_threshold_dgamma <- sapply(scale_parameter_dgamma,function(y){stats::quantile(extraDistr::rdgamma(n = 100000,shape = 1,scale = y),prob = alpha_right)})
    exceedances <- which(sample>=q_threshold_dgamma)
    scale_parameter_dgp <- exp(tcrossprod(as.matrix(x[exceedances,]),matrix(betas_scale_dgp,nrow=1)))

    scale_exceed <- sapply(scale_parameter_dgp,function(y){sample(x = 0:500, size = 1,prob = ddgp(x = 0:500,shape = xi,scale = y))})

    sample[exceedances] <- sample[exceedances] + scale_exceed

    return(list( data = data.table(x[,.SD,.SDcols = c("x1","x2")], y = sample),
                 exceedances = exceedances,
                 scale_parameter_dgp = scale_parameter_dgp))
}



#' Simulate Data from a Discrete Generalized Pareto Process
#'
#' Generates data based solely on a DGP model where the scale parameter depends on covariates. Useful as a pure heavy-tail generator.
#'
#' @param n Integer. Number of samples to generate.
#' @param xi Numeric. Shape parameter for the DGP distribution.
#'
#' @return A list with:
#' \describe{
#'   \item{data}{A \code{data.table} with covariates \code{x1}, \code{x2}, and simulated outcome \code{y}.}
#'   \item{scale_parameter_dgp}{Vector of scale parameters for the DGP model per observation.}
#' }
#'
#' @details The covariates \code{x1} and \code{x2} are generated from scaled beta distributions and used to model the scale of the DGP distribution via a log-linear formulation.
#'
#' @examples
#' out <- sim_dgp()
#' plot(table(out$data$y), xlab = "y", ylab = "Frequency")
#'
#' @export

sim_dgp <- function(n = 100,
                    xi =  0.3){

  scale_x <- function(x,a,b){x*(b-a)+a}

  x1 <- scale_x(x = stats::rbeta(n = n,shape1 = 5/3,shape2 = 5/3),
                a = -5/3,
                b = 5/3)
  x2 <- scale_x(x = stats::rbeta(n = n,shape1 = 5/3,shape2 = 5/3),
                a = -5/3,
                b = 5/3)

  x <- data.table(x1 = x1, x2 = x2)
  x <- cbind(1,x)
  betas_scale_dgamma <- c(1,0.5,0)
  betas_scale_dgp <- c(1,0,1)

  scale_parameter_dgp <- exp(tcrossprod(as.matrix(x),matrix(betas_scale_dgp,nrow=1)))

  sample_ <- sapply(scale_parameter_dgp,function(y){sample(x = 0:500, size = 1,prob = ddgp(x = 0:500,shape = xi,scale = y))})


  return(list( data = data.table(x[,.SD,.SDcols = c("x1","x2")], y = sample_),
               scale_parameter_dgp = scale_parameter_dgp))
}


#' Simulate Data for Scenario One (Scale-Invariant Tail)
#'
#' Generates synthetic count data where the bulk is modeled by a discrete gamma distribution
#' and the tail by a discrete generalized Pareto (DGP) with constant scale.
#'
#' @param n Integer. Number of samples to generate.
#' @param phi Numeric. Exceedance probability.
#' @param xi Numeric. Shape parameter of the DGP tail.
#' @param scale_dgp Numeric. Scale of the DGP tail.
#' @param shape_dgamma Numeric. Shape of the discrete gamma bulk distribution.
#' @param weibull_shape Numeric. Unused. Kept for backwards compatibility.
#' @param weibull_scale Numeric. Unused. Kept for backwards compatibility.
#' @param beta_scale_dgamma Numeric vector. Coefficients for the gamma scale regression.
#' @param probs Numeric vector. Probability levels to compute true quantiles.
#' @param index Optional index to be returned for tracking.
#'
#' @return A list with:
#'   \item{data}{Data.table with covariates and generated response.}
#'   \item{pmf_dggp}{Matrix of the simulated probability mass function.}
#'   \item{true_quantiles}{Matrix of true quantiles at specified probabilities.}
#'
#' @details
#' This simulation corresponds to Section 3.1.1 of the main manuscript for the xflex4cast model.
#'
#' @export
simulation_scenario_one <- function(n = 5000,
                         phi = 0.1,
                         xi = 0.3,
                         scale_dgp = 2.5,
                         shape_dgamma = 1.5,
                         weibull_shape = 2.76,
                         weibull_scale = 5.90,
                         beta_scale_dgamma = c(1,1),
                         probs = c(0.001,0.01,0.05,0.25,0.5,0.75,0.8,0.9,0.95,0.99,0.999,0.9999),
                         index = NULL){

  scale_to_range <- function(x, a, b) {
    (x - min(x)) / (max(x) - min(x)) * (b - a) + a
  }

  sim_data <- x1_x2_df_simulation
  sim_data_sample <- sample(1:nrow(sim_data), size = n, replace = TRUE)
  x1 <- scale_to_range(sim_data$x1[sim_data_sample], a = -0.5, b = 5/3)
  x2 <- sim_data$x2[sim_data_sample]

  x <- data.table::data.table(x1 = x1, x2 = x2)
  x <- cbind(1, x)
  betas_scale_dgamma <- c(1, 2, 0)
  scale_parameter_dgamma <- exp(tcrossprod(as.matrix(x), matrix(betas_scale_dgamma, nrow = 1)))

  sample <- extraDistr::rdgamma(n = n, shape = shape_dgamma, scale = scale_parameter_dgamma)

  dgamma_cdf <- do.call(rbind, lapply(scale_parameter_dgamma, function(x) {
    extraDistr::pdgamma(q = 0:1000, shape = shape_dgamma, scale = x)
  }))
  dgamma_pmf <- do.call(rbind, lapply(scale_parameter_dgamma, function(x) {
    extraDistr::ddgamma(x = 0:1000, shape = shape_dgamma, scale = x)
  }))

  index_quantile <- apply(dgamma_cdf, 1, function(y) {
    which.min(abs(y[(y - (1 - phi)) < 0] - (1 - phi)))
  })

  update_phi <- numeric()
  for (i in 1:nrow(dgamma_cdf)) {
    update_phi[i] <- 1 - dgamma_cdf[i, index_quantile[[i]]]
  }

  xx_max <- ncol(dgamma_cdf) - 1
  cdf_dggp <- pmf_dggp <- matrix(0, nrow = nrow(x), ncol = ncol(dgamma_cdf))
  colnames(cdf_dggp) <- colnames(pmf_dggp) <- as.character(0:xx_max)
  true_quantiles <- matrix(0, nrow = nrow(x), ncol = length(probs))
  y <- numeric(nrow(x))

  for (j in 1:nrow(x)) {
    pmf_dggp[j, 1:(index_quantile[j])] <- dgamma_pmf[j, 1:(index_quantile[j])]
    if (index_quantile[j] < (xx_max + 1)) {
      pmf_dggp[j, (index_quantile[j] + 1):(xx_max + 1)] <- update_phi[j] * ddgp(x = 0:(xx_max - index_quantile[j]), shape = xi, scale = scale_dgp)
    } else {
      pmf_dggp[j, index_quantile[j]] <- update_phi[j] * ddgp(x = 0, shape = xi, scale = scale_dgp[j])
    }
    cdf_dggp[j, ] <- cumsum(pmf_dggp[j, ])
    true_quantiles[j, ] <- sapply(probs, function(p) which((cdf_dggp[j, ] - p) > 0)[1] - 2)
    y[j] <- sample(x = 0:xx_max, size = 1, prob = pmf_dggp[j, ], replace = TRUE)
  }

  true_quantiles <- ifelse(true_quantiles < 0, 0, true_quantiles)
  colnames(true_quantiles) <- probs


  true_quantiles <- ifelse(true_quantiles < 0, 0, true_quantiles)
  colnames(true_quantiles) <- probs

  list(
    data = data.table::data.table(fit_data_row = 1:n, x1 = x$x1, x2 = x$x2, y = y),
    pmf_dggp = pmf_dggp,
    true_quantiles = true_quantiles
  )
}


#' Simulate Data for Scenario Two (Non-Stationary Tail)
#'
#' Extends Scenario One by allowing the DGP scale parameter
#' to vary with an additional predictor.
#'
#' @inheritParams simulation_scenario_one
#' @return Same structure as \code{scenario_one}.
#' @details
#' This simulation corresponds to Section 3.1.1 of the main manuscript for
#' the xflex4cast model.
#'
#' @export
simulation_scenario_two <- function(n = 5000,
                         phi = 0.1,
                         xi = 0.3,
                         shape_dgamma = 1.5,
                         weibull_shape = 2.76,
                         weibull_scale = 5.90,
                         beta_scale_dgamma = c(1,1),
                         probs = c(0.001,0.01,0.05,0.25,0.5,0.75,0.8,0.9,0.95,0.99,0.999,0.9999),
                         index = NULL){

  scale_to_range <- function(x, a, b) {
    (x - min(x)) / (max(x) - min(x)) * (b - a) + a
  }

  sim_data <- x1_x2_df_simulation # This was calculated using the real data, through empirical CDF
  sim_data_sample <- sample(1:nrow(sim_data), size = n, replace = TRUE)

  x1 <- scale_to_range(sim_data$x1[sim_data_sample], a = -0.5, b = 5/3)
  x2 <- sim_data$x2[sim_data_sample]

  x <- data.table::data.table(x1 = x1, x2 = x2)
  x <- cbind(1, x)

  betas_scale_dgamma <- c(1, 2, 0)
  betas_scale_dgp <- c(-0.05, 0, 0.85)

  scale_parameter_dgamma <- exp(tcrossprod(as.matrix(x), matrix(betas_scale_dgamma, nrow = 1)))
  scale_dgp <- exp(tcrossprod(as.matrix(x), matrix(betas_scale_dgp, nrow = 1)))

  sample <- extraDistr::rdgamma(n = n, shape = shape_dgamma, scale = scale_parameter_dgamma)

  dgamma_cdf <- do.call(rbind, lapply(scale_parameter_dgamma, function(x) {
    extraDistr::pdgamma(q = 0:1000, shape = shape_dgamma, scale = x)
  }))
  dgamma_pmf <- do.call(rbind, lapply(scale_parameter_dgamma, function(x) {
    extraDistr::ddgamma(x = 0:1000, shape = shape_dgamma, scale = x)
  }))

  index_quantile <- apply(dgamma_cdf, 1, function(y) {
    which.min(abs(y[(y - (1 - phi)) < 0] - (1 - phi)))
  })

  update_phi <- numeric()
  for (i in 1:nrow(dgamma_cdf)) {
    update_phi[i] <- 1 - dgamma_cdf[i, index_quantile[[i]]]
  }

  xx_max <- ncol(dgamma_cdf) - 1
  cdf_dggp <- pmf_dggp <- matrix(0, nrow = nrow(x), ncol = ncol(dgamma_cdf))
  colnames(cdf_dggp) <- colnames(pmf_dggp) <- as.character(0:xx_max)
  true_quantiles <- matrix(0, nrow = nrow(x), ncol = length(probs))
  y <- numeric(nrow(x))

  for (j in 1:nrow(x)) {
    pmf_dggp[j, 1:(index_quantile[j])] <- dgamma_pmf[j, 1:(index_quantile[j])]
    if (index_quantile[j] < (xx_max + 1)) {
      pmf_dggp[j, (index_quantile[j] + 1):(xx_max + 1)] <- update_phi[j] * ddgp(x = 0:(xx_max - index_quantile[j]), shape = xi, scale = scale_dgp[j])
    } else {
      pmf_dggp[j, index_quantile[j]] <- update_phi[j] * ddgp(x = 0, shape = xi, scale = scale_dgp[j])
    }
    cdf_dggp[j, ] <- cumsum(pmf_dggp[j, ])
    true_quantiles[j, ] <- sapply(probs, function(p) which((cdf_dggp[j, ] - p) > 0)[1] - 2)
    y[j] <- sample(x = 0:xx_max, size = 1, prob = pmf_dggp[j, ], replace = TRUE)
  }

  true_quantiles <- ifelse(true_quantiles < 0, 0, true_quantiles)
  colnames(true_quantiles) <- probs


  list(
    data = data.table::data.table(fit_data_row = 1:n, x1 = x$x1, x2 = x$x2, y = y),
    pmf_dggp = pmf_dggp,
    true_quantiles = true_quantiles
  )
}




