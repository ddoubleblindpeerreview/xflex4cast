#' @importFrom data.table data.table :=
#' @importFrom stats predict quantile sigma
#' @importFrom utils data flush.console
NULL

#' Cumulative distribution function for the Discrete Generalised Pareto distribution
#'
#' @param q quantiles
#' @param shape shape parameter
#' @param scale scale parameter
#'
#' @details Internal use.
#'
pdgp <- function(q, shape, scale){
  if (any(q<0)) {
    stop('Invalid support for x. I.e: x = 0,1,2,... \n')
  }

  if(any(scale<0)){
    stop('Invalid scale parameter.')
  }

  if(abs(shape) < 5e-2 ){
    return(1 - exp(-(1+c(q))/c(scale)))
  }
  # See page 890 from Ranjbar, Setareh, et al. 2022 to see this Equation.
  return ( 1 - ( 1+ (c(shape)*(1+c(q)))/(c(scale)))^(-1/c(shape)) )
}

#' Cumulative distribution function for the continuous Generalised Pareto distribution
#'
#' @param q quantiles
#' @param shape shape parameter
#' @param scale scale parameter
#'
#' @details Internal use.
#'
pgp <- function(q, shape, scale){
  if (any(q<0) ) {
    stop('Invalid support for x. I.e: x = 0,1,2,... \n')
  }

  if (any(scale<0)) {
    stop('Invalid scale parameter.')
  }

  p <- numeric(length(q))
  yu <- pmax(q,0)/scale

  if(abs(shape) < 1e-6){
    p <- 1-exp(-yu)
    return(p)
  }

  syu <- pmax(1+shape*yu,0)
  p <- 1 - syu^(-1/shape)

  return(p)
}



#' Calculate the quantile for a Generalised Pareto Distribution (Continuous or
#'  Discrete) fitted through a GAMLSSS model.
#'
#' @param p probability level associated with the quantile
#' @param newdata Input data to be used to generate new predictions.
#' @param p4r_tail  A `p4r_tail` object that has been fit.
#'
#' @return a quantile associated with the `p` probability level for exceedances
#' from a `p4r_tail` fit.
#'
#' @export
#'
qgp <- function(p,
                newdata,
                p4r_tail,
                sigma_new = NULL,
                shape = NULL){

  est_coef <- stats::coefficients(p4r_tail$tail_mod)

  if(!length(est_coef) %in% c(1,2,3)){
    stop("Model was not defined for other formulas rather than constant shape and univariate scale")
  }

  if(is.null(shape)){
    shape <- if(p4r_tail$tail_distribution %in% c("GPII","DGPII")){
      exp(est_coef[1])
    } else {

      if(p4r_tail$tail_distribution == "DGP0"){
        0.0
      } else {
        est_coef[1]
      }
    }
  }
  if(is.null(sigma_new)){
    if(length(est_coef)==3){
      sigma_new <- c(exp(est_coef[2]+crossprod(est_coef[3:length(est_coef)],
                                               newdata[,get(all.vars(p4r_tail$scale_formula))])))
    } else if (length(est_coef)==2){
      sigma_new <- if(p4r_tail$tail_distribution=="DGP0"){
        c(exp(est_coef[1]+crossprod(est_coef[2],
                                    newdata[,get(all.vars(p4r_tail$scale_formula))])))
      } else {
        exp(est_coef[2])
      }
    } else if(length(est_coef)==1){
      sigma_new <- exp(est_coef[1])
    } else {
      stop("Model was only defined for univariate scale dependence.")
    }
  }
  if(nrow(newdata)!=length(sigma_new) && length(sigma_new)!=1){
    stop("Estimated scale parameter is not the same dimension of newdata")
  }

  shape0bool <- abs(shape) < 1e-6

  if(shape0bool){
    if(p4r_tail$tail_distribution%in% c("GP","GPII","GP0")){
      q = -sigma_new*log(1-p)
    } else if(p4r_tail$tail_distribution%in% c("DGP","DGPII","DGP0")){
      q = ceiling(-sigma_new * log(1 - p)) - 1
    }

  } else {

    if(p4r_tail$tail_distribution%in% c("GP","GPII")){
      q = (sigma_new/shape)*((1-p)^(-shape)-1)
    } else if(p4r_tail$tail_distribution%in% c("DGP","DGPII")){
      q = ceiling((sigma_new/shape)*((1 - p)^(-shape)-1)) - 1
    }
  }

  return(q)
}


#' Function to get the quantiles of a Generalised Pareto distribution
#' as a function of the shape and scale parameters.
#'
#' @param p probability level associated with the quantile
#' @param p4r_mqgam A `p4r_tail` object that has been fit.
#' @param sigma_new a numeric value or vector for the scale parameter
#' @param shape a constant value for the shape parameter
#'
#' @return  a quantile associated with the `p` probability levelfor the
#' `p4r_mqgam` adjusted model as function of shape and scale parameter.
#'
get_qgp_using_parameters <- function(p,
                                     p4r_tail,
                                     sigma_new = NULL,
                                     shape = NULL){


  if(length(sigma_new)!=1){ # If sigma is constant;
    if(length(p)!=length(sigma_new)){
      stop("Estimated scale parameter is not the same dimension of newdata")
    }
  }

  shape0bool <- abs(shape) < 1e-6

  if(shape0bool){
    if(p4r_tail$tail_distribution %in% c("GP","GPII","GP0")){
      q = -sigma_new*log(1-p)
    } else if(p4r_tail$tail_distribution %in% c("DGP","DGPII","DGP0")){
      q = ceiling(-sigma_new * log(1 - p)) - 1
    }

  } else {

    if(p4r_tail$tail_distribution %in% c("GP","GPII")){
      q = (sigma_new/shape)*((1-p)^(-shape)-1)
    } else if(p4r_tail$tail_distribution %in% c("DGP","DGPII")){
      q = ceiling((sigma_new/shape)*((1 - p)^(-shape)-1)) - 1
    }
  }

  return(q)
}


#' Compute Quantiles for Generalized Pareto Tail Models
#'
#' Calculates the quantile function values for various continuous or discrete Generalized Pareto (GP) tail models based on given probabilities and tail parameters.
#'
#' @param p Numeric vector of probabilities for which quantiles are to be computed. Values must lie in [0,1).
#' @param tail_model A list or object containing tail model parameters, specifically a \code{margins} attribute indicating the model type. Supported margins: "GP", "GPII", "GP0", "DGP", "DGPII", "DGP0".
#' @param sigma_new Numeric scalar or vector. Scale parameter(s) for the tail model. If vector, must match length of \code{p}.
#' @param shape Numeric scalar. Shape parameter of the tail distribution.
#'
#' @return Numeric vector of quantiles corresponding to probabilities \code{p} according to the specified tail model and parameters.
#'
#' @details
#' The function implements quantile computations for both continuous GP models and their discrete analogues. For shape parameters near zero, the function uses the exponential distribution quantile form. Discrete models apply ceiling and shift adjustments to map continuous quantiles to discrete support.
#'
#' @examples
#' tail_model <- list(margins = "DGP")
#' get_qgp_using_parameters_from_tail_model(p = c(0.9, 0.95), tail_model = tail_model, sigma_new = 2, shape = 0.3)
#'
#' @export

get_qgp_using_parameters_from_tail_model <- function(p,
                                     tail_model,
                                     sigma_new = NULL,
                                     shape = NULL){


  if(length(sigma_new)!=1){ # If sigma is constant;
    if(length(p)!=length(sigma_new)){
      stop("Estimated scale parameter is not the same dimension of newdata")
    }
  }

  shape0bool <- abs(shape) < 1e-6

  if(shape0bool){
    if(tail_model$margins[1] %in% c("GP","GPII","GP0")){
      q = -sigma_new*log(1-p)
    } else if(tail_model$margins[1] %in% c("DGP","DGPII","DGP0")){
      q = ceiling(-sigma_new * log(1 - p)) - 1
    }

  } else {

    if(tail_model$margins[1] %in% c("GP","GPII")){
      q = (sigma_new/shape)*((1-p)^(-shape)-1)
    } else if(tail_model$margins[1] %in% c("DGP","DGPII")){
      q = ceiling((sigma_new/shape)*((1 - p)^(-shape)-1)) - 1
    }
  }

  return(q)
}

#' Probability Mass Function of the Discrete Generalized Pareto Distribution
#'
#' Computes the probability mass function (PMF) values of the Discrete Generalized Pareto (DGP) distribution at given points.
#'
#' @param x Numeric vector of non-negative integers at which to evaluate the PMF.
#' @param shape Numeric scalar. Shape parameter of the DGP distribution.
#' @param scale Numeric scalar. Scale parameter of the DGP distribution.
#'
#' @return Numeric vector of PMF values corresponding to \code{x}.
#'
#' @details
#' The DGP distribution is a discrete analogue of the continuous Generalized Pareto distribution. When the shape parameter is near zero, the PMF simplifies to differences of exponential terms; otherwise, it follows the standard DGP formula.
#'
#' The support of the distribution is \eqn{x \in \{0,1,2,...\}}.
#'
#' @examples
#' ddgp(x = 0:5, shape = 0.3, scale = 2)
#'
#' @export

ddgp <- function(x, shape, scale){

  if(abs(shape)<1e-6){
    return( exp(-x/scale)-(exp(-(x+1)/scale)) )
  } else {
    return((1+(shape*x)/scale)^(-1/shape)-(1+((shape*(x+1))/scale))^(-1/shape))
  }

}


