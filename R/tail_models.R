#' @title Compute Exceedances Using Jittered Model
#' @description
#' Internal function to compute exceedances using a jittered model prediction based on a quantile threshold.
#'
#' @details
#' This function uses jittered quantile predictions to compute exceedances, defined as the difference between the observed values and the predicted quantile. Observations with exceedances greater than or equal to zero are stored.
#'
#' @param model A list representing the model object. Must contain elements like `data`, `quantile_threshold`, and `mqgam_formula`.
#' @param mqgam_jitter A list of jittered MQGAM models, indexed by quantile thresholds.
#'
#' @return The input `model` object with `exceedance_data` and `alpha_right` elements updated.
#'
#' @keywords internal
#' @noRd
get_exceedances_jitter <- function(model, mqgam_jitter){

  # Select the jittered models with respect to the quantile_threshold
  jittered_quantile_threshold <- mqgam_jitter[[as.character(model$quantile_threshold)]]

  quantile_prediction <- predict.mqgam_jitter(mqgam_jitter = mqgam_jitter,newdata = data,qu = model$quantile_threshold)$mean_jittering

  if(nrow(data)!=length(quantile_prediction)){
    stop("Not compatible mqgam_jitter and data ")
  }

  model$exceedance_data <- model$data[,exceedance:=get(model$mqgam_formula[[2]])-ceiling(quantile_prediction)][exceedance>=0,]
  model$alpha_right <- quantile_prediction

  return(model)

}

#' @title Compute Exceedances Using qgam Prediction
#' @description
#' Internal function to compute exceedances using the `qgam::qdo()` method for quantile regression prediction.
#'
#' @details
#' Computes exceedances as the difference between observed values and the quantile predictions. Keeps only those observations where the exceedance is non-negative.
#'
#' @param model A list representing the model object. Must contain elements like `data`, `quantile_threshold`, `mqgam_fit`, and `mqgam_formula`.
#'
#' @return The input `model` object with updated `exceedance_data` and `alpha_right`.
#'
#' @keywords internal
#' @noRd
get_exceedances_qgam <- function(model){

  # Select the jittered models with respect to the quantile_threshold
  quantile_prediction <- qgam::qdo(model$mqgam_fit,
                                  qu = as.numeric(model$quantile_threshold),
                                  newdata = model$data,fun = predict)

  if(nrow(model$data)!=length(quantile_prediction)){
    stop("Not compatible mqgam_jitter and data ")
  }

  model$exceedance_data <- model$data[,exceedance:=get(model$mqgam_formula[[2]])-ceiling(quantile_prediction)][exceedance>=0,]
  model$alpha_right <- quantile_prediction

  return(model)

}

#' @title Fit Tail Distribution Model
#' @description
#' Internal function to fit a tail distribution model to exceedance data using GJRM.
#'
#' @details
#' Based on the user-defined `tail_distribution` and `scale_formula`, this function fits a tail model to previously computed exceedance data using `GJRM::gamlss`.
#'
#' @param model A list representing the model object. Must contain `exceedance_data`, `scale_formula`, `tail_distribution`, and `quantile_threshold`.
#'
#' @return The input `model` object with the additional `tail_mod` element containing the fitted model.
#'
#' @keywords internal
#' @noRd
fit_tail <- function(model) {

  if (!(all(all.vars(model$scale_formula[[2]]) %in% colnames(model$data)) || model$scale_formula[[2]]==1)) {
    stop(paste0(model$scale_formula[[2]], " formula is invalid for the p4r_mqgam$data object."))
  }

  if(is.null(model$quantile_threshold)){
    stop("Specify a value for the quantile value to be used as threshold.")
  }

  if (model$quantile_threshold > max(model$mqgam_quantiles)) {
    stop("Quantile threshold not properly defined for the model model.")
  }


  gjrm_formula <- if(model$tail_distribution == "DGP0"){
    list(stats::update(model$scale_formula,exceedance ~ .))
  } else {
    list(exceedance ~ 1, model$scale_formula)
  }


  gjrm_gamlss <- GJRM::gamlss
  gjrm_call <- call("gjrm_gamlss",
                    formula = gjrm_formula,
                    family = model$tail_distribution,
                    data = model$exceedance_data
  )

  model$tail_mod <- eval(gjrm_call)

  return(model)
}
