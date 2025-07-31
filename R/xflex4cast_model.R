#' @importFrom data.table data.table :=
#' @importFrom stats predict quantile sigma
#' @importFrom utils data flush.console
NULL

#' @title Define xflex4cast Model Parameters
#' @description
#' Internal function to define and initialize parameters for an xflex4cast model.
#'
#' @param mqgam_formula A formula for the MQGAM model.
#' @param mqgam_quantiles Quantile levels to be fitted. Default: c(0.05, 0.25, 0.5, 0.75, 0.95). If `fit_tail=TRUE` values
#' above the quantile_threshold will be modelled by the DGP, while below to it the `mqgam`.
#' @param fit_tail Logical, whether to fit a tail model for high quantiles. Default: TRUE.
#' @param scale_formula Formula used in tail model (GJRM). Default: ~ 1.
#' @param quantile_threshold Threshold quantile for tail modeling. Default: 0.95.
#' @param tail_distribution String indicating the family to be used in tail modeling. Default: "DGP".
#' @param maximum_quantile_level Maximum quantile level to be used. Default: 0.999.
#'
#' @return A list with class `xflex4cast`.
#'
#' @keywords internal
#' @export
define_xflex4cast <- function(mqgam_formula,
                                mqgam_quantiles = c(0.05,0.25,0.5,0.75,0.95),
                                fit_tail = TRUE,
                                scale_formula = ~ 1,
                                quantile_threshold = 0.95,
                                tail_distribution = "DGP",
                                maximum_quantile_level = 0.999
                                ){

    # Updating the maximum quantile level accordingly
    maximum_quantile_level <- max(c(mqgam_quantiles,maximum_quantile_level))

    # Including the quantile_threshold at the mqgam_quantiles
    mqgam_quantiles <- unique(sort(c(mqgam_quantiles,quantile_threshold)))

    model <- list(mqgam_formula = mqgam_formula,
                  mqgam_quantiles = mqgam_quantiles,
                  fit_tail = fit_tail,
                  scale_formula = scale_formula,
                  quantile_threshold = quantile_threshold,
                  tail_distribution = tail_distribution,
                  maximum_quantile_level = maximum_quantile_level)

    class(model) <- "xflex4cast"
    return(model)
}


#' @title Fit xflex4cast Model
#' @description
#' Internal method to fit a `xflex4cast` model using MQGAM and optionally fit a tail model.
#'
#' @param model A `xflex4cast` object.
#' @param data A `data.frame` or `data.table` with the training data.
#' @param index Column name for unique identifier used in cross-validation. If `NULL`, uses row numbers.
#'
#' @return A fitted `xflex4cast` object.
#'
#' @export

fit.xflex4cast <- function(model,
                           data,
                           index = NULL) {


  # Add training data with index
  if( is.null(index) ) {
    data[,fit_data_row:=.I]
    index <- "fit_data_row"
  } else {

    if(!(index %in% colnames(data))){
      stop('The inserted "index" is not a column name from the data.')
    }

    if( nrow(data) != length(unique(data[[index]])) ){
      stop("All elements of data[[index]] must be unique." )
    }
  }


  if(!(all(all.vars(model$scale_formula[[2]]) %in% colnames(data)) || model$scale_formula[[2]]==1)){
    stop('Not all terms from p4r_tail$scale_formula are in the data.')
  }

  model$data_index <- index
  formula_mqgam <- stats::update(model$mqgam_formula,paste0("~.+", index))

  if(length(all.vars(model$scale_formula[[2]]))!=0){
    formula_mqgam_tail <- stats::update(formula_mqgam,paste0("~.+",paste0(all.vars(model$scale_formula[[2]]),collapse = "+")))
  } else {
    formula_mqgam_tail <- formula_mqgam
  }

  model$data <- data.table::as.data.table(stats::get_all_vars(formula_mqgam_tail,data))



  model$mqgam_fit <- qgam::mqgam(model$mqgam_formula,
                                 data = data, qu = model$mqgam_quantiles,
                                 control = list(verbose = FALSE, progress = FALSE))


  if(model$fit_tail){
    # Getting exceed data
    model <- get_exceedances_qgam(model)
  }

  # Fitting the tail model
  if(model$fit_tail){
    model <- fit_tail(model = model)
    if( sum(as.numeric(summary(model$tail_mod)$indx == FALSE)) != 0) warning("Condition violated. You need to check the fit to the tail model.", call. = FALSE)
  }



  return(model)

}



#' @title Predict from xflex4cast Model
#' @description
#' Predict quantiles from a fitted `xflex4cast` model. Returns long-format table.
#'
#' @param model A fitted `xflex4cast` model.
#' @param newdata A `data.frame` or `data.table` of new data to predict on.
#' @param indexes Names of identifying columns to track predictions. Default: `index` with row numbers.
#'
#' @return A `data.table` in long format with columns: indexes, quantile, and value.
#'
#' @export
predict.xflex4cast <- function(model,
                               newdata,
                               indexes = NULL){

  if (is.null(indexes)) {
    index <- newdata[,.I]
    indexes <- "index"
  } else {
    index <- newdata[,.SD,.SDcols=indexes]
  }

  predictions_long <- data.table()

  pred_mqgam <- qgam::qdo(obj = model$mqgam_fit,
                          qu = model$mqgam_quantiles,
                          fun = predict,newdata = newdata)

  for( i in 1:length(model$mqgam_quantiles) ) {
    predictions_long <- rbind(
      predictions_long,
      cbind(index,
            data.table(quantile = model$mqgam_quantiles[i],
                        value = pred_mqgam[[i]])
      )
    )
  }




  predictions_long <- quantile_crossing_correction(predictions = predictions_long,
                                                   indexes = indexes)
  if(model$fit_tail){
    predictions_long <- predictions_long[quantile<=model$quantile_threshold,]
  }

  if(!is.null(model$tail_mod) ){

    tail_predictions_long <- data.table()
    if(is.null(model$tail_quantiles)){
      tail_quantiles <- model$mqgam_quantiles[model$mqgam_quantiles>model$quantile_threshold]
    } else {
      tail_quantiles <- model$tail_quantiles
    }
    tail_quantiles_predictions <- vector("list",length(tail_quantiles))


    quantile_threshold <- predictions_long[,.(quantile_threshold = stats::approx(y = quantile,
                                                    x = value,
                                                    xout = floor(max(value)),
                                                    yleft = 0, yright = max(quantile),
                                                    ties = min)$y), by = indexes]


    dgp_quantile_level <- lapply(tail_quantiles,
                                 function(quantile_level){ (quantile_level-quantile_threshold$quantile_threshold)/(1-quantile_threshold$quantile_threshold) })

    names(dgp_quantile_level) <- tail_quantiles


    tail_quantiles_predictions <- lapply(dgp_quantile_level, function(prob){
      GJRM::pred.gp(x = model$tail_mod,p = prob, newdata = newdata,n.sim = 2)$qp
    })

    quantile_threshold_values <- predictions_long[quantile==model$quantile_threshold,][["value"]]



    for( j in 1:length(tail_quantiles_predictions) ) {
      tail_predictions_long <- rbind(
        tail_predictions_long,
        cbind(index,
              data.table(quantile = tail_quantiles[j],
                          value = tail_quantiles_predictions[[j]] + ceiling(quantile_threshold_values)-1))
        )
    }




    predictions_long <- rbind(predictions_long,tail_predictions_long)
  }

  predictions_long <- predictions_long[value<0,value:=0]
  return(predictions_long)
}


#' @title Forecast with xflex4cast Model
#' @description
#' Generate forecast from an `xflex4cast` model, returning a probability mass function.
#'
#' @param model A fitted `xflex4cast` model.
#' @param newdata Data to use for forecasting.
#' @param indexes Identifier columns in `newdata` to keep.
#' @param max_quantile_predicted Maximum quantile for tail extrapolation. Default: NULL.
#'
#' @return A long-format `data.table` with columns: index, x, cdf, and pmf.
#'
#' @export
forecast <- function(model,
                     newdata,
                     indexes,
                     max_quantile_predicted = NULL){


  if (is.null(indexes)) {
    index <- newdata[,.I]
    indexes <- "index"
  } else {
    index <- newdata[,.SD,.SDcols=indexes]
  }

  predictions_long <- data.table()

  pred_mqgam <- qgam::qdo(obj = model$mqgam_fit,
                          qu = model$mqgam_quantiles,
                          fun = predict,newdata = newdata)

  for( i in 1:length(model$mqgam_quantiles) ) {
    predictions_long <- rbind(
      predictions_long,
      cbind(index,
            data.table(quantile = model$mqgam_quantiles[i],
                        value = pred_mqgam[[i]])
      )
    )
  }




  predictions_long <- quantile_crossing_correction(predictions = predictions_long,
                                                   indexes = indexes)

  cdf <- if(!is.null(max_quantile_predicted)){
        predictions_long[,
                     .(x = 0:floor(min(max(value),max_quantile_predicted)),cdf = stats::approx(
                       y = quantile, x = value,
                       xout = 0:floor(min(max(value),max_quantile_predicted)),
                       yleft = 0, yright = max(quantile),ties = min)$y),
                     by = indexes]
  } else {
        predictions_long[,
                         .(x = 0:floor(min(max(value),max_quantile_predicted)),cdf = stats::approx(
                           y = quantile, x = value,
                           xout = 0:floor(min(max(value),max_quantile_predicted)),
                           yleft = 0, yright = max(quantile),ties = min)$y),
                         by = indexes]

  }

  if(!is.null(model$tail_mod) ){

    alpha_right <- unique(cdf[,.(cdf_alpha_right = max(cdf),q_alpha_right = max(x)), by = indexes][,.SD,.SDcols = c(indexes,"cdf_alpha_right","q_alpha_right")])
    alpha_tilde <- unique(cdf[,alpha_tild_right := (model$maximum_quantile_level-max(cdf))/(1-max(cdf)), by = indexes][,.SD,.SDcols = c(indexes,"alpha_tild_right")])

    # Results from the tail model
    tail_result <- pred.gp(p = model$maximum_quantile_level,x = model$tail_mod,newdata = newdata,n.sim = 2)
    if(is.null(tail_result$xi)){
      tail_result$xi <- rep(0,length(tail_result$sigma))
    }

    alpha_tilde[,sigma := tail_result$sigma]
    alpha_tilde[,shape := tail_result$xi]
    alpha_tilde[,max_quantile := get_qgp_using_parameters_from_tail_model(p = alpha_tild_right,
                                                                          sigma_new = sigma,
                                                                          shape = shape, tail_model = model$tail_mod),
                by = indexes]
    cdf_tail <- if(is.null(max_quantile_predicted)){
      alpha_tilde[,.(x_tail = 0:max_quantile,
                     cdf_tail = pdgp(q = 0:max_quantile,shape = shape,scale = sigma)),
                  by = indexes]
    } else {
     alpha_tilde[,.(x_tail = 0:max_quantile_predicted,
                              cdf_tail = pdgp(q = 0:max_quantile_predicted,shape = shape,scale = sigma)),
                           by = indexes]
    }

    cdf_long_tail <- merge(cdf_tail,alpha_right,by = indexes)[,.(x_tail=x_tail+q_alpha_right+1,cdf_tail=cdf_alpha_right+(1-cdf_alpha_right)*cdf_tail), by = indexes]
    setnames(cdf_long_tail,c("x_tail","cdf_tail"),c("x","cdf"))

    all_cdf <- rbind(cdf[,.SD,.SDcols = c(indexes,"x","cdf")],cdf_long_tail)

    if( stats::coefficients(model$tail_mod)[1] < 0) {
      all_cdf[,cdf:=ifelse(is.nan(cdf),max(cdf,na.rm = TRUE),cdf),by = eval(model$data_index)]
    }


  } else {

    # This is used to calculate the pmf beyond the right-limit from the quantile
    if(!is.null(max_quantile_predicted)){

        cdf_long_tail <- data.table()
        # Avoid duplicating before
        for(i in 1:nrow(newdata)){
            if(max(cdf[fit_data_row==i,"x"])<max_quantile_predicted){
              cdf_long_tail_aux <- unique(cdf[fit_data_row==i,.(cdf_tail = max(cdf),x_tail = (max(x)+1):max_quantile_predicted), by = indexes][,.SD,.SDcols = c(indexes,"cdf_tail","x_tail")])
              setnames(cdf_long_tail_aux,c("x_tail","cdf_tail"),c("x","cdf"))
              cdf_long_tail <- rbind(cdf_long_tail,cdf_long_tail_aux)
            }
        }

        all_cdf <- rbind(cdf[,.SD,.SDcols = c(indexes,"x","cdf")],cdf_long_tail)

      } else {
          all_cdf <- cdf
      }


  }


  # Ordering according to the CDF
  all_cdf <- all_cdf[order(all_cdf,get(indexes))]

  all_cdf[, pmf := c(cdf[1], diff(cdf)), by = indexes]

  # Returning values up to the upper limit
  if(!is.null(max_quantile_predicted)){
    all_cdf <- all_cdf[x<=max_quantile_predicted,]
  }

  return(all_cdf)
}


#' @title Deprecated: Predict from Poisson Model
#' @description
#' Deprecated function to generate Poisson quantile predictions.
#'
#' @param newdata A `data.table` containing prediction input.
#' @param indexes Identifier columns in the input.
#' @param quantiles Quantile levels to compute.
#'
#' @return A `data.table` in long format with quantile predictions.
#'
#' @keywords internal
#' @noRd
predict_poisson <- function(newdata,
                            indexes,
                            quantiles){

  if (is.null(indexes)) {
    index <- newdata[,.I]
    indexes <- "index"
  } else {
    index <- newdata[,.SD,.SDcols=indexes]
  }

  lambda_new <- get_lambda_poisson(newdata)

  predictions_long <- data.table()

  for(i in 1:length(quantiles)){
    predictions_long <- rbind(
      predictions_long,
      cbind(index,
            data.table(quantile = quantiles[i],
                       value = stats::qpois(p = quantiles[i],lambda = lambda_new)))
    )
  }

  return(predictions_long)

}

#' @title Deprecated: Poisson Forecasting Function
#' @description
#' Generate Poisson forecasts for a given dataset, returning the full pmf.
#'
#' @param newdata A `data.table` with input data.
#' @param indexes Identifier column names.
#' @param max_quantile_predicted Upper bound of quantiles. Default: NULL.
#'
#' @return A long-format `data.table` with columns: index, x, cdf, and pmf.
#'
#' @keywords internal
#' @noRd
poisson_forecast <- function(newdata,
                             indexes,
                             max_quantile_predicted = NULL){


  if (is.null(indexes)) {
    index <- newdata[,.I]
    indexes <- "index"
  } else {
    index <- newdata[,.SD,.SDcols=indexes]
  }

  predictions_long <- data.table()

  lambda_new <- get_lambda_poisson(x = newdata)

  for( i in 1:length(lambda_new) ) {
    predictions_long <- rbind(
      predictions_long,
      cbind(index[i],
            data.table(x = 0:max_quantile_predicted,
                       cdf = stats::ppois(q = 0:max_quantile_predicted,lambda = lambda_new[i]),
                       pmf = stats::dpois(x = 0:max_quantile_predicted,lambda = lambda_new[i]))
      )
    )
  }

  # Ordering according to the CDF
  all_cdf <- predictions_long

  return(all_cdf)
}
