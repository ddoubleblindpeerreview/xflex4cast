#' @importFrom data.table data.table :=
#' @importFrom stats predict quantile sigma
#' @importFrom utils data flush.console
NULL

#' Join prediction and actual data
#'
#' @param predictions probabilistic fault forecasts, quantiles in long
#' format `data.table`.
#' @param data actual fault data with columns `indexes` and column `actual`
#' containing outcome corresponding to predictions.
#' @param valid_time_index common column name containing the index on which
#' to merge `predictions` and `data`. If `NULL`, predictions and data will be
#' merged using a common column name if only one exists.
#'
#' @description This function joints predictions and actual data in a single
#' data table.
#'
#' @returns A `data.table` of merged predictions with new column of
#' corresponding "actuals".
#'
#' @keywords internal
#' @noRd

merge_pred_actuals <- function(predictions, data, valid_time_index) {
  # Check if 'predictions' and 'data' are data.tables
  if (!is.data.table(predictions) || !is.data.table(data)) {
    stop("Input must be data.table")
  }

  # Check if the required columns are present in 'predictions'
  if (!all(c(valid_time_index,"quantile", "value") %in% names(predictions))) {
    stop("Missing required columns in 'predictions', see ?merge_pred_actuals()")
  }

  # Check if the required columns are present in 'data'
  if (!"actual" %in% names(data)) {
    stop("Missing required column 'actual' in 'data'")
  }

  if( is.null(valid_time_index) ){
    valid_time_index <- intersect(colnames(predictions), colnames(data))
    if( length(valid_time_index)!=1 ){
      stop("valid_time_index not provided and no single common column name to merge on.")
    }
  }

  if( !all(predictions[,get(valid_time_index)] %in% data[,get(valid_time_index)]) ) {
    warning("Corresponding actual/outturn not found in data for some predictions. Some evaluation results may be NA.")
  }

  predictions <- merge(predictions,
                       data[,.SD,.SDcols = c(valid_time_index,"actual")],
                       by = paste0(valid_time_index),
                       all.x = T)

  return(predictions)
}

#' Pinball Loss
#'
#' @param predictions probabilistic fault forecasts, quantiles in long
#' format `data.table`.
#' @param data actual fault data with columns `indexes` and column `actual`
#' containing outcome corresponding to predictions.
#' @param valid_time_index common column name containing the index on which
#' to merge `predictions` and `data`. If `NULL`, a merge will be attempted on
#' if there a single common column name exists.
#' @param indexes column names in `predictions` used to group average pinball
#' scores by.
#'
#' @description This function takes as input a data.table of predictions
#' (with columns `quantile` and `value`) and a `data.table` of
#' actual values (with column `actual`). Both should contain an index column
#' to be used for merging, which may be specified using `valid_time_index`. It computes
#' the average pinball loss grouped by quantiles and `indexes` and returns the
#' results as a `data.table`.
#'
#' @returns a `data.table` of pinball scores with the same indexes as `predictions`
#' (dtm(s) and quantiles).
#'
#' @export
pinball <- function(predictions,
                    data,
                    valid_time_index=NULL,
                    indexes = NULL) {

  # Check if 'predictions' and 'data' are data.tables
  if (!is.data.table(predictions) || !is.data.table(data)) {
    stop("Input must be data.table")
  }

  # Check if the required columns are present in 'predictions'
  if (!all(c("quantile", "value") %in% names(predictions))) {
    stop("Missing required columns 'quantile' and 'value' in 'predictions'")
  }

  # Check if the required columns are present in 'data'
  if (!"actual" %in% names(data)) {
    stop("Missing required column 'actual' in 'data'")
  }

  predictions = merge_pred_actuals(predictions,data,valid_time_index)

  predictions[,pinball := (actual-value)*quantile*(actual>=value) +
                (value-actual)*(1-quantile)*(actual<value)]

  return(predictions[,.(pinball=mean(pinball)),by=c("quantile",indexes)])

}


#' Reliability/Calibration
#'
#' @param predictions probabilistic fault forecasts, quantiles in long format
#' `data.table`.
#' @param data actual fault data with columns `indexes` and column `actual`
#' containing outcome corresponding to predictions.
#' @param valid_time_index common column name containing the index on which
#' to merge `predictions` and `data`.
#' @param indexes column names in `predictions` used to group result by,
#' defaults to `NULL` (no grouping)
#' @param exclude_zero_pred exclude predicted quantiles that equal zero. Recommended
#' of the target variable is zero-inflated. NB: other treatment for this should
#' be explored.
#'
#' @description This function computes empirical coverage for each quantile
#' in `predictions` and given observations `data` that should include a common
#' column of indexes for merging, which may be specified by `valid_time_index`.
#'
#' @returns a `data.table` containing the empirical coverage for each quantile
#' in `predictions`.
#'
#' @export
reliability <- function(predictions,
                        data,
                        valid_time_index=NULL,
                        indexes=NULL,
                        exclude_zero_pred = T) {

  # Check if 'predictions' and 'data' are data.tables
  if (!is.data.table(predictions) || !is.data.table(data)) {
    stop("Input must be data.table")
  }

  # Check if the required columns are present in 'predictions'
  if (!all(c("quantile", "value") %in% names(predictions))) {
    stop("Missing required columns 'quantile' and 'value' in 'predictions'")
  }

  # Check if the required columns are present in 'data'
  if (!"actual" %in% names(data)) {
    stop("Missing required column 'actual' in 'data'")
  }

  predictions = merge_pred_actuals(predictions,data,valid_time_index)


  reliability <- predictions[if(exclude_zero_pred) {value!=0} else {1:.N},
                             .(empirical = mean(actual<=value)),
                             by=c("quantile",indexes)]

  return(reliability)

}

#' Sharpness (average interval width)
#'
#' @param predictions probabilistic fault forecasts, quantiles in long format
#' `data.table`.
#' @param indexes column names in `predictions` used to group result by,
#' defaults to `NULL` (no grouping).
#'
#' @description This function computes average interval width for all symmetric
#' pairs of quantiles (e.g. average width of 80% interval based on 10% and 90%
#' quantiles) for `predictions`, optionally grouped by `indexes`. The result is
#' returned as a `data.table`.
#'
#' @returns a `data.table` with the same time indexes as `predictions` and the
#' interval width for all pairs of quantiles symmetric about q50.
#'
#' @export
sharpness <- function(predictions,valid_time_index, indexes=NULL) {

  # Check if 'predictions' and 'data' are data.tables
  if (!is.data.table(predictions)) {
    stop("Input must be data.table")
  }

  # Check if the required columns are present in 'predictions'
  if (!all(c("quantile", "value") %in% names(predictions))) {
    stop("Missing required columns in 'predictions'. 'quantile' and 'value' are required.")
  }

  unique_quantiles <- predictions[,unique(quantile)]
  symetric_intervals <- c()
  for( q in  unique_quantiles[unique_quantiles<0.5] ) {
    if( (1-q) %in% unique_quantiles ) {
      symetric_intervals <- c(symetric_intervals,1-2*q)
    }
  }

  if(length(symetric_intervals)==0){
    warning("No symmetric quantiles available to construct intervals. Cannot compute sharpness.")
    return(NULL)
  }

  interval_width <- data.table()
  for( int in symetric_intervals ) {
    int_width <- predictions[quantile %in% c(0.5+int/2,0.5-int/2),
                             .(int_width = range(value)),by=c(valid_time_index,indexes)]

    interval_width <- rbind(interval_width,
                            int_width[,.(interval = int,
                                         width = mean(int_width)),
                                      by=indexes])
  }

  return(interval_width)

}


#' Convert to probabilistic category forecasts
#'
#' @param predictions probabilistic fault forecasts, quantiles in long format
#' `data.table`
#' @param categories category definitions, a vector of length 2
#' `c(green-amber,amber-red)`.
#' @param indexes column names to group results by
#' @param data (optional) actual fault data with date-time column(s) and column
#' `actual` containing outcome corresponding to predictions. A merge will be
#' attempted between `predictions` and `data` on all column names in both
#' `indexes` and `data`.
#'
#' @description this function converts quantile predictions into three
#' categorical predictions.
#'
#' @returns a `data.table` of with the same times indexes as `predictions` and
#' probabilistic forecasts of each category being realised (columns `category`
#' and `probability`).
#'
#' @keywords internal
#' @noRd
quantiles_to_categories <- function(
    predictions,
    categories,
    indexes = c("dtm"),
    data = NULL) {

  # Check if 'predictions' is a data.tables
  if (!is.data.table(predictions)) {
    stop("Input must be data.table")
  }
  # Check that 'categories' has length 2
  if (length(categories) != 2) {
    stop("Input must have length 2")
  }

  cat_preds <- predictions[,.(cat_1 = stats::approx(
    y=quantile,x = value,
    xout = categories[1],
    yleft = 0,yright = max(quantile),ties=min)$y),
    by=indexes]

  cat_preds <- merge(cat_preds,
                     predictions[,
                                 .(cat_2_cdf = stats::approx(
                                   y=quantile,x = value,
                                   xout = categories[2],
                                   yleft = 0,yright = max(quantile),ties=min)$y),
                                 by=indexes],
                     by=indexes)


  cat_preds[,cat_2 := cat_2_cdf-cat_1]
  cat_preds[,cat_3 := 1 - cat_2_cdf]


  # Tidy-up predictions at limits of model (max estimated quantile)
  max_quantile <- max(predictions$quantile)
  cat_preds[cat_1==max_quantile, cat_2:=1-max_quantile]
  cat_preds[cat_1==max_quantile, cat_3:=0]
  cat_preds[,cat_2_cdf:=NULL]

  cat_preds <- data.table::melt(data = cat_preds, id.vars = indexes,
                    variable.name = "category",
                    value.name = "probability")

  if( !is.null(data) ) {

    cat_preds <- merge(cat_preds, data,
                       by = indexes[indexes%in%colnames(data)],
                       all.x=T)

    if( cat_preds[,sum(is.na(actual))]>0 ) {
      warning("Corresponding actual/outturn not found in data for some cat_preds. Some evaluation results may be NA.")
    }

    cat_preds[,actual_cat:=F]
    cat_preds[actual<=categories[1] & category=="cat_1",actual_cat:=T]
    cat_preds[actual>categories[1] & actual<=categories[2] & category=="cat_2",
              actual_cat:=T]
    cat_preds[actual>categories[2] & category=="cat_3",actual_cat:=T]
  }

  return(cat_preds)

}



#' Convert probability mass function to probabilistic category forecasts
#'
#' @param predictions probabilistic mass function from probabilistic fault forecast
#' with a mixture model in long format `data.table`.
#' @param categories category definitions, a vector of length 2
#' `c(green-amber,amber-red)`.
#' @param indexes column names to group results by
#' @param data (optional) actual fault data with date-time column(s) and column
#' `actual` containing outcome corresponding to predictions. A merge will be
#' attempted between `predictions` and `data` on all column names in both
#' `indexes` and `data`.
#'
#' @description this function converts the probability massa function description
#' into categorical predictions.
#'
#' @returns a `data.table` with the same times indexes as `predictions` and
#' probabilistic forecasts of each category being realised (columns `category`
#' and `probability`).
#' @keywords internal
#' @noRd
pmf_to_categories <- function(predictions,
                              categories,
                              indexes = NULL,
                              data = NULL){

  # Check if 'predictions' is a data.tables
  if (!is.data.table(predictions)) {
    stop("Input must be data.table")
  }

  # Check that 'categories' has length 2
  if (length(categories) != 2) {
    stop("Input must have length 2")
  }


  predictions[,max_quantiles:=ifelse(max(cumsum(quantile))>1,
                                     1,max(cumsum(quantile))),
              by = indexes]

  predictions[value<=categories[1],cat_1 := sum(quantile),
              by = indexes]
  predictions[value<=categories[2],cat_2_cdf := sum(quantile),
              by = indexes]

  cat_preds <- unique(predictions[,.SD,.SDcols = c(indexes,'cat_1','cat_2_cdf','max_quantiles')],by = indexes)
  cat_preds[,cat_2:=cat_2_cdf-cat_1]
  cat_preds[,cat_3 := 1 - cat_2_cdf]

  cat_preds[cat_1==max_quantiles, cat_2:=1-max_quantiles]
  cat_preds[cat_1==max_quantiles, cat_3:=0]
  cat_preds[,cat_2_cdf:=NULL]
  cat_preds <- cat_preds[,.SD,.SDcols = c(indexes,'cat_1','cat_2','cat_3')]

  cat_preds <- data.table::melt(data = cat_preds, id.vars = indexes,
                    variable.name = "category",
                    value.name = "probability")

  if( !is.null(data) ) {

    cat_preds <- merge(cat_preds, data,
                       by = indexes[indexes%in%colnames(data)],
                       all.x=T)

    if( cat_preds[,sum(is.na(actual))]>0 ) {
      warning("Corresponding actual/outturn not found in data for some cat_preds. Some evaluation results may be NA.")
    }

    cat_preds[,actual_cat:=F]
    cat_preds[actual<=categories[1] & category=="cat_1",actual_cat:=T]
    cat_preds[actual>categories[1] & actual<=categories[2] & category=="cat_2",
              actual_cat:=T]
    cat_preds[actual>categories[2] & category=="cat_3",actual_cat:=T]
  }

  return(cat_preds)

}

#' Get confusion matrix from categorical predictions and actuals
#'
#' @param category_predictions a `data.table` containing probabilistic
#' predictions and actuals, as produced by `quantiles_to_categories()`.
#' @param index name of index/timestamp column
#' @param cat_1_above UI displays Green when probability of `cat_1` is >=
#' `cat_1_above`. Default = 0.95.
#' @param cat_3_above UI displays Red when `P(cat_3)>=P(cat_2)` (always) and
#' when `P(cat_3)` is >= `cat_3_above`. Default = 1, i.e. only first criteria
#' used.
#'
#' @description Calculates a stylised confusion matrix for categorical
#' predictions.
#'
#' @returns A 3x3 confusion matrix with a row for each primary prediction
#' (as displayed on the user interface) and column for each outcome.
#' @keywords internal
#' @noRd
get_confusion_matrix <- function(category_predictions,
                                 index="dtm",
                                 cat_1_above=0.8,
                                 cat_3_above=0.2){


  actual_category <- category_predictions[actual_cat==T,.(index=get(index),actual_category=category)]
  category_predictions <- data.table::dcast(category_predictions,
                                            formula = paste0(index,"~category"),
                                            value.var = "probability",
                                            fun.aggregate = sum)

  category_predictions <- merge(category_predictions,actual_category,
                                by.x = index,
                                by.y="index")

  category_predictions[cat_1>=cat_1_above,ui_category:="cat_1"]
  category_predictions[cat_1<cat_1_above & cat_2>cat_3 & cat_3 < cat_3_above, ui_category:="cat_2"]
  category_predictions[cat_1<cat_1_above & (cat_3>=cat_2 | cat_3 >= cat_3_above), ui_category:="cat_3"]

  valid_category_levels <- sort(union(unique(category_predictions$ui_category),
                                      levels(category_predictions$actual_category)))

  output <- table(factor(category_predictions$ui_category,
                         levels = valid_category_levels),
                  factor(category_predictions$actual_category,
                         levels = valid_category_levels),dnn = c("UI Display","Actual"))


  return(output)
}


#' Default set of evaluate metrics for P4R models
#'
#' @param predictions probabilistic fault forecasts, quantiles in long format
#' `data.table`
#' @param data Validation data (same format as training data)
#' @param valid_time_index Must be provided for results to include sharpness.
#' @param indexes indexes to group by
#'
#' @description function to perform default evaluation, e.g., as part of
#' `fit.p4r_mqgam()`.
#'
#' @export
#' @keywords internal
#' @noRd
default_evaluation <- function(predictions,
                               data,
                               valid_time_index = NULL,
                               indexes = NULL) {

  evaluation_results <- list(
    pinball = pinball(predictions, data, valid_time_index, indexes),
    reliability = reliability(predictions, data, valid_time_index, indexes),
    sharpness = try(sharpness(predictions, valid_time_index, indexes)))


  return(evaluation_results)
}

#' Categorical forecast evaluation for P4R models
#'
#' @param predictions probabilistic fault forecasts, quantiles in long format
#' `data.table`
#' @param data actual fault data with date-time column(s) and column
#' `actual` containing outcome corresponding to predictions. A merge will be
#' attempted between `predictions` and `data` on all column names in both
#' `indexes` and `data`.
#' @param categories category definitions, a vector of length 2
#' `c(green-amber,amber-red)`.
#' @param indexes a list with elements `valid_time` and `base_time`.
#' `base_time` my be `NULL`.
#' @param lead_time_units passed to `difftime()` to specify units for forecast
#' lead-time in results
#' @param ... Additional arguments passed to `get_confusion_matrix()`
#'
#' @details Calculates confusion matrix by regulatory year and lead-time, and
#' ROC by lead-time.
#'
#' @returns A list containing confusion matrices and `pROC` objects.
#'
#' @export
default_category_evaluation <- function(
    predictions, data, categories,
    indexes=list(valid_time="dtm",
                 base_time="dtm_basetime"),
    lead_time_units = "hours",
    ...){


  results <- list()

  if(is.null(indexes$base_time)){
    predictions[,lead_time := "Reanalysis"]
  }else{
    predictions[,lead_time := difftime(get(indexes$valid_time),get(indexes$base_time),units = lead_time_units)]
  }

  predictions[,regulatory_year:=get_regulatory_year(get(indexes$valid_time))]

  valid_regulatory_years <- intersect(predictions[,unique(regulatory_year)], unique(get_regulatory_year(predictions[,get(indexes$base_time)])))

  for( lt in predictions[,unique(lead_time)] ) {

    category_predictions <- quantiles_to_categories(
      predictions[lead_time==lt,],
      categories,
      indexes = indexes$valid_time,
      data)

    results[["all_data"]][[paste0(lt,lead_time_units)]][["confusion_matrix"]] <- get_confusion_matrix(
      category_predictions = category_predictions,
      index = indexes$valid_time,...)

    for( i in 1:3 ){
      results[["all_data"]][[paste0(lt,lead_time_units)]][[paste0("roc_", i)]] <- try(pROC::roc(
        response=category_predictions[category==paste0("cat_",i),actual_cat],
        predictor=category_predictions[category==paste0("cat_",i),probability],
        quiet=T))
    }

    for( regyear in valid_regulatory_years ) {

      if( nrow(predictions[lead_time==lt & regulatory_year==regyear,])==0 ) {
        next
      }

      category_predictions <- quantiles_to_categories(
        predictions[lead_time==lt & regulatory_year==regyear,],
        categories,
        indexes = indexes$valid_time,
        data)

      results[[paste0("RY_",regyear)]][[paste0(lt,lead_time_units)]][["confusion_matrix"]] <- get_confusion_matrix(
        category_predictions = category_predictions,index = indexes$valid_time,...)
    }

  }

  return(results)
}

#' Check feature bias between training feature and forecast features
#'
#' @param model description
#' @param ensemble_nwp description
#' @param variables List of variables to check (columns in both `model$data` and
#' `ensemble_nwp`)
#' @param indexes List with entries `valid_time` and `base_time`, both should be in
#' `ensemble_nwp`, and `valid_time` in `model$data`.
#' @param lead_time_units passed to `difftime()` to specify units for forecast
#' lead-time in results
#'
#' @description
#' Calculates bias of forecast features (used in forecasting) relative to
#' features used in training, RMSE of ensemble mean, and spread (standard
#' deviation) of ensemble forecasts. Results used in backtesting report.
#'
#' @export
#' @keywords internal
#' @noRd
evaluate_feature_bias <- function(
    model,
    ensemble_nwp,
    variables=NULL,
    indexes=list(valid_time="dtm",
                 base_time="dtm_basetime"),
    lead_time_units = "hours"
){

  if( is.null(variables)) {
    variables <- all.vars(model$formula[-2])
  }

  results <- list()

  nrow_ensemble_nwp_in <- nrow(ensemble_nwp)

  ensemble_nwp <- merge(ensemble_nwp,model$data,by=indexes$valid_time,suffixes = c("","_actual"))
  ensemble_nwp[,lead_time:=difftime(get(indexes$valid_time),get(indexes$base_time),units = lead_time_units)]

  nrow_ensemble_nwp_eval <- nrow(ensemble_nwp)

  results[["nwp_data_fraction"]] <- nrow_ensemble_nwp_eval/nrow_ensemble_nwp_in
  results[["eval_basetime_range"]] <- ensemble_nwp[,range(get(indexes$base_time))]
  results[["n_unique_basetimes"]] <- ensemble_nwp[,length(unique(get(indexes$base_time)))]

  if( results$nwp_data_fraction < 1 ) {
    warning(paste0("Only ", floor(100*results$nwp_data_fraction),
                   "% of ensemble_nwp matched with observations and included in evaluation. There are ",
                   results$n_unique_basetimes, " unique basetimes with range ",
                   results$eval_basetime_range[1]," to ",results$eval_basetime_range[2],
                   " (",difftime(results$eval_basetime_range[2],
                                 results$eval_basetime_range[1],
                                 units = "days")+1,
                   " days)"))
  }

  for( variable in variables ){
    ensemble_nwp[,error := get(paste0(variable,"_actual")) - get(variable)]


    results[[variable]] <- rbind(
      ensemble_nwp[,.(value=mean(error),
                      Metric="Bias"),by=lead_time],
      ensemble_nwp[,.(ens_mean=mean(get(variable)),
                      actual = get(paste0(variable,"_actual"))[1]),by=c(indexes$base_time,"lead_time")][
                        ,.(value=sqrt(mean((ens_mean-actual)^2)),
                           Metric="RMSE"),by=lead_time],
      ensemble_nwp[,stats::sd(get(variable)),by=c(indexes$base_time,"lead_time")][
        ,.(value=mean(V1),
           Metric="Spread"),by=lead_time]
    )

    ensemble_nwp[,error:=NULL]
  }
  return(results)
}


#' Default set of metrics for the assessment of count data models
#'
#' @param p4r_mqgam  A `p4r_mqgam` object that has been fit where `p4r_mqgam$cv=TRUE`.
#'
#' @description function to perform default evaluation, e.g., as part of
#' `fit.p4r_mqgam()`.
#'
#' @details Internal use.
#' @keywords internal
#' @noRd
default_count_evaluation <- function(p4r_mqgam) {

  if(isFALSE(p4r_mqgam$cv)){
    stop("No CV was performed. Cannot compute count metrics evaluation.")
  }

  folds <- p4r_mqgam$cv_results$predictions[!is.na(cv_fold), unique(cv_fold)]

  count_evaluation <- data.table()

  # Addressing each fold to the original sample;
  data <- unique(merge(p4r_mqgam$data,
                       p4r_mqgam$cv_results$predictions[,.SD,.SDcols = c(p4r_mqgam$data_valid_time_index,"cv_fold")],
                       by = p4r_mqgam$data_valid_time_index,all.x = TRUE), by = p4r_mqgam$data_valid_time_index)

  for(fold in folds){

    cat("\n")
    cat(paste0("Iterating over the cv_fold: ", fold ," \n\n"))


    fold_evaluation <- scoringRules_general(y = data[cv_fold == fold,get(p4r_mqgam$formula[[2]])],
                                            newdata = data[cv_fold == fold,],
                                            margin = p4r_mqgam$p4r_tail$tail_distribution,
                                            p4r_mqgam_fold = p4r_mqgam$cv_results$models[[fold]])

    count_evaluation <- rbind(count_evaluation,
                              cbind(data[cv_fold == fold,.SD,.SDcols = c(p4r_mqgam$data_valid_time_index,"cv_fold")],
                                    do.call(cbind,fold_evaluation)))
    flush.console()

  }

  return(count_evaluation)

}


#' Calculate scoring rules for discrete distributions (including DGP)
#'
#' @param A `p4r_mqgam` object that has been fit.
#' @param y the observed response
#' @param newdata Input data to be used to generate new predictions.
#' @param margin distribution used to calculate the distributions. It is available all from the GAMLSS:: package
#' @param params list of parameters adjusted for a given distribution.
#' Not need to specify if is a Discrete Generalised Pareto (DGP) or Generalised
#' Pareto (GP) distribution approximation.
#'
#' @return a list with the calculated logarithm score (`logs`), the quadratic
#' score (`qs`), spherical score (`sphs`), and ranked probability score (`rps`).
#' For more details see (Czado, 2009.).
#' @details Internal use.
#' @keywords internal
#' @noRd
scoringRules_general <- function(y,flexForecast,newdata,max_table,forecast_table = NULL,reference = FALSE,true_pmf = NULL){

  # xx relates to the support of the distribution, 500 is a reasonable choice given the nature of the problem.
  max_global_value <- max(max_table[["max_faults"]])
  if(!is.null(forecast_table)){
    forecast_table <- (merge(forecast_table,max_table,by = "fit_data_row"))[x<=max_faults,,]
  }
  Pdist <- pdens <- numeric(length(y))
  Px <- px <- ind <- matrix(nrow = length(y),ncol = max_global_value+1)
  if(length(y) != nrow(newdata)){
    stop("Incompatible response and dataset.")
  }

  # Params size
  params_size <- nrow(newdata)
  rps_matrix <- matrix(NA,nrow= params_size,ncol = 4)
  colnames(rps_matrix) <- c('logs','qs','sphs','rps')

  count <- 0

  # mu <- get_lambda_poisson(newdata)

  # Computing the case for either fixed shape and scale or single observation
  for(i in 1:params_size){

    # Setting up the progress bar
    count <- count + 1
    sr_progress_bar(current = count,total = params_size,bar_length = 50)

    rps.all <- numeric(length(y))

    if(reference){

        if(is.null(true_pmf)){
            xx <- 0:(max_table[fit_data_row==i,][["max_faults"]])
            pdens[i] <- stats::dpois(x = y[i],lambda = mu[i])
            Px[i,1:(xx+1)] <- stats::ppois(q = xx,lambda = mu[i])
            px[i,1:(xx+1)] <- stats::dpois(x = xx,lambda = mu[i])
        } else {
            if(!is.matrix(true_pmf)){
                xx <- max(as.numeric(names(true_pmf)))
                pdens[i] <- true_pmf[as.character(y[i])]
                Px[i,1:(xx+1)] <- cumsum(true_pmf)
                px[i,1:(xx+1)] <- true_pmf
            } else {
                xx <- max(as.numeric(colnames(true_pmf)))
                pdens[i] <- true_pmf[i,as.character(y[i])]
                Px[i,1:(xx+1)] <- cumsum(true_pmf[i,])
                px[i,1:(xx+1)] <- true_pmf[i,]
            }

        }
    } else {

        if(inherits(flexForecast,"xflex4cast")){

          xx <- 0:(max_table[fit_data_row==i,][["max_faults"]])

          forecast_obs <- if(is.null(forecast_table)){
            forecast(model = flexForecast,
                                                 newdata = newdata[i,],
                                                 indexes = flexForecast$data_index,
                                                 max_quantile_predicted = max(xx))
          } else {
            forecast_table[fit_data_row==i,]
          }


          pdens[i] <- forecast_obs[x == y[i],][["pmf"]]
          Px[i,(xx+1)] <- forecast_obs[["cdf"]]
          px[i,(xx+1)] <- forecast_obs[["pmf"]]

        } else if(is.null(flexForecast$tail_model)){

          forecast_obs <- forecast(flexForecast,
                                   newdata = newdata[i,],
                                   indexes = flexForecast$data_index,
                                   max_quantile_predicted = max(xx))
          pdens[i] <- forecast_obs[x == y[i],][["pmf"]]
          Px[i,(xx+1)] <- forecast_obs[["cdf"]]
          px[i,(xx+1)] <- forecast_obs[["pmf"]]

        }


    }

    ind[i,(xx+1)] <- as.numeric(y[i]<=xx)

  }

  p2 <- apply((px^2),1,sum,na.rm = TRUE)
  ## Not necessary to keep the same length for logs
  # if(any(pdens==0)){
  #   stop("Pdens shouldn't be equal to zero.")
  # }

  logs <- -(log(pdens))

  logs <- ifelse(is.infinite(logs), 10000,logs)
  qs <- (-2*pdens+p2)
  sphs <- (-pdens/sqrt(p2))
  rps <- ( rowSums((Px - ind)^2, na.rm = TRUE) )

  # Returning the matrix
  return(list(logs = logs, qs = qs, sphs = sphs, rps = rps))


}


#' Calculate scoring rules for discrete distributions (including DGP)
#'
#' @param A `p4r_mqgam` object that has been fit.
#' @param y the observed response
#' @param newdata Input data to be used to generate new predictions.
#' @param margin distribution used to calculate the distributions. It is available all from the GAMLSS:: package
#' @param params list of parameters adjusted for a given distribution.
#' Not need to specify if is a Discrete Generalised Pareto (DGP) or Generalised
#' Pareto (GP) distribution approximation.
#'
#' @return a list with the calculated logarithm score (`logs`), the quadratic
#' score (`qs`), spherical score (`sphs`), and ranked probability score (`rps`).
#' For more details see (Czado, 2009.).
#' @details Internal use.
#'
#' @keywords internal
#' @noRd
get_scoringRules_ss_from_eval <- function(model_eval, reference_eval,perfect_eval){

  ss_list <- list(ss_logs = (reference_eval[["logs"]]-model_eval[["logs"]])/(reference_eval[["logs"]]-perfect_eval[["logs"]]),
                   ss_qs = (reference_eval[["qs"]]-model_eval[["qs"]])/(reference_eval[["qs"]]-perfect_eval[["qs"]]),
                   ss_sphs =  (reference_eval[["sphs"]]-model_eval[["sphs"]])/(reference_eval[["sphs"]]-perfect_eval[["sphs"]]),
                   ss_rps = (reference_eval[["rps"]]-model_eval[["rps"]])/(reference_eval[["rps"]]-perfect_eval[["rps"]]))

  # Returning the list
  return(ss_list)

}

#' @title Compute Scoring Rules-Based Skill Scores
#' @description
#' Internal function to compute skill scores using scoring rules (logarithmic score, quadratic score, spherical score, and ranked probability score) for count data.
#'
#' @details
#' The function compares a model and a reference forecast against the true data-generating process using various proper scoring rules. The function assumes Poisson distributions for now and operates on a fixed grid (support) for count outcomes.
#'
#' The skill scores are computed as:
#' \deqn{SS = \frac{Score_{ref} - Score_{model}}{Score_{ref} - Score_{true}}}
#' for each metric.
#'
#' @param y Numeric vector of observed counts.
#' @param model A fitted model of class `xflex4cast`.
#' @param reference A reference model of class `xflex4cast`.
#' @param newdata A `data.table` or `data.frame` with the covariates.
#' @param simulation_distribution Distribution assumed for simulation. Default is "poisson".
#'
#' @return A list with skill scores:
#' \describe{
#'   \item{ss_logs}{Logarithmic score skill score}
#'   \item{ss_qs}{Quadratic score skill score}
#'   \item{ss_sphs}{Spherical score skill score}
#'   \item{ss_rps}{Ranked probability score skill score}
#' }
#'
#' @keywords internal
#' @noRd

scoringRules_ss <- function(y, model, reference,newdata, simulation_distribution = "poisson"){

  # xx relates to the support of the distribution, 500 is a reasonable choice given the nature of the problem.
  xx <- 0:3000
  n <- length(y)
  pdens_true <-  pdens_model <- pdens_reference <- numeric(length(y))
  Px_true <- px_true <- Px_reference <- px_reference <- Px_model <- px_model <- ind <- matrix(nrow = length(y),ncol = length(xx))

  if(length(y) != nrow(newdata)){
    stop("Incompatible response and dataset.")
  }

  # Params size
  params_size <- nrow(newdata)

  count <- 0
  if(simulation_distribution == "poisson"){
    mu <- get_lambda_poisson(newdata)
  }

  # Creating the list to store all ements
  p2 <- logs <- qs <- sphs <- rps <- list(model = numeric(n),reference  = numeric(n), true = numeric(n))

  # Computing the case for either fixed shape and scale or single observation
  for(i in 1:params_size){

    # Setting up the progress bar
    count <- count + 1
    sr_progress_bar(current = count,total = params_size,bar_length = 50)

    rps.all <- numeric(length(y))

    if(simulation_distribution == "poisson"){

      pdens_true[i] <- stats::dpois(x = y[i],lambda = mu[i])
      Px_true[i,] <- stats::ppois(q = xx,lambda = mu[i])
      px_true[i,] <- stats::dpois(x = xx,lambda = mu[i])

    }

    if(inherits(model,"xflex4cast")){

      forecast_obs_model <- forecast(model = model,
                                                 newdata = newdata[i,],
                                                 indexes = model$data_index,
                                                 max_quantile_predicted = max(xx))


      pdens_model[i] <- forecast_obs_model[x == y[i],][["pmf"]]
      Px_model[i,] <- forecast_obs_model[["cdf"]]
      px_model[i,] <- forecast_obs_model[["pmf"]]

    }

    if(inherits(model,"xflex4cast")) {


      forecast_obs_reference <- forecast(model = reference,
                                                     newdata = newdata[i,],
                                                     indexes = reference$data_index,
                                                     max_quantile_predicted = max(xx))


      pdens_reference[i] <- forecast_obs_reference[x == y[i],][["pmf"]]
      Px_reference[i,] <- forecast_obs_reference[["cdf"]]
      px_reference[i,] <- forecast_obs_reference[["pmf"]]

    }



    ind[i,] <- as.numeric(y[i]<=xx)

  }

  # Calculating all metrics

  p2[["true"]] <- apply((px_true^2),1,sum)
  p2[["model"]] <- apply((px_model^2),1,sum)
  p2[["reference"]] <- apply((px_reference^2),1,sum)

  logs[["true"]] <- -(log(pdens_true))
  logs[["model"]] <- -(log(pdens_model))
  logs[["reference"]] <- -(log(pdens_reference))


  logs$true <- ifelse(is.infinite(logs$true),10000,logs$true)
  logs$model <- ifelse(is.infinite(logs$model),10000,logs$model)
  logs$reference <- ifelse(is.infinite(logs$reference),10000,logs$reference)


  qs[["true"]] <- (-2*pdens_true+p2[["true"]])
  qs[["model"]] <- (-2*pdens_model+p2[["model"]])
  qs[["reference"]] <- (-2*pdens_reference+p2[["reference"]])

  sphs[["true"]] <- (-pdens_true/sqrt(p2[["true"]]))
  sphs[["model"]] <- (-pdens_model/sqrt(p2[["model"]]))
  sphs[["reference"]] <- (-pdens_reference/sqrt(p2[["reference"]]))

  rps[["true"]] <- ( rowSums((Px_true - ind)^2) )
  rps[["model"]] <- ( rowSums((Px_model - ind)^2) )
  rps[["reference"]] <- ( rowSums((Px_reference - ind)^2) )

  # Returning the matrix
  return(list(ss_logs = (logs[["reference"]]-logs[["model"]])/(logs[["reference"]]-logs[["true"]]),
              ss_qs = (qs[["reference"]]-qs[["model"]])/(qs[["reference"]]-qs[["true"]]),
              ss_sphs =  (sphs[["reference"]]-sphs[["model"]])/(sphs[["reference"]]-sphs[["true"]]),
              ss_rps = (rps[["reference"]]-rps[["model"]])/(rps[["reference"]]-rps[["true"]])) )


}

#' @title Compute Simple Skill Score
#' @description
#' Internal function to compute a simple skill score (`ss_score`) comparing model performance against a reference and a perfect model.
#'
#' @details
#' The score is calculated using the formula:
#' \deqn{SS = \frac{ref - eval}{ref - perf}}
#' where `eval` is the model being evaluated, `ref` is the reference model, and `perf` is the perfect model.
#' A higher skill score indicates better performance relative to the reference.
#'
#' @param eval_model Numeric vector or value representing performance of the evaluated model.
#' @param ref_model Numeric vector or value representing performance of the reference model.
#' @param perf_model Numeric vector or value representing performance of the perfect model.
#'
#' @return A numeric value or vector containing the computed skill score.
#'
#' @examples
#' get_simple_skill_score(eval_model = 0.2, ref_model = 0.5, perf_model = 0.0)
#'
#' @keywords internal
#' @noRd

get_simple_skill_score <- function(eval_model,
                                   ref_model,
                                   perf_model){

  return((ref_model-eval_model)/(ref_model-perf_model))

}

#' Calculate the bootstrap means for the metrics used for the predictive model
#' assessment for count data
#'
#' @param p4r_mqgam a fitted `p4r_mqgam` model.
#' @param n_replications number of bootstrap means to be calculated.
#' @param metrics a character vector to select which metrics will be calculated.
#' The possible options are: logarithm score (`logs`), the quadratic
#' score (`qs`), spherical score (`sphs`), and ranked probability score (`rps`).
#' Default is `c("logs","qs","sphs","rps")`. For more details see (Czado, 2009.).
#'
#' @return a numerical matrix where each column correspond to each metric, and
#' each row contains the mean corresponding to a bootstraping sample.
#'
#' @details Internal use.
#' @keywords internal
#' @noRd
bootstrap_metrics <- function(p4r_mqgam,
                              n_replications = 1000,
                              metrics = c("logs","qs","sphs","rps")){


  sample_size <- nrow(p4r_mqgam$cv_results$evaluation$count_eval)

  bootstrap_samples <- replicate(n = n_replications,
                                 sample(1:sample_size,size = sample_size,replace = TRUE))

  bootstrap_means <- matrix(NA,nrow = n_replications, ncol = length(metrics))
  colnames(bootstrap_means) <- metrics

  for(i in 1:n_replications){
    bootstrap_means[i,] <- colMeans(p4r_mqgam$cv_results$evaluation$count_eval[bootstrap_samples[i,],.SD,.SDcols = metrics])
  }

  return(bootstrap_means)
}
