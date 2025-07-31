################################################################################
# Simulation Code: Section 3.1.1 – First Scenario: Scale-Invariance in the Tail
#
# Description:
# ------------
# This simulation evaluates the predictive performance of three approaches:
#     - Standard QGAM
#     - X-flexForecast (with tail modeling)
#     - Extended Discrete Generalized Pareto (EDGP) (https://doi.org/10.1177/1471082X24126672)
#
# The simulation assumes:
#     - A discrete Gamma bulk distribution with predictor-dependent scale.
#     - Tail observations following a Discrete Generalized Pareto with shape ξ = 0.3, and \sigma = 2.5
#     - φ = 0.1 as the threshold exceedance proportion.
#
# Setup:
#  - Wind-speed derived covariates (x1) are sampled empirically.
#  - Penalized splines (k = 10) are used for modeling conditional quantiles.
#
# The goal is to compare the RMSE of predicted quantiles across replications.
#
# Note:
#  - The number of replications here is reduced (n = 10) for illustration.
#  - The original manuscript uses 1000 replications.
#
# Output:
#  - A summary of average RMSE ± SD across replications for each method.
#  - Output is saved as RDS for reproducibility and visualization.
#
#

rm(list=ls())

# Installing only if necessary
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

library(devtools)
set.seed(42)

# Loading all package dependencies -- you can also install the xflex4cast package instead,
# and load it : library(xflex4cast)
devtools::load_all()
load_degpd() # Loading all functions necessary to run the competing model DEGPD.
             # For more details: (Ahmad, T., Gaetan, C., & Naveau, P. (2024).
             # An extended generalized Pareto regression model for count data.
             # *Statistical Modelling*, 0(0).
             #https://doi.org/10.1177/1471082X241266729)

# Setting the parameters of the simulation
n <- 5000
phi <- 0.1
n_replications <- 10 # For illustration purposes we reduced the number of replications to 100, while in the main manuscript is 1000
xi <- 0.3 # Shape parameter
probs <- c(0.05,0.25,0.5,0.75,0.9,0.95,0.99,0.999,0.9999) # Quantile levels to calculate on

# Iterating over the number of replications
store_replications_results <- vector("list",length = n_replications)
names(store_replications_results) <- paste0("replication",1:n_replications)


for( i in 1:n_replications ){
  train_data <- simulation_scenario_one(n = n,
                                        phi = phi,
                                        xi = xi,
                                        probs = probs,
                                        index = "fit_data_row")

  test_data <- simulation_scenario_one(n = n,
                                       phi = phi,
                                       xi = xi,
                                       probs = probs,
                                       index = "fit_data_row")


# Defining the qgam model only
std_qgam_def <- define_xflex4cast(mqgam_formula = y ~ s(x1,k=10) ,
                                    mqgam_quantiles = probs,
                                    quantile_threshold = (1-phi),
                                    fit_tail = FALSE)

# Using the flex4cast approach
xflex4cast_def <- define_xflex4cast(mqgam_formula = y ~ s(x1,k=10) ,
                                          scale_formula = ~ 1,
                                          mqgam_quantiles = probs,
                                          quantile_threshold = (1-phi),
                                          fit_tail = TRUE)

# Defining the formula for the DEGP approach
formula_degp <- list(lsigma = y ~ s(x1,k=10), lxi = ~1, lkappa = ~1)


# Fitting the models
qgam_fit <- fit(model = std_qgam_def,
                data = train_data$data,
                index = "fit_data_row")

xflex4cast_fit <- fit(model = xflex4cast_def,
                data = train_data$data,
                index = "fit_data_row")

degpd_fit <- devgam(formula_degp,data = as.data.frame(train_data$data),
                    family  = "degpd",
                    degpd.args = list(m=1),
                    trace = 0)

# Obtaining the predictions
oos_pred_qgam <- predict(qgam_fit,
                         newdata = test_data$data[,fit_data_row:=1:.N],
                         indexes = "fit_data_row")

oos_pred_xflex4cast <- predict(xflex4cast_fit,
                               newdata = test_data$data[,fit_data_row:=1:.N],
                               indexes = "fit_data_row")

# Obtaining out-of-sample observations
oos_pred_degpd <- as.data.table(predict(degpd_fit,as.data.frame(test_data$data),prob = probs))[,fit_data_row := 1:.N]
oos_pred_degpd <- data.table::melt(oos_pred_degpd,id.vars = "fit_data_row",variable.name = "quantile",value.name = "value")
oos_pred_degpd[, quantile := as.numeric(as.character(quantile))]


# Creating the table for the true quantiles in a long format
if(is.vector(test_data$true_quantile)){
  true_qt_matrix <- test_data$true_quantiles
  colnames(true_qt_matrix) <- probs
  true_quantile <- data.table::melt(as.data.table(true_qt_matrix,keep.rownames = TRUE)[,fit_data_row:=.I],
                        id.vars = "fit_data_row", variable.name = "quantile",value.name = "value")
} else {
  true_quantile <-  data.table::melt(as.data.table(test_data$true_quantiles,keep.rownames = TRUE)[,fit_data_row:=.I],
                         id.vars = "fit_data_row", variable.name = "quantile",value.name = "value")
}

true_quantile <- true_quantile[,value:=ifelse(is.na(value),1000,value)]
true_quantile <- true_quantile[,quantile:=as.numeric(as.character(quantile))]


all_pred <- merge(true_quantile,oos_pred_qgam, by = c("fit_data_row","quantile"),suffixes = c(".true",".qgam"))
all_pred <- merge(all_pred, oos_pred_xflex4cast, by = c("fit_data_row","quantile"))
all_pred <- setnames(all_pred,old = "value","value.tail")
all_pred <- merge(all_pred,oos_pred_degpd, by = c("fit_data_row","quantile"))
all_pred <- setnames(all_pred,old = "value","value.degpd")

# Generating the new columns (doing the same for the forecast)
rmse_summary <- all_pred[,.(rmse_qgam = sqrt(mean((value.true-round(value.qgam))^2)),
                            rmsexflex4cast = sqrt(mean((value.true-round(value.tail))^2)),
                            rmse_degpd = sqrt(mean((value.true-round(value.degpd))^2)),
                            mae_qgam = mean(abs(value.true-round(value.qgam))),
                            maexflex4cast = mean(abs(value.true-round(value.tail))),
                            mae_degpd = mean(abs(value.true-round(value.degpd)))),
                         by = "quantile"]
store_replications_results[[i]] <- rmse_summary[,replication_number:=i]

# Progress bar
progress <- round((i / n_replications) * 100)
bar <- paste0(rep("=", progress / 2), collapse = "")
cat(sprintf("\rProgress: [%-50s] %3d%%", bar, progress))
flush.console()

}

# Combine all data.tables into one
combined <- data.table::rbindlist(store_replications_results)

# Computing mean and sd by quantile
summary_rmse <- combined[, .(
  mean_rmse_qgam  = mean(rmse_qgam),
  sd_rmse_qgam    = stats::sd(rmse_qgam),
  mean_rmsexflex4cast  = mean(rmsexflex4cast),
  sd_rmsexflex4cast    = stats::sd(rmsexflex4cast),
  mean_rmse_degpd = mean(rmse_degpd),
  sd_rmse_degpd   = stats::sd(rmse_degpd)
), by = quantile]

# Visualizing the results as in the paper
summary_formatted <- summary_rmse[, .(
  quantile,
  rmsexflex4cast  = sprintf("%.2f ± %.2f", mean_rmsexflex4cast, sd_rmsexflex4cast),
  rmse_qgam  = sprintf("%.2f ± %.2f", mean_rmse_qgam, sd_rmse_qgam),
  rmse_degpd = sprintf("%.2f ± %.2f", mean_rmse_degpd, sd_rmse_degpd)
)]

print(summary_formatted)


# saveRDS(summary_formatted,
#         file = paste0("inst/simulations/summary_formatted_",n_replications,"_replications_simulatio_3_1_1.rds"))
