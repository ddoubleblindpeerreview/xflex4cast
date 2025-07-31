#' @importFrom data.table data.table :=
#' @importFrom stats predict quantile sigma
#' @importFrom utils data flush.console
NULL

#' Load DEGPD and ZIDEGPD R and C++ Functions
#'
#' Dynamically loads all R and C++ source files required for the DEGPD and ZIDEGPD models.
#'
#' This function sources all R files from `inst/degpd-and-zidegpd/R/` and attempts to
#' compile and load all C++ files via `Rcpp::sourceCpp()` from `inst/degpd-and-zidegpd/src/`.
#' It is designed to be used during package development or advanced workflows where
#' direct loading of source files is needed.
#'
#' @details
#' The function assumes the existence of two subdirectories within the `inst/degpd-and-zidegpd/` path:
#' - `R/`: containing R source files to be loaded via `source()`.
#' - `src/`: containing C++ source files to be compiled and loaded via `Rcpp::sourceCpp()`.
#'
#' If any C++ file fails to compile, the function will suppress the error using `try()`, and
#' continue executing. This may result in incomplete loading of required functionality.
#'
#' This function is intended for internal use, particularly useful in development environments. For more details
#' check the original repository: https://github.com/touqeerahmadunipd/degpd-and-zidegpd; and mail the original
#' author at touqeer.ahmad8960@gmail.com.
#'
#' @return Invisibly returns `NULL`. Used for its side effects.
#'
#' @seealso [Rcpp::sourceCpp()], [base::source()], [list.files()]
#'
#' @importFrom Rcpp sourceCpp
#' @references
#' Ahmad, T., Gaetan, C., & Naveau, P. (2024). An extended generalized Pareto regression model for count data. *Statistical Modelling*, 0(0). https://doi.org/10.1177/1471082X241266729
#'
#' @examples
#' \dontrun{
#' # Load all DEGPD and ZIDEGPD functions (R and C++)
#' load_degpd()
#' }
#' @export
#'
load_degpd <- function(){
  # Loading all files for the DEGPD fit
  dest <- "inst/degpd-and-zidegpd/R/"      # this function all R function
  files = list.files(dest, full.names = T)
  for (i in 1:length(files)) {
    source(files[i])
  }

  dest <- "inst/degpd-and-zidegpd/src/"  # This function will call all C++ files
  files = list.files(dest, full.names = T)
  try({
    for (i in 1:length(files)) {
      Rcpp::sourceCpp(files[i])
    }
  })
  return(invisible(NULL))
}
