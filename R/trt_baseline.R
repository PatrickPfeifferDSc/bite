#' Simulated data on baseline subject level
#'
#' A simulated dataset (from the framework of the shared factor model,
#' see ?select_trt_sharedfac) to
#' demonstrate the use of the SR-model. The file contains
#' 5000 subjects with their respective feature information at baseline time.
#' Consider these features to be relevant for confounding and influencing the
#' treatment intake, not the outcomes.
#' First of two data tables given to function \code{\link{bayesTrtEffects}}.
#'
#' @docType data
#'
#' @format A data frame with subject related feature columns:
#' \describe{
#'    \item{ID: the column which uniquely identifies each subject}
#'    \item{Ti: panel times, the number of subsequent meassures for a subject}
#'    \item{trt: treatment, usually coded 0 or 1}
#'    \item{V1: variable 1, normally distributed}
#'    \item{V2: variable 2, 0/1 distributed}
#'    \item{V3: variable 3, 0/1 distributed}
#' }
#' @source This dataset stems from a simulation process and represents fictive data.

"trt_baseline"
