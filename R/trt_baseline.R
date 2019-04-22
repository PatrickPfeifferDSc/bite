#' Example dataset to demonstrate use and structure of bate package
#'
#' A simulated dataset (from the framework of the shared factor model) to
#' demonstrate the use of the SR model in package 'bite'. The file contains
#' 5000 subjects with their respective feature information at baseline time.
#' First of two data tables given to \code{\link{bayesTrtEffects}}.
#'
#' @format A list containing numerous objects:
#' \itemize{
#'    \item{ID: the column which uniquely identifies each subject}
#'    \item{Ti: panel times, the number of subsequent meassures for a subject}
#'    \item{trt: treatment, usually coded 0 or 1}
#'    \item{V1: variable 1, for this dataset normally distributed}
#'    \item{V2: variable 2, for this dataset binary}
#'    \item{V3: variable 3, for this dataset binary}
#' }
#'
#' @source This dataset stems from a simulation process and represents fictive data.

"trt_baseline"
