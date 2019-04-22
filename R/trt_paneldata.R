#' Example dataset containing the panel information
#'
#' A simulated dataset (from the framework of the shared factor model) to
#' demonstrate the use of the SR model in package 'bite'. The file contains
#' 5000 subjects, each of which has 4 panel observations on simulated variables.
#' The second of 2 argument datasets for\code{\link{bayesTrtEffects}}.
#'
#' @format A list containing numerous objects:
#' \describe{
#' \item{ID}{ID of subject/item}
#' \item{panelT}{timepoint t of meassured features and dependent variable}
#' \item{y}{dependent variable, outcome}
#' \item{V1}{variable 1, same as V1 in baseline}
#' \item{V2}{variable 2, same as V2 in baseline}
#' \item{t2}{dummy for panel time = 2}
#' \item{t3}{dummy for panel time = 3}
#' \item{t4}{dummy for panel time = 4}
#' }
#'
#' @source This dataset stems from a simulation process and represents fictive data.

"trt_paneldata"
