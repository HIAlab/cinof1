#' Simulated N-of-1 study data with 300 patients
#'
#' The data set containing exposure and outcome of around 300 patients over time.
#'
#' @format A data frame with 33600 rows and 7 variables:
#' \describe{
#'   \item{X}{row id}
#'   \item{date}{date of observation}
#'   \item{day}{day in study}
#'   \item{Activity}{Steps per day}
#'   \item{Treatment_1}{Identifies whether the patient gets Treatment 1 or not}
#'   \item{Uncertain_Low_back_Pain}{Observation of uncertain low back pain (outcome) on a scale of 0 to 15}
#' }
#' @source Simulate Data with \url{https://github.com/thogaertner/n-of-1-simulation}
"simpatdat"


#' Dag
#'
#' This is an example dag for the data set `simpatdat`. It contains an _exposure_ variable (`treatment`), the _outcome_ variable (`Uncertain_Low_Bacl_Pain`) and one confounding variable (`Activity`).
#'
#' @format A dagitty string with an example dag for the `simpatdat` data set.
#' \describe{
#'   \item{Dag}{dag string}
#' }
#' @source This dag represents the causal effects in the example data set `simpatdat`
"simpatdag"
