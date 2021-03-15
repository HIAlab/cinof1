#' Simulated N-of-1 study data with 300 patients
#'
#' The dataset containing exposure and outcome of around 300 patients over time.
#'
#' @format A data frame with 42000 rows and 12 variables:
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


