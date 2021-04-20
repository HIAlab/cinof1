##############################################################
# G-Estimation
# 4/4/2021
# created by T. Gärtner (thomas.gaertner@student.hpi.de)
##############################################################
#
# This file contains functions to analyze the treatment effect for each patient.
#



#' Fit LM with Adjusted Treatment Effect
#' @description This function fits a model with a given treatment effect.
#' @param data Data Frame
#' @param exposure
#' @param outcome
#' @param variables
#' @param effects
#' @return fitted linear model
#' @example
#' outcome <- "Uncertain_Low_Back_Pain"
#' exposure <- "treatment"
#' variables <- c("Activity")
#' id <- "patient_id"
#' time_col <- "day"
#'
#' result <- estimate_gamma_tau(data = simpatdat, outcome = outcome, exposure = exposure, variables = variables, bound = 3, symmetric = TRUE, id=id, time_col = time_col)
#' fit.adj.lm(data = simpatdat, outcome = outcome, exposure = exposure, variables = variables, effects = result$best, id = id, time_col = time_col)
#' @export

fit.adj.lm <- function(data, outcome, exposure, variables, effects, id, time_col, one.hot=FALSE){
  if(one.hot==FALSE){
    res <- prep.onehot(data, exposure)
    data <- res$data
    exposure.columns <- res$names
  }else{
    exposure.columns <- exposure
  }


  for (exposure.column in exposure.columns){
    gamma <- effects[[paste(exposure.column, "gamma", sep = ".")]]
    tau <- effects[[paste(exposure.column, "tau", sep = ".")]]
    data[,exposure.column] <- gen_treatment_effect(data, exposure=exposure.column, gamma=gamma, tau=tau, id=id, time_col = time_col)
  }

  str_formula <- sprintf("%s ~ %s", outcome, paste(exposure.columns, variables, sep=" + ", collapse = " + "))
  lin_m <- lm(formula(str_formula), data = data, na.action = na.omit)
  return(lin_m)
}


#' Generate Treatment Effect
#' @description  This function generates an expected treatment effect for each time point based on gamma and tau.
#' @param treatment vector of 1 (treated) and 0 (not treated)
#' @param gamma wash-in effect (default 1)
#' @param tau wash-out effect (default 1)
#' @return vector with expected treatment effect
#' @examples
#' # Generate dummy data
#' treatment <- rep(c(rep(0,7), rep(1,7)),4)
#' # run example
#' gamma <- 5
#' tau <- 7
#' gen_treatment_effect(treatment, gamma, tau)

gen_treatment_effect <- function(data, exposure, gamma, tau, id, time_col) {

  est_effect <- rep(NA, nrow(data))
  for(v in c(1:length(est_effect))){

    data_point <- data[data[,id] == data[v,id] & data[,time_col] == data[v,time_col]-1,]
    if (!nrow(data_point)==0){
      x_j <- data_point[1,exposure]
    }else{
      x_j <- 0
    }

    est_effect[v] = (x_j + ((1 - x_j) / tau) * data[v,exposure] - (x_j / gamma) * (1 - data[v,exposure]))
  }
  return(est_effect)
}

#' Estimate Gamma and Tau
#' @description This function analyzes a patient with time varying treatment.
#' It analyzes the patient with several wash-in and wash-out effects and returns the best gamma and tau. It is implemented for 2 different treatments.
#' @param data a data frame holding treatment variables and outcome.
#' @param outcome defines the outcome column
#' @param exposure defines the exposure variable (level)
#' @param bound maximum value for gamma and tau
#' @param symmetric identifies if wash-in effect is equal to wash-out effect
#' @return a vector with best gamma and tau.
#' @examples
#' # Define Variables
#' outcome <- "Uncertain_Low_Back_Pain"
#' exposure <- "treatment"
#' variables <- c("Activity")
#' id <- "patient_id"
#' time_col <- "day"
#' # Estimate best gamma and tau
#' result <- estimate_gamma_tau(data = simpatdat, outcome = outcome, exposure = exposure, variables = variables, bound = 3, symmetric = TRUE, id=id, time_col = time_col)
#' # Fit Linear Model
#' fit.adj.lm(data = simpatdat, outcome = outcome, exposure = exposure, variables = variables, effects = result$best, id = id, time_col = time_col)
#'
#' @export
estimate_gamma_tau <- function(data, outcome, exposure, variables, bound=10, symmetric=TRUE, id="id", time_col="day", one.hot=FALSE) {

  # Prepare One Hot encodiing
  if(one.hot==FALSE){
    res <- prep.onehot(data, exposure)
    data <- res$data
    exposure.names <- res$names
  }else{
    exposure <- NA
    exposure.names <- exposure
  }



  grid <- list()

  # define array for later evaluation
  if(symmetric){
    # iterate over all gamma values
    for(i in c(1:length(exposure.names))){
      grid[[paste(exposure.names[i], "gamma", sep=".")]] <- c(1:bound)
    }

    #expand
    grid <- expand.grid(grid)

    # copy values for tau to have symmetric
    for(i in c(1:length(exposure.names))){
      grid[,paste(exposure.names[i], "tau", sep=".")] <- grid[,paste(exposure.names[i], "gamma", sep=".")]
    }


  }else{
    for(i in c(1:length(exposure.names))){
      grid[[paste(exposure.names[i], "gamma", sep=".")]] <- c(1:bound)
      grid[[paste(exposure.names[i], "tau", sep=".")]] <- c(1:bound)
    }
    #expand
    grid <- expand.grid(grid)
  }

  r2 <- c()

  #iterate over values
  for(row_id in c(1:nrow(grid))){
    test.effects = grid[row_id,]
    res <- summary(fit.adj.lm(data, outcome, exposure.names, variables, effects=test.effects, id=id, time_col = time_col, one.hot=TRUE))
    r2[row_id] <- res$adj.r.squared
  }

  # get the best model by the biggest r2 value
  best <- which(r2 == max(r2))
  best.values <- grid[best,]

  return(list(data = cbind(grid, data.frame(r2)), best=best.values))
}


