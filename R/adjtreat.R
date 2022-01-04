##############################################################
# G-Estimation
# last edit 6/5/2021
# created by T. Gärtner (thomas.gaertner@student.hpi.de)
##############################################################
#
# This file contains functions to analyze the treatment effect for each patient.
#



#' Fit LM with Adjusted Treatment Effect
#' @description This function fits a linear model with a given treatment effect.
#' For that, it pre processes the data with the given effect and returns the
#' fitted linear model. An unadjusted linear model can be generated with setting
#' `effects` to `NA`.
#' @param data Data Frame with patient data
#' @param exposure name of the exposure variable column
#' @param outcome name of the outcome variable column
#' @param variables vector of names of other variables
#' @param effects list of effects for
#' @param one.hot boolean weather the exporure
#' variable is on hot encoded or not-
#' @return fitted linear model
#' @example
#' # Define Variables
#' outcome <- "Uncertain_Low_Back_Pain"
#' exposure <- "treatment"
#' variables <- c("Activity")
#' id <- "patient_id"
#' time_col <- "day"
#'
#' # use the estimate.gamma.tau function to estimate the best gamma tau values
#' result <- estimate.gamma.tau(data = simpatdat, outcome = outcome, exposure = exposure, variables = variables, bound = 3, symmetric = TRUE, id=id, time_col = time_col)
#' # fit the adjusted linear model
#' fit.adj.lm(data = simpatdat, outcome = outcome, exposure = exposure, variables = variables, id = id, time_col = time_col, effects = result$best)
#' @export

fit.adj.lm <- function(data, outcome, exposure, variables, id, time_col, effects=NA, one.hot=FALSE){
  # Prepare One Hot encoding
  if(one.hot==FALSE){
    res <- prep.onehot(data, exposure)
    data <- res$data
    exposure.columns <- res$names
  }else{
    exposure.columns <- exposure
  }

  # If no carry-over effects should be adjusted
  if(is.na(effects)){
    # generate string formula
    str_formula <- sprintf("%s ~ %s", outcome, paste(exposure.columns, variables, sep=" + ", collapse = " + "))
    #fit model
    lin_m <- lm(formula(str_formula), data = data, na.action = na.omit)
    # return model
    return(lin_m)
  }else{

    # empty list of adjusted exposure column name
    adj.exposure.columns <- c()

    # iterate over exposure columns
    # generate for each column the adjusted effect and finally add the name
    # to adj.exposure.columns
    for(i in c(1:length(exposure.columns))){
      exposure.column <- exposure.columns[i]
      # get gamma tau values from the effect list
      gamma <- effects[[paste(exposure.column, "gamma", sep = ".")]]
      tau <- effects[[paste(exposure.column, "tau", sep = ".")]]
      # prepare data
      data <- gen.treatment.effect.col(data,
                                       exposure=exposure.column,
                                       gamma=gamma,
                                       tau=tau,
                                       id=id,
                                       time_col = time_col)
      adj.exposure.columns[i] <- paste(exposure.column,
                                       "gamma",
                                       gamma,
                                       "tau",
                                       tau,
                                       sep=".")
    }

    # generate formula
    str_formula <- sprintf("%s ~ %s",
                           outcome,
                           paste(
                             adj.exposure.columns,
                             variables,
                             sep=" + ",
                             collapse = " + "))

    # fit model
    lin_m <- lm(formula(str_formula), data = data, na.action = na.omit)
    # return fitted model
    return(lin_m)
  }
}

#' Generate Treatment Effect
#' Calculates the treatment effect for given values
est.effect <- function(x_i, exposure_j, tau, gamma){
  return(x_i + ((1 - x_i) / tau) * exposure_j - (x_i / gamma) * (1 - exposure_j))
}



#' Generate Treatment Effect Column (not exported)
#'
#' @description  This function generates an expected treatment effect for each time point based on gamma and tau value.
#' @param data data frame with patients
#' @param exposure vector of 1 (treated) and 0 (not treated)
#' @param gamma wash-in effect
#' @param tau wash-out effect
#' @param id defines the unique identifier for each patient
#' @param time_col column name, which defines the time
#' @return vector with expected treatment effect
#' @examples
#' # Generate dummy data
#' treatment <- rep(c(rep(0,7), rep(1,7)),4)
#' # run example
#' gamma <- 5
#' tau <- 7
#' gen.treatment.effect(treatment, gamma, tau)
gen.treatment.effect.col <- function(data, exposure, gamma, tau, id, time_col) {

  unique.pat.ids <- unique(data[,id])

  res.data <- foreach::foreach(i = 1:length(unique.pat.ids) , .combine ="rbind") %do% {
    patient.id <- unique.pat.ids[i]
    pat.data <- data[data[,id]==patient.id,]

    pat.data[,paste(exposure,"gamma",gamma,"tau",tau, sep=".")] <- rep(NA, nrow(pat.data))

    pat.data <- pat.data[order(pat.data[,time_col]),]

    for(v in c(1:nrow(pat.data))){

      e_j <- pat.data[v,exposure]
      if (is.na(e_j)){
        res <- NA
      }else{
        if (v>1){
          x_i <- pat.data[v-1,paste(exposure,"gamma",gamma,"tau",tau, sep=".")]
          if (is.na(x_i)){
            x_i <- e_j
          }
        }else{
          x_i <- 0
        }
        res <- est.effect(x_i, e_j, tau, gamma)
      }
      pat.data[v,paste(exposure,"gamma",gamma,"tau",tau, sep=".")] <- res
    }
    pat.data
  }

  return(res.data)
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
#'
#' # use the estimate.gamma.tau function to estimate the best gamma tau values
#' result <- estimate.gamma.tau(data = simpatdat, outcome = outcome, exposure = exposure, variables = variables, bound = 3, symmetric = TRUE, id=id, time_col = time_col)
#' # fit the adjusted linear model
#' fit.adj.lm(data = simpatdat, outcome = outcome, exposure = exposure, variables = variables, id = id, time_col = time_col, effects = result$best)
#' @export
estimate.gamma.tau <- function(data, outcome, exposure, variables, bound=10, symmetric=TRUE, id="id", time_col="day", one.hot=FALSE) {

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


  #iterate over values
  result <- foreach::foreach(row_id = c(1:nrow(grid)), .combine = rbind) %do% {
    test.effects <- grid[row_id,]
    res <- summary(fit.adj.lm(data, outcome, exposure.names, variables, effects=test.effects, id=id, time_col = time_col, one.hot=TRUE))


    # create return values
    res_vec <- list()
    res_vec$r2 <- res$adj.r.squared

    for(exposure.column in exposure.names){
      gamma <- test.effects[[paste(exposure.column, "gamma", sep = ".")]]
      tau <- test.effects[[paste(exposure.column, "tau", sep = ".")]]
      exp <- paste(exposure.column,"gamma",gamma,"tau",tau, sep=".")

      if(exp %in% row.names(res$coefficients)){
        res_vec[[paste(exposure.column, "Estimate")]] <- res$coefficients[exp,"Estimate"]
        res_vec[[paste(exposure.column, "Std. Error")]] <- res$coefficients[exp,"Std. Error"]
      }else{
        res_vec[[paste(exposure.column, "Estimate")]] <- NA
        res_vec[[paste(exposure.column, "Std. Error")]] <-NA
      }
    }
    data.frame(res_vec)
  }

  r2 <- result$r2
  # get the best model by the biggest r2 value
  best <- which(r2 == max(r2))
  best.values <- grid[best,]

  return(list(data = cbind(grid, result), best=best.values))
}


