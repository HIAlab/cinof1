##############################################################
# G-Estimation
# 4/4/2021
# created by T. Gärtner (thomas.gaertner@student.hpi.de)
##############################################################
#
# This file contains functions to analyze the treatment effect for each patient.
#

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

gen_treatment_effect <- function(treatment, gamma=1, tau=1) {
  est_effect <- rep(NA, length(treatment))
  for(v in c(1:length(treatment))){
    if (v==1){
      x_j<-0
    }else{
      x_j<-est_effect[v-1]
    }
    est_effect[v] = (x_j + ((1 - x_j) / tau) * treatment[v] - (x_j / gamma) * (1 - treatment[v]))
  }
  return(est_effect)
}

#' Estimate Gamma and Tau
#' @description This function analyzes a patient with time varying treatment.
#' It analyzes the patient with several wash-in and wash-out effects and returns the best gamma and tau.
#' @param data a data frame holding treatment variables and outcome.
#' @param outcome defines the outcome column
#' @param exposure defines the exposure variable (level)
#' @param bound maximum value for gamma and tau
#' @param symmetric identifies if wash-in effect is equal to wash-out effect
#' @return a vector with best gamma and tau.
#' @examples
#' load(simpatdat)
#' bound <- 5
#' estimate_gamma_tau(patient, "Uncertain_Low_Back_Pain", exposure = "treatment")
estimate_gamma_tau <- function(patients_data, outcome, exposure, bound=10, symmetric=TRUE) {

 res <- prep.onehot(patients_data, exposure)
 patients_data <- res$data
 exposure.columns <- res$names


  # defines a helper function to evaluate the gamma and tau
  evaluate_gamma_tau <- function(data, exposure.columns, outcome, gamma_1, tau_1, gamma_2, tau_2){
    for (exposure.column in exposure.columns){
      data[exposure.column] <- gen_treatment_effect(patients_data[exposure.column], gamma=gamma_1, tau=tau_1)
    }
    str_formula <- sprintf("%s ~ %s", outcome, paste(exposure.column, sep=" + ", collapse = " + "))
    lin_m <- lm(formula(str_formula), data = patients_data, na.action = na.omit)
    return(lin_m)
    }

  # define array for later evaluation
  if(symmetric){
    r2 <- array(NA,dim=rep(bound,2))
    #iterate over values
    for(gamma_1 in c(1:bound)){
        for(gamma_2 in c(1:bound)){
            r2[gamma_1, gamma_2] <- summary(evaluate_gamma_tau(patients_data, exposure.columns, outcome, gamma_1, gamma_1, gamma_2, gamma_2))$adj.r.squared
          }
    }
    # get the best model by the smallest r2 value
    best <- which(r2 == max(r2), arr.ind = TRUE)
    gamma_1 <- best[1]
    tau_1 <- best[1]
    gamma_2<- best[2]
    tau_2 <- best[2]
  }else{
    r2 <- array(NA,dim=rep(bound,4))
    #iterate over values
    for(gamma_1 in c(1:bound)){
      for(tau_1 in c(1:bound)){
        for(gamma_2 in c(1:bound)){
          for(tau_2 in c(1:bound)){
            r2[gamma_1, tau_1, gamma_2, tau_2] <- summary(evaluate_gamma_tau(patients_data, exposure.columns, outcome,gamma_1, tau_1, gamma_2, tau_2))$adj.r.squared
          }
        }
      }
    }
    # get the best model by the smallest r2 value
    best <- which(r2 == max(r2), arr.ind = TRUE)
    gamma_1 <- best[1]
    tau_1 <- best[2]
    gamma_2<- best[3]
    tau_2 <- best[4]
  }
  return(c(gamma_1, gamma_2, tau_1, tau_2))
}



#' #' Analyse Study
#' #' @description This function is a general function to estimate gamma and tau for a data frame.
#' #' @param study r data frame with patients with variables Treatment_1, Treatment_2 and outcome
#' #' @param symmetric identifies if wash-in effect is equal to wash-out effect
#' #' @return experiments with several gammas and taus.
#' #' @examples
#' #' load(simpatdat)
#' #' analyse_study(simpatdat, true)
#' #' @export
#' analyse_study <- function(study, symmetric, outcome, exposure) {
#'   # defines a helper function
#'   gen_treatment_effect <- function(variable, gamma, tau) {
#'     est_effect <- rep(NA, length(variable))
#'     for(v in c(1:length(variable))){
#'       if (v==1){
#'         x_j<-0
#'       }else{
#'         x_j<-est_effect[v-1]
#'       }
#'       est_effect[v] = (x_j + ((1 - x_j) / tau) * variable[v] - (x_j / gamma) * (1 - variable[v]))
#'     }
#'     return(est_effect)
#'   }
#'   columns <- c("Pat_ID",
#'                "Gamma_1",
#'                "Tau_1",
#'                "Gamma_2",
#'                "Tau_2",
#'                "est_t1",
#'                "est_t2",
#'                "est_t1_2.5",
#'                "est_t1_97.5",
#'                "est_t2_2.5",
#'                "est_t2_97.5",
#'                "Treatment_1",
#'                "Treatment_2",
#'                "Treatment_1_2.5",
#'                "Treatment_1_97.5",
#'                "Treatment_2_2.5",
#'                "Treatment_2_97.5")
#'
#'
#'   experiments <- data.frame(matrix(ncol=length(columns), nrow=0))
#'   names(experiments) <- columns
#'
#'   iterations <- 300
#'   pb <- txtProgressBar(max = iterations, style = 3)
#'   progress <- function(n) setTxtProgressBar(pb, n)
#'   opts <- list(progress = progress)
#'
#'   experiments <- foreach(patient_id=0:299, .combine=rbind,
#'                          .options.snow = opts) %dopar% {
#'                            # Select Patient based on that ID
#'                            patient <- study[study$patient_id==patient_id,]
#'                            best <- estimate_gamma_tau(patient, symmetric)
#'
#'                            gamma_1 <- best[1]
#'                            tau_1 <- best[2]
#'                            gamma_2<- best[3]
#'                            tau_2 <- best[4]
#'
#'                            patient$t1 <- gen_treatment_effect(patient$Treatment_1, gamma=gamma_1, tau=tau_1)
#'                            patient$t2 <- gen_treatment_effect(patient$Treatment_2, gamma=gamma_2, tau=tau_2)
#'
#'                            l_reg <- lm(outcome~Treatment_1+Treatment_2, data = patient, na.action = na.omit)
#'                            interval <- confint(l_reg)
#'
#'                            l_reg_est <- lm(outcome~t1+t2, data = patient, na.action = na.omit)
#'                            interval_est <- confint(l_reg_est)
#'
#'                            add_row <- data.frame(patient_id,
#'                                                  gamma_1,
#'                                                  tau_1,
#'                                                  gamma_2,
#'                                                  tau_2,
#'                                                  l_reg_est$coefficients["t1"],
#'                                                  l_reg_est$coefficients["t2"],
#'                                                  interval_est["t1","2.5 %"],
#'                                                  interval_est["t1","97.5 %"],
#'                                                  interval_est["t2","2.5 %"],
#'                                                  interval_est["t2","97.5 %"],
#'                                                  l_reg$coefficients["Treatment_1"],
#'                                                  l_reg$coefficients["Treatment_2"],
#'                                                  interval["Treatment_1","2.5 %"],
#'                                                  interval["Treatment_1","97.5 %"],
#'                                                  interval["Treatment_2","2.5 %"],
#'                                                  interval["Treatment_2","97.5 %"])
#'
#'                            names(add_row) <- columns
#'                            add_row
#'                          }
#'   close(pb)
#'   return(experiments)
#' }
