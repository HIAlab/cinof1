##############################################################
# G-Estimation
# 4/4/2021
# created by T. Gärtner (thomas.gaertner@student.hpi.de)
##############################################################

#' Nofgest
#' @description This function use standard IV estimator via additive marginal structural models to analyse a treatment effect
#' @param data
#' @param outcome
#' @param exposure
#' @param confounder
#' @param id
#' @param number_it
#' @param upper_bound_psi
#' @param lower_bound_psi
#' @param max_number_it
#' @example
#' load(simpatdat)
#' nofgest()
#' @export
nofgest <- function(
  data,
  outcome,
  exposure,
  confounder,
  id,
  upper_bound_psi,
  lower_bound_psi,
  max_number_it){

  if(missing(epsilon)) {
    epsilon <- .Machine$double.eps
  }

  if(missing(max_number_it)) {
    max_number_it <- NULL
    epsilon <- 2.220446e-16
  }

  if(missing(number_it)) {
    number_it <- 1
  }

  if(missing(upper_bound_psi)){
    upper_bound_psi <- 10
  }

  if(missing(lower_bound_psi)){
    lower_bound_psi <- 10
  }

  if(missing(upper_bound_beta)) {
    upper_bound_beta <-  calc_beta(data = data, psi = upper_bound_psi, outcome, exposure, confounder, id="patient_id")[1]
  }

  if(missing(lower_bound_beta)) {
    lower_bound_beta <-  calc_beta(data = data, psi = lower_bound_psi, outcome, exposure, confounder, id="patient_id")[1]
  }

  return(find_optimum_rec(data,
                          outcome,
                          exposure,
                          confounder,
                          id,
                          upper_bound_psi,
                          upper_bound_beta,
                          lower_bound_psi,
                          lower_bound_beta,
                          number_it,
                          max_number_it,
                           epsilon))
}



#' Find Optimum Recursively
#' @description This function finds the best psi value recursively
#' @param data
#' @param outcome
#' @param exposure
#' @param confounder
#' @param id
#' @param number_it
#' @param upper_bound_psi
#' @param lower_bound_psi
#' @param max_number_it
find_optimum_rec <- function(
  data,
  outcome,
  exposure,
  confounder,
  id,
  upper_bound_psi,
  upper_bound_beta,
  lower_bound_psi,
  lower_bound_beta,
  number_it,
  max_number_it,
  epsilon){


  calc_beta <- function(data, psi, outcome, exposure, confounder, id){
    data$Hpsi <- as.numeric(data[,outcome])-psi*as.numeric(data[,exposure])

    str_formula <- sprintf("%s ~ %s", exposure, paste(confounder, "Hpsi", sep=" + ", collapse = " + "))

    g.est <- geeglm(as.formula(str_formula), family=gaussian, data=data, id = data[[id]], corstr="ar1")

    beta <- summary(g.est)$coefficients["Hpsi","Estimate"]

    return(beta)
  }

  process_return_values <- function(upper_bound_psi, upper_bound_beta, lower_bound_psi,lower_bound_beta){
    if(abs(upper_bound_beta)>=abs(lower_bound_beta)){
      return(data.frame(lower_bound_psi))
    }else{
      return(data.frame(lower_bound_psi))
    }
  }

  number_it <- number_it+1
  if(abs(upper_bound_psi-lower_bound_psi) <=epsilon){
    print("Converged!")
    return(process_return_values(upper_bound_beta, lower_bound_beta, upper_bound_psi, lower_bound_psi))
  }

  if(number_it == max_number_it){
    print("Max Iterations")
    return(process_return_values(upper_bound_beta, lower_bound_beta, upper_bound_psi, lower_bound_psi))
  }
  if(upper_bound_beta<= lower_bound_beta | upper_bound_beta<0 | lower_bound_beta>0){
    print("Error")
    return(process_return_values(upper_bound_beta, lower_bound_beta, upper_bound_beta, lower_bound_beta))
  }
  psi <- mean(c(upper_bound_psi, lower_bound_psi))
  res <- calc_beta(data = data, psi = psi, outcome, exposure, confounder, id)
  new_beta <- res[1]
  if ((new_beta >=0) & (new_beta <upper_bound)){
    upper_bound_psi <- psi
    upper_bound_beta <- new_beta
    return(find_optimum_rec(data, outcome, exposure, confounder, id, upper_bound_psi, upper_bound_beta, lower_bound_psi,lower_bound_beta, number_it, max_number_it))
  }else if((new_beta <=0) & (new_beta >lower_bound)){
    lower_bound_psi <- psi
    lower_bound_beta <- new_beta
    return(find_optimum_rec(data, outcome, exposure, confounder, id, upper_bound_psi, upper_bound_beta, lower_bound_psi,lower_bound_beta, number_it, max_number_it))
  }else{
    return(process_return_values(upper_bound_beta, lower_bound_beta, upper_bound_beta, lower_bound_beta))
  }

}
