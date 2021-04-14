##############################################################
# G-Estimation
# 4/4/2021
# created by T. Gärtner (thomas.gaertner@student.hpi.de)
##############################################################

#' Calc beta
#' @export
calc_beta <- function(data, psi, outcome, exposure, confounder, id){
  data$Hpsi <- as.numeric(data[,outcome])-psi*as.numeric(data[,exposure])
  summary(data)
  str_formula <- sprintf("%s ~ %s", exposure, paste(confounder, "Hpsi", sep=" + ", collapse = " + "))
  g.est <- geepack::geeglm(as.formula(str_formula), family=gaussian, data=data, id = data[[id]], corstr="ar1")
  beta <- summary(g.est)$coefficients["Hpsi","Estimate"]
  return(beta)
}




#' Nofgest
#' @description This function use standard IV estimator via additive marginal structural models to analyse a treatment effect
#' @param data data frame
#' @param outcome outcome variable name
#' @param exposure exposure variable name
#' @param confounder list of confounders
#' @param id patient id
#' @param max_number_it max umber of iterations
#' @examples
#' load(simpatdat)
#' nofgest()
#' @export
nofgest <- function(
  data,
  outcome,
  exposure,
  confounder,
  id,
  method="iterate",
  upper_bound_psi=2,
  lower_bound_psi=-2,
  max_number_it=NULL,
  steps=10){

  # handle missing values
  if(missing(max_number_it)) {
    max_number_it <- NULL
  }
  epsilon <- .Machine$double.eps

  # prepare data: y must be numeric
  data[,exposure] <- as.numeric(data[,exposure])


  # Default values
  number_it <- 1

  if(method=="iterate"){
    return(find_optimum_iterate(data,
                                outcome,
                                exposure,
                                confounder,
                                id,
                                upper_bound_psi,
                                lower_bound_psi,
                                steps))
  }else{
    upper_bound_beta <-  calc_beta(data = data, psi = upper_bound_psi, outcome, exposure, confounder, id=id)
    lower_bound_beta <-  calc_beta(data = data, psi = lower_bound_psi, outcome, exposure, confounder, id=id)

    # return rec search
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
                            epsilon,
                            method=method))
  }
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
  epsilon,
  method){


  process_return_values <- function(){
    print(paste("Number of iterations:", number_it))
    return(list(upper_psi = upper_bound_psi, upper_beta =upper_bound_beta, lower_psi = lower_bound_psi, lower_beta = lower_bound_beta, n_it=number_it))
  }


  if(abs(upper_bound_beta) <epsilon){
    print(paste("Converged! Optimal Psi: ", upper_bound_psi))
    return(process_return_values())
  }else if(abs(lower_bound_beta) < epsilon){
    print(paste("Converged! Optimal Psi: ", lower_bound_psi))
    return(process_return_values())
  }

  if(!is.null(max_number_it)){
    if(number_it == max_number_it){
    print("Max Iterations")
    return(process_return_values())
    }
  }

  number_it <- number_it+1

  if(method=="rec_mean"){
    new_psi <- mean(c(upper_bound_psi, lower_bound_psi))
  }else{
    beta_1 <- lower_bound_beta
    beta_2 <- upper_bound_beta
    psi_1 <- lower_bound_psi
    psi_2 <- upper_bound_psi
    m <- (beta_1 - beta_2)/(psi_1 - psi_2)
    a <- beta_1 -m*psi_1
    new_psi <- -a/m
  }

  new_beta <- calc_beta(data = data, psi = new_psi, outcome, exposure, confounder, id)

  if (((upper_bound_beta >= 0) & (new_beta >=0) ) | ((upper_bound_beta <= 0) & (new_beta <=0))){
    upper_bound_psi <- new_psi
    upper_bound_beta <- new_beta
    return(find_optimum_rec(data, outcome, exposure, confounder, id, upper_bound_psi, upper_bound_beta, lower_bound_psi,lower_bound_beta, number_it, max_number_it, epsilon, method=method))
  }else if(((lower_bound_beta >= 0) & (new_beta >=0)) | ((lower_bound_beta <= 0) & (new_beta <=0))){
    lower_bound_psi <- new_psi
    lower_bound_beta <- new_beta
    return(find_optimum_rec(data, outcome, exposure, confounder, id, upper_bound_psi, upper_bound_beta, lower_bound_psi,lower_bound_beta, number_it, max_number_it, epsilon, method=method))
  }else{
    print(paste("Error at iteration ", number_it))
    return(process_return_values())
  }

}



#' Find Optimum Brute Force
#' @description This function finds the best psi value brute forced
#' @param data
#' @param outcome
#' @param exposure
#' @param confounder
#' @param id
#' @param number_it
#' @param upper_bound_psi
#' @param lower_bound_psi
#' @param max_number_it
find_optimum_iterate <- function(data,
                                 outcome,
                                 exposure,
                                 confounder,
                                 id,
                                 upper_bound_psi,
                                 lower_bound_psi,
                                 steps){

  psis <- seq(upper_bound_psi, lower_bound_psi, length = steps)
  betas <- c()

  for(i in 1:length(psis)){
    psi <- psis[[i]]
    beta <- calc_beta(data = data, psi = psi, outcome, exposure, confounder, id)
    betas[i]<- beta
  }
  df <- data.frame(psis, betas)
  colnames(df) <- c("PSI","Beta")
  return(df)
}
