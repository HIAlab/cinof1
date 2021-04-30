##############################################################
# G-Estimation
# 4/4/2021
# created by T. Gärtner (thomas.gaertner@student.hpi.de)
##############################################################

#' Calc beta
#' @description This function fits a geeglm to the data for a given psi. it returns a list of values
#' @param data data frame with values
#' @param psi number, given psi
#' @param outcome string, defines the outcome variable
#' @param exposure string, defines the exposure variable
#' @param confounder string vector, defines the confounder variables
#' @param id string, identify the unique identifier
calc_beta <- function(data, psi, outcome, exposure, confounder, id){
  data$Hpsi <- as.numeric(data[,outcome])-psi*as.numeric(data[,exposure])
  summary(data)
  str_formula <- sprintf("%s ~ %s", exposure, paste(confounder, "Hpsi", sep=" + ", collapse = " + "))
  g.est <- geepack::geeglm(as.formula(str_formula), family=gaussian, data=data, id = data[[id]], corstr="ar1")
  beta <- summary(g.est)$coefficients["Hpsi","Estimate"]
  se <- summary(g.est)$coefficients["Hpsi","Std.err"]
  return(list(beta = beta, se = se))
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
  steps=10,
  verbose=TRUE){

  # handle missing values
  if(missing(max_number_it)) {
    max_number_it <- NULL
  }
  epsilon <- .Machine$double.eps

  # prepare data: y must be numeric
  data[,exposure] <- as.numeric(data[,exposure])


  # Default values
  number_it <- 0

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
    res <- calc_beta(data = data, psi = upper_bound_psi, outcome, exposure, confounder, id=id)
    upper_bound_beta <- res$beta
    upper_se <- res$se

    res <- calc_beta(data = data, psi = lower_bound_psi, outcome, exposure, confounder, id=id)
    lower_bound_beta <-  res$beta
    lower_se <- res$se

    # return rec search
    df <- find_optimum_rec(data,
                     outcome,
                     exposure,
                     confounder,
                     id,
                     upper_bound_psi,
                     upper_bound_beta,
                     upper_se,
                     lower_bound_psi,
                     lower_bound_beta,
                     lower_se,
                     number_it+1,
                     max_number_it,
                     epsilon,
                     method=method,
                     verbose=verbose)

    df[nrow(df) + 1,] <- c(upper_bound_psi, upper_bound_beta, upper_se, lower_bound_psi, lower_bound_beta, lower_se, number_it)
    return(df)
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
  upper_se=0,
  lower_bound_psi,
  lower_bound_beta,
  lower_se=0,
  number_it,
  max_number_it,
  epsilon,
  method,
  verbose=TRUE){


  process_return_values <- function(){
    df <- data.frame(upper_psi = double(),
               upper_beta = double(),
               upper_se = double(),
               lower_psi = double(),
               lower_beta = double(),
               lower_se = double(),
               iteration = double())
    df[1,] <- c(upper_bound_psi, upper_bound_beta,  upper_se, lower_bound_psi, lower_bound_beta,lower_se, number_it)
    return(df)
  }

  if(abs(upper_bound_beta) < epsilon){
    print(paste("Converged! Optimal Psi: ", upper_bound_psi))
    return(process_return_values())
  }else if(abs(lower_bound_beta) < epsilon){
    print(paste("Converged! Optimal Psi: ", lower_bound_psi))
    return(process_return_values())
  }

  if(!is.null(max_number_it)){
    if(number_it == max_number_it){
    print("Max Iterations:", number_it)
    return(process_return_values())
    }
  }

  if (verbose){
    vec <- c(paste("Iteration:", number_it),
             paste("Upper Psi:", upper_bound_psi),
             paste("Upper Beta:", upper_bound_beta),
             paste("Lower Psi:", lower_bound_psi),
             paste("Lower Beta:", lower_bound_beta),
             "")
    writeLines(paste(vec, sep="\n", collapse = "\n"))
  }


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

  res <- calc_beta(data = data, psi = new_psi, outcome, exposure, confounder, id)

  new_beta <- res$beta
  new_se <- res$se

  if (((upper_bound_beta >= 0) & (new_beta >=0) ) | ((upper_bound_beta <= 0) & (new_beta <=0))){
    upper_bound_psi <- new_psi
    upper_bound_beta <- new_beta
    upper_se <- new_se
  }else if(((lower_bound_beta >= 0) & (new_beta >=0)) | ((lower_bound_beta <= 0) & (new_beta <=0))){
    lower_bound_psi <- new_psi
    lower_bound_beta <- new_beta
    lower_bound_se <- new_se
  }else{
    print(paste("Error at iteration ", number_it))
    return(process_return_values())
  }

  df <- find_optimum_rec(data = data,
                         outcome = outcome,
                         exposure = exposure,
                         confounder = confounder,
                         id = id,
                         upper_bound_psi = upper_bound_psi,
                         upper_bound_beta = upper_bound_beta,
                         upper_se = upper_se,
                         lower_bound_psi = lower_bound_psi,
                         lower_bound_beta = lower_bound_beta,
                         lower_se = lower_se,
                         number_it = number_it + 1,
                         max_number_it = max_number_it,
                         epsilon = epsilon,
                         method=method,
                         verbose = verbose)
  df[nrow(df) + 1,] <-  c(upper_bound_psi, upper_bound_beta, upper_se, lower_bound_psi, lower_bound_beta, lower_se, number_it)
  return(df)

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

  library(doParallel)

  cores=detectCores()
  cl <- makeForkCluster(cores[1]-1) #not to overload your computer
  registerDoParallel(cl)


  df <- foreach::foreach(i = c(1:length(psis)), .combine ="rbind") %dopar% {
    psi <- psis[[i]]
    res <- calc_beta(data = data, psi = psi, outcome, exposure, confounder, id)
    beta <- res$beta
    se <- res$se
    c(psi, beta, se)
  }

  # Stop Cluster
  stopCluster(cl)

  colnames(df) <- c("PSI","Beta", "Std.Err")
  return(data.frame(df))
}
