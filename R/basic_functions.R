##############################################################
# Basic functions to analyse N-of-1
# 4/4/2021
# created by T. Gärtner (thomas.gaertner@student.hpi.de)
##############################################################


#' Comparative Plot
#' @description This function creates a comparative within a tretment variable.
#' @param data defines a data frame.
#' @param exposure defines the exposure variable
#' @param outcome defines the outcome variable
#' @examples
#' comparative.plot(data)
#' @export
comparative.plot <- function(data, exposure, outcome){
  # Create Formula
  str_formula <- sprintf("%s ~ %s", outcome, exposure)
  # plot results
  boxplot(as.formula(str_formula), data=simpatdat)
}

#' Wilcox test
#' @description This function performs a wilcox test to identify, if there is an effect of the exposure on the outcome.
#' @param data dataframe with trail data
#' @param exposure identifies the column of exposure variable
#' @param outcome identifies the column of outccome variable
#' @example
#' load(simpatdat)
#' wilcox.nofone(simpatdat, "treatment", "Uncertain_Low_Back_Pain)
#' @references``
#' This function uses the method wilcox.test from stats package. Run '?wilcox.test' for more information.
#' @export
wilcox.nofone <- function(data, exposure, outcome){
  # Create Formula
  str_formula <- sprintf("%s ~ %s", outcome, exposure)
  # Return Result
  return(wilcox.test(as.formula(str_formula), data = data, exact = FALSE, correct = FALSE, conf.int = FALSE))
}



