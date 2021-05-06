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
#' @param method one of `boxplot` (default) or  `denseplot`
#' @param adjustment adjust parameter for ggplot2::geom_density
#' @param title title string for plot
#' @examples
#' # Define outcome and exposure column
#' outcome <- "Uncertain_Low_Back_Pain"
#' exposure <- "treatment"
#' # Plot outcome among different exposures
#' comparative.plot(simpatdat, exposure = exposure, outcome = outcome)
#' @export
comparative.plot <- function(data, exposure, outcome, method="boxplot", adjustment=3, title = "Comparative Plot"){
  # plot results

  # Create mean
  labels <- levels(data[,exposure])

  mu <- data.frame(row.names = c(1:length(labels)))
  mu[,exposure] <- labels
  mu[,"mean"] <- rep(NA, length(labels))


  for(i in c(1:length(labels))){
    mu[i,"mean"] <- (mean(na.omit(data[data[,exposure]==labels[i],outcome])))
  }


  if(method=="denseplot"){
    p <- ggplot2::ggplot(data=data, mapping=ggplot2::aes_string(x=outcome, color = exposure , fill=exposure)) +
      ggplot2::geom_density( color="#e9ecef", alpha=0.6, position = 'identity', adjust = adjustment ) +
      ggplot2::labs(fill="") +
      ggplot2::ggtitle(label=title, subtitle = paste("Mean difference: ", round(abs(mu[1,"mean"]-mu[2,"mean"]), digits = 3))) +
      ggplot2::geom_vline(data=mu, ggplot2::aes(xintercept=mean),
                 linetype="dashed") +
      ggplot2::scale_color_brewer(palette="Dark2") +
      ggplot2::theme(legend.position="top")

    return(p)
  }else{
    p<-ggplot2::ggplot(data=data, mapping=ggplot2::aes_string(x=exposure, y=outcome, color=outcome)) +
      ggplot2::geom_boxplot() +
      ggplot2::ggtitle(label=title, subtitle = paste("Mean difference: ", round(abs(mu[1,"mean"]-mu[2,"mean"]), digits = 3)))

    return(p)
  }
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



