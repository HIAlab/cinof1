##############################################################
# Helper to generate data
# 4/4/2021
# created by T. Gärtner (thomas.gaertner@student.hpi.de)
##############################################################

prep.onehot <- function(data, var){
  data_frame <- data[var]
  dummy <- dummyVars(" ~ .", data=data_frame)
  newdata <- data.frame(predict(dummy, newdata = data_frame))
  return(list(data = cbind(data, newdata), names = colnames(newdata)))
}
