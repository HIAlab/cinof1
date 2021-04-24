##############################################################
# Helper to generate data
# 4/4/2021
# created by T. Gärtner (thomas.gaertner@student.hpi.de)
##############################################################

prep.onehot <- function(data, var){
  data[,var] <- as.factor(data[,var])
  levels <- levels(data[,var])
  new_variables <- c()

  for(i in c(1:length(levels))){
    if(levels[i]!="" && !is.na(levels[i])){
      data[,paste(var, levels[i], sep=".")] <- ifelse(data[,var]==levels[i], 1,0)
      new_variables[i] <- paste(var, levels[i], sep=".")
    }
  }
  return(list(data = data, names = new_variables))
}


prep.normalize <- function(data){
  return((data-mean(data))/sd(data))
}
