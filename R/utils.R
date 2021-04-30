##############################################################
# Helper to generate data
# 4/4/2021
# created by T. Gärtner (thomas.gaertner@student.hpi.de)
##############################################################

prep.onehot <- function(data, var){
  data[,var] <- as.factor(data[,var])
  levels <- levels(data[,var])
  new_variables <- c()

  j <- 1

  for(i in c(1:length(levels))){
    if(levels[i]!="" && !is.na(levels[i])){
      data[,paste(var, levels[i], sep=".")] <- ifelse(is.na(data[,var]), NA, ifelse(data[,var]==levels[i], 1,0))
      new_variables[j] <- paste(var, levels[i], sep=".")
      j <- j + 1
    }
  }
  return(list(data = data, names = new_variables))
}


prep.normalize <- function(data, scale="sd"){
  if(scale=="sd"){
    return((data-mean(data))/sd(data))
  }else if(scale=="uni"){
    return((data-min(data))/(max(data)-min(data)))
  }else{
    return(data())
  }
}
