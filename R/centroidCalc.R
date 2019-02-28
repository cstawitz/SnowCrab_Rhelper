#########################################
# Calculate centroid
centroidCalc <- function(dataframe){
  yrs <- unique(dataframe$year)
  nyr <- length(yrs)
  require(dplyr)
  centMat <- matrix(rep(0, nyr*2), nrow=nyr)
  for(i in 1:nyr){
    annual <- filter(dataframe, year==yrs[i])
    
    centMat[i,1] <- sum(dataframe$number*dataframe$posLongs, na.rm=TRUE)/sum(dataframe$number, na.rm=TRUE)
    centMat[i,2] <- sum(dataframe$number*dataframe$horizPos2, na.rm=TRUE)/sum(dataframe$number, na.rm=TRUE)
  }
  return(centMat)
}