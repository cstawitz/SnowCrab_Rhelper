#' @description function to get centroid differences between the first and second data set for a subset of years and sizes
getCentroidDifference <- function(TheData1, TheData2, years, sizes){
  
  if(length(sizes)<=1){
    centroidLat <- centroidLong <- matrix(nrow=length(years)-1, ncol=2)} else{
      centroidLat <- centroidLong <- array(data=NA, dim=c(length(years)-1, length(sizes)-1,2))
    }
  for(i in 1:(length(years)-1)){
    early <- filter(TheData1, YEAR==years[i])
    late <- filter(TheData2, YEAR==years[i+1])
    if(length(sizes)>1){
      for(j in 1:(length(sizes)-1)){
        earlyCurr <- filter(early, SIZE==sizes[j])
        lateCurr <- filter(late, SIZE==sizes[j+1])
        if((nrow(earlyCurr)>0)&(nrow(lateCurr)>0)){
          centroidLat[i,j,1]<- sum(earlyCurr$numCPUE*earlyCurr$LATITUDE, na.rm=TRUE)/sum(earlyCurr$numCPUE, na.rm=TRUE)
          centroidLat[i,j,2]<- sum(lateCurr$numCPUE*lateCurr$LATITUDE, na.rm=TRUE)/sum(lateCurr$numCPUE, na.rm=TRUE)
          centroidLong[i,j,1] <- sum(earlyCurr$numCPUE*earlyCurr$LONGITUDE, na.rm=TRUE)/sum(earlyCurr$numCPUE, na.rm=TRUE)
          centroidLong[i,j,2] <- sum(lateCurr$numCPUE*lateCurr$LONGITUDE, na.rm=TRUE)/sum(lateCurr$numCPUE, na.rm=TRUE)
        }
      }
    } else{
      
      earlyCurr <- filter(early, SIZE==sizes)
      lateCurr <- filter(late, SIZE==(sizes+5))
      if((nrow(earlyCurr)>0)&(nrow(lateCurr)>0)){
        centroidLat[i,1]<- sum(earlyCurr$numCPUE*earlyCurr$LATITUDE, na.rm=TRUE)/sum(earlyCurr$numCPUE, na.rm=TRUE)
        centroidLat[i,2]<- sum(lateCurr$numCPUE*lateCurr$LATITUDE, na.rm=TRUE)/sum(lateCurr$numCPUE, na.rm=TRUE)
        centroidLong[i,1] <- sum(earlyCurr$numCPUE*earlyCurr$LONGITUDE, na.rm=TRUE)/sum(earlyCurr$numCPUE, na.rm=TRUE)
        centroidLong[i,2] <- sum(lateCurr$numCPUE*lateCurr$LONGITUDE, na.rm=TRUE)/sum(lateCurr$numCPUE, na.rm=TRUE)
      }
    }
  }
  return(list("lat"=centroidLat,"long"=centroidLong))
}