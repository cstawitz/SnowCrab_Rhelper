GetDist <- function(TheData){
  distance <- angle <- TheData$lat[,,1]
  for(i in 1:(nrow(TheData$lat[,,1]))){
    for(j in 1:(ncol(TheData$lat[,,1])-1)){
      x <- rbind(c(TheData$lat[i,j,1], TheData$long[i,j,1]),
                 c(TheData$lat[i,j,2], TheData$long[i,j,2]))
      distance[i,j]  <- dist(x)
      dot.prod <- x[1,]%*%x[2,] 
      if(!any(is.na(x))){
        theta <- atan((x[2,1]-x[1,1])/(x[2,2]-x[1,2]))
        angle[i,j]  <- as.numeric(theta)
      } else{
        angle[i,j] <- 0
      }
    }
  }
  return(list("dist"=distance, "angle"=angle))}