#` compare derived maturity distribution with the one observed
checkMaturity <-  function(df1, df2, moltmod, sd, mean){
  require(dplyr)
  sizes <- unique(c(df1$size, df2$size))
  pvec <- rep(0,length(sizes))
  for(i in 1:length(sizes)){
    df1 <- filter(df1, size==sizes[i])
    df2 <- filter(df2, size==sizes[i])
    if((nrow(df1)>0)){
      if(nrow(df2)>0){
        pvec[i] <- sum(df2$number)/(sum(df1$number)+sum(df2$number))
      } else{
        pvec[i] <- 0
      }
    } else{
      pvec[i] <- 1
    }
  }
  sdPre <- (sd)
  meanPre <- (mean)
  CW <- seq(20,120, length.out=100)
  plot(pnorm(CW,mean=meanPre,sd=sdPre)~CW, type="l")
  points(pvec~sizes, pch=19)
  print(pvec)
}