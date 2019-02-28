latGrid <- seq(66,53.69231, length.out=41)
longGrid <- seq(-179,-154.3846,length.out=41)
require(dplyr)
require(readr)
latGrid <- latGrid[1:40]
longGrid <- longGrid[1:40]
proj_years <- 2017:2022
filter
settleMap <- lapply(1:length(proj_years),makeMatrices, latGrid, longGrid, proj_years)



mort <- matrix(1,nrow=40,ncol=40)
time<-matrix(NA,nrow=40,ncol=40)
surv <- bigsettled%>% select(id,horizPos1, horizPos2,number, time)
survInd <- which((surv$horizPos1>min(longGrid))&(surv$horizPos1<max(longGrid))&
                   (surv$horizPos2>min(latGrid))&(surv$horizPos2<max(latGrid)))
longs <- surv[survInd,"horizPos1"]
lats <- surv[survInd,"horizPos2"]
times <- surv[survInd,"time"]
gridcell_long <- gridcell_lat <- rep(0,length(longs))
for(i in 1:length(lats)){
  gridcell_long[i] <- last(which(longs[i]>longGrid))
  gridcell_lat[i] <- last(which(lats[i]<latGrid))
  mort[gridcell_lat[i],gridcell_long[i]]<-1-surv[survInd,"number"][i]
  time[gridcell_lat[i],gridcell_long[i]]<- substr(surv[survInd,"time"][i],1,4)
}
write.table(mort,
            file=paste("MortProportion2018.csv"),
            col.names=as.character(round(longGrid,2)), 
            row.names=as.character(round(latGrid,2)), sep=",")

write.table(time,
            file=paste("YearSettled2018.csv"),
            col.names=as.character(round(longGrid,2)), 
            row.names=as.character(round(latGrid,2)), sep=",")
}
