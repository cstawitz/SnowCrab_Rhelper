setwd("C:/Users/Christine.Stawitz/Documents/SurveyData")
require(dplyr)
CPUEYear<- read.csv("CPUE.csv")
years <- unique(CPUEYear$YEAR)
for(i in 1:length(years)){
  eachYr<- filter(CPUEYear, YEAR==years[i]) %>%
    filter(numCPUE>0)
  newshellmale <- eachYr %>% filter(SEX=="MALE") %>%
    filter(SHELL_CONDITION=="NEW_SHELL")
newshellmale %>% group_by(LATITUDE, LONGITUDE) %>%
  summarise(tot=sum(numCPUE), meansize = average(WIDTH))
}

require(dplyr)
mat<- filter(CPUEYear, MATURITY=="MATURE")
immat <- filter(CPUEYear, MATURITY=="IMMATURE")

yrs <- seq(2006,2012)

compareCentroids <- function(mat, immat, survey, yrs){
SurveyFem <- modelFem <- SurveyImm <- modelImm <- vector("list")
meanLong <- meanLat <- meanModLong <- meanModLat <- meanLongImm <- meanLatImm <- meanModImmLong <- meanModImmLat <- rep(0,6)
for(i in 1:length(yrs)){
  SurveyFem[[i]] <- filter(mat, YEAR==yrs[i])
  SurveyImm[[i]] <- filter(immat, YEAR==yrs[i])
  meanLong[i] <- sum(SurveyFem[[i]]$numCPUE*SurveyFem[[i]]$LONGITUDE)/sum(SurveyFem[[i]]$numCPUE)
  meanLat[i] <- sum(SurveyFem[[i]]$numCPUE*SurveyFem[[i]]$LATITUDE)/sum(SurveyFem[[i]]$numCPUE)
  meanLongImm[i] <- sum(SurveyImm[[i]]$numCPUE*SurveyImm[[i]]$LONGITUDE)/sum(SurveyImm[[i]]$numCPUE)
  meanLatImm[i] <- sum(SurveyImm[[i]]$numCPUE*SurveyImm[[i]]$LATITUDE)/sum(SurveyImm[[i]]$numCPUE)
  
  modelFem[[i]] <- filter(junJulyFem, substr(time, 1,4)==yrs[i])
  meanModLong[i] <- sum( modelFem[[i]]$number*modelFem[[i]]$horizPos1)/sum(modelFem[[i]]$number)
  meanModLat[i] <- sum( modelFem[[i]]$number*modelFem[[i]]$horizPos2)/sum(modelFem[[i]]$number)
  
  modelImm[[i]] <- filter(junJulyImm, substr(time, 1,4)==yrs[i])
  meanModImmLong[i] <- sum( modelImm[[i]]$number*modelImm[[i]]$horizPos1)/sum(modelImm[[i]]$number)
  meanModImmLat[i] <- sum( modelImm[[i]]$number*modelImm[[i]]$horizPos2)/sum(modelImm[[i]]$number)
}
allLong <- c(meanModImmLong,meanModLong,meanLongImm,meanLong)
allLat <- c(meanModImmLat, meanModLat, meanLatImm,meanLat)
minLat <- min(allLat, na.rm=T)
maxLat <- max(allLat, na.rm=T)
minLong <- min(allLong, na.rm=T)
maxLong <- max(allLong, na.rm=T)
plot(x=meanLat, y=meanLong, xlim=c(minLat,maxLat),
     ylim=c(minLong,maxLong), pch=as.character(c(6,7,8,9,0,1,2)))
points(x=meanModLat, y=meanModLong, col="red", pch=as.character(c(6,7,8,9,0,1,2)))
points(x=meanModImmLat, y=meanModImmLong, col="blue", pch=as.character(c(6,7,8,9,0,1,2)))
points(x=meanLatImm, y=meanLongImm, col="green", pch=as.character(c(6,7,8,9,0,1,2)))
}

pos <- filter(SurveyImm[[1]], numCPUE>0)
sizelatlong<- pos %>% group_by(LONGITUDE,LATITUDE) %>% summarise(mean(SIZE), sum(numCPUE))
write.csv(sizelatlong, "ImmatFemaleSize2006.csv")
unique(pos$LONGITUDE)
