#Get netcdf files for batch runs
install.packages("raster")
install.packages("ncdf4")
library(ncdf4)
r <- nc_open("F:/ROMS.BeringSeaForecasts/Bering_grid_10k.nc")
str(r)
#By Year
latGrid <- seq(66,53.69231, length.out=41)
longGrid <- getLongGrid(latGrid)
longGrid <- seq(-179,-154.3846,length.out=41)
midLat <- ((latGrid[2:41]-latGrid[1:40])/2)+latGrid[1:40]
midLong <- ((longGrid[2:41]-longGrid[1:40])/2)+longGrid[1:40]
crabLoc <- expand.grid(midLong, midLat)

lat <- seq(66,54,length.out=40)
lon <- seq(-179,-155, length.out=40)


library(gdistance)
library(raster)
library(maps)
library(maptools)

data(wrld_simpl)


pts <- SpatialPoints(crabLoc, proj4string=CRS(proj4string(wrld_simpl)))

## Find which points fall over land
land <- matrix(!is.na(over(pts, wrld_simpl)$FIPS), nrow=1600, ncol=2)
ind <- which(land[,1])
crabLoc[-ind,]
pts2 <- SpatialPoints(crabLoc[-ind,], proj4string=CRS(proj4string(wrld_simpl)))
## Check that it worked
plot(wrld_simpl,xlim = c(min(longGrid), max(longGrid)), ylim = c(min(latGrid),max(latGrid)))
points(pts, col=1+land, pch=16)
points(pts2, col=3, pch=2)
write.csv(crabLoc[-ind,],"CrabLatLong.csv")

# ncrab <- rep(0,nrow(longGrid))
# for(i in 1:nrow(longGrid)){
#   ncrab[i] <- sum(longGrid[i,]<(-155))
# }
# inds <- seq(1,39,length.out=20)
# whichLongs <- longGrid[inds,]
# whichLats <- latGrid[inds]
# 
# latlong <- matrix(nrow=1600,ncol=2)
# for(i in 1:nrow(whichLongs)){
#   latlong[(1+(i-1)*80):(80+(i-1)*80),] <- cbind(whichLats, whichLongs[i,])
# }
# 
# latLong<-latlong[latlong[,2]<(-155),]
write.csv(latLong, "LatLongGrid.csv")


start <- read.csv("initAttsAllLoc.csv")
max(start$Horiz..position.1)
max(start$Horiz..position.2)



read.csv("MortProportion2018.csv", header=T)

setwd("C:\\Users\\Christine.Stawitz\\Documents\\SnowCrab_InputFiles\\")
save(settletime, file="Settletimes.RData")
load("Settletimes.RData")
settletime[[5]]
k <- 2020
lab <- rep("",33)
for(i in 1:33){
  lab[i] <- paste(k, "- ", k+1)
  k <- k+1
}

par(oma=c(1,3,0,3))
plot(NA, ylim=c(1,33), xlim=c(1,24), axes=F, xlab="Month", ylab="")
axis(side=1, labels=c("Jan","Feb", "May", "July", "Sept", "|" ,"Jan", "Feb", "May", "July","Sept"), at=c(1,3,5,7,9,12.5,13,15,17,19,21))
#axis(side=4, labels=rev(lab), at=seq(1,33), las=1, cex.axis=.8)
axis(side=2, labels=2050:2018, at=1:33, las=1)
for(i in c(1,3:33)){
  yr <- substr(settletime[[i]],1,4)
  mn <- substr(settletime[[i]],6,7)
  if(length(unique(yr))==2){
    val <- ifelse(substr(settletime[[i]],1,4)==min(yr), 0,1)
  } else{
    val <-1
  }
  points(x=jitter(val*12+as.numeric(mn),1), y=rep(length(proj)-i+1, length(settletime[[i]])))
  text(unique(yr)[1], x=23, y=length(proj)-i+1+.3)
  text(unique(yr)[2], x=24, y=length(proj)-i+1+.3)
}
mtext(side=2, "Birth Year", line=4)
text("Rec Year(s)", x=23.5, y=34)
mtext("Crab instar 8 timing", side=3)

plot(as.numeric(substr(time,1,4))~horizPos1, data=bigsettled)

longGrid <- longGrid[1:40]
latGrid <- latGrid[1:40]
ind<-matrix(NA, nrow=2, ncol=length(tofind))
probBylat <- vector("list")
for(j in 1:6){
  probBylat[[j]] <- array("",dim=c(length(longGrid),length(latGrid),nrow(tofind)))
}
maxLong <- longGrid[40]+(longGrid[40]-longGrid[39])

k<-0
mort <- rep(0,6)
tofind <- filter(zooea1,ageInStage==0)
bigsettled <- filter(immat_male, ageInStage==0)
for(i in 1:nrow(tofind)){
  
  yrind <- switch(substr(tofind[i,"time"], 1,4), "2017"=1, "2018"=2, "2019"=3)
  if(length(which(bigsettled$id==tofind[i,'id']))==1){
    ind[i,1] <- which(bigsettled$id==tofind[i,'id'])
    long2 <- filter(bigsettled,id==tofind[i,'id'])$horizPos1
    lat2 <- filter(bigsettled,id==tofind[i,'id'])$horizPos2
    long1 <- filter(zooea1_Start,id==tofind[i,'id'])$horizPos1
    lat1 <- filter(zooea1_Start,id==tofind[i,'id'])$horizPos2
  } else if(length(which(bigsettled$id==tofind[i,'id']))>1){
    long <- filter(bigsettled,id==tofind[i,'id'])$horizPos1
    lat <- filter(bigsettled,id==tofind[i,'id'])$horizPos2
    long1 <- filter(zooea1_Start, id==tofind[i,'id'])$horizPos1
    lat1 <- filter(zooea1_Start, id==tofind[i,'id'])$horizPos2
    if((long>min(longGrid))&&(long<max(longGrid))&&(lat>min(latGrid))&&(lat<max(latGrid))){
      y <- last(which(long1>longGrid))+1
      x <- first(which(lat1>latGrid))
      y2 <- last(which(long>longGrid))+1
      x2 <- first(which(lat>latGrid))
      
      probBylat[[yrind]][x,y,i] <- paste(x2,y2, sep=",")
    } else{
      mort[yrind] = mort[yrind] + filter(tofind,id==tofind[i,"id"])$number
      ind[i]<-NA
    }
  } else{
    mort[yrind] = mort[yrind] + filter(tofind,id==tofind[i,"id"])$number
    ind[i]<-NA
  }
}

sumProb <- probBylat[,,1]
for(i in 1:(dim(probBylat)[3])){
  tmp <- probBylat[,,i][probBylat[,,i]!=""]
  if(length(tmp)>0){
    sumProb[probBylat[,,i]!=""]<-paste(sumProb[probBylat[,,i]!=""], probBylat[,,i][probBylat[,,i]!=""])
  }
}
settlementDates<- table(settled$time)
timeofsettlement <- substr(as.character(settled$time),6,7)
plyr::count(timeofsettlement)
str(settled)
max(settled$age)
plot(settlementDates, axes=F, ylab="Frequency")
axis(1, labels=substr(names(settlementDates)[seq(1,2117, length.out=15)],1,7), at=seq(1,2117, length.out=15),
     las=2)
axis(2)

zooea1SettlementMap <- array(dim=c(40,40,1600))
for(i in 1:(dim(zooea1SettlementMap)[3])){
  if(sumProb[i]!=""){
    transitionList<- unlist(strsplit(sumProb[i], " "))
    transitionList <- transitionList[transitionList!=""]
    map <- table(transitionList)
    coords<- as.numeric(unlist(strsplit(names(map),",")))
    if(length(coords)==2){
      zooea1SettlementMap[coords[1], coords[2],i] <- 1
    } else if(length(coords)>2){
      coords <- matrix(coords, ncol=2, byrow=T) 
      
      for(j in 1:nrow(coords)){
        zooea1SettlementMap[coords[j,1], coords[j,2],i] <- map[j]/sum(map)
      }
    }
  }
}
save(zooea1SettlementMap,file="Zooea1toSettledArray.RData")
load("Zooea1toSettledArray.RData")
zooea1SettlementMap[is.na(zooea1SettlementMap)]<-rep(0, length(is.na(zooea1SettlementMap)))
transitionList<- strsplit(sumProb[sumProb!=""], " ")


for(i in 1:length(transitionList)){
  places <- table(transitionList[[i]])
}




transitions <- as.data.frame(sumProb, row.names=as.character(latGrid))
names(transitions) <- as.character(longGrid)
write.csv(transitions,"Transition.csv")

locTransition<- apply(X=probBylat, MARGIN=c(1,2), FUN=paste, na.rm=T)

probBylat <- as.numeric(probBylat)

oceanTimes <- getOceanTimes(topdir = "F:/ROMS.BeringSeaForecasts/miroc-rcp85")
names(oceanTimes)
oceanTimes$date
batchmode <- read.csv("BatchModeRuns.csv")
starttimes <- batchmode$startTime
dir <- rep("", length(starttimes))
for(i in 1:length(starttimes)){
  rowind <- which(oceanTimes$ocean_time> starttimes[i])[1]-1
  print(oceanTimes$date[rowind])
  dir[i] <- oceanTimes$name[rowind]
}
rowind <- which((substr(oceanTimes$date,6,9)=="03-0"))
dates <- oceanTimes[rowind,]
batchmode <- data.frame(rep(NA, 119))
batchmode$subDirectory<-substr(dates$date,1,4)
batchmode$startTime <- dates$ocean_time
batchmode$file_ROMSDataset <- file.path("F:/ROMS.BeringSeaForecasts/miroc-rcp85",dates$name)
write.csv(batchmode, "BatchModeRuns.csv")
