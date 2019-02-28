#'@description Function to create transition matrices from a number of DisMELS batch runs and write them to output files
#'@param i - index of each batch run, where i is within the dimensions of the proj_years input vector
#'@param latGrid - A vector of latitudes specifying the minimum latitude for each grid cell
#'@param longGrid - A vector of longitudes specifying the minimum longitude for each grid cell
#'@param proj_years - vector of years output is available for; these are assumed to correspond 
#'@param dir.batchruns - the parent directory all of the yearly batch run outputs are stored in
#'to the folder names of the batch runs, where each folder proj_years[i] contains a ConnRes.csv file
#'@return a list; the first element is the mortality rate and the second is the year of settlement

makeMatrices <- function(i, latGrid, longGrid,  proj_years, dir.batchruns){
  
  require(readr)
  require(dplyr)
  require(data.table)
  
  #Set to directory that output files are in; assumes each batch run output is its own folder with the name of the year
  setwd(paste(dir.batchruns,proj_years[i],sep=""))
  
  #Read in data
  allresults <- fread("ConnRes.csv", fill=T)
  
  #Filter results to look at immature males at time of settlement
  settled <- filter(allresults, typeName =="ImmatureMale") %>% filter(alive=="TRUE") %>% filter(ageInStage==0)
  
  #Filter out those born after year i
  settled2017 <- slice(settled, which(substr(settled$startTime,1,4)==proj_years[i]))
  
  #IDs of crab that survive to settlement born in year i
  tofind <- unique(settled2017$id)
  
  #Find these crab as zoea
  zoeaSurv <- filter(allresults, typeName == "Zooea1") %>% filter(id %in% tofind) %>% filter(ageInStage==0)
  
  #Find all zoea born in 2017
  allZooea <- filter(allresults, typeName == "Zooea1") 
  Zooea2017 <- slice(allZooea, which(substr(allZooea$startTime,1,4)==proj_years[i])) %>% filter(ageInStage==0)
  
  #Cumulative mortality proportion is the number of zoea that survived to settle over number born
  recCumMort <- 1-(nrow(zoeaSurv)/nrow(Zooea2017))
  totspawn <- nrow(Zooea2017)
  
  #Initialize objects
  settletime <- settlelat <- settlelong <- startlat <- startlong <- num <-  rep(0,nrow(zoeaSurv))
  probByLat <- sumProb <- zooea1SettlementMap <- vector("list")
  yrSettled <- unique(substr(settled2017$time,1,4))
  for(k in 1:length(yrSettled)){
    probByLat[[k]]<- array("",dim=c(length(longGrid),length(latGrid),nrow(zoeaSurv)))
  }
  
  
  #Loop to get time or settlement and path of each settled crab
  for(j in 1:length(tofind)){
    #For each crab, find it in the settled matrix
    thiscrab <- filter(settled2017, id==tofind[j])
    
    #The time it settled is the first time it shows as an immature crab
    settletime[j] <- thiscrab[1,"time"]
     yr<- as.integer(substr(thiscrab[1,"time"], 1,4))
     yrind <- which(yr==yrSettled)
     
     #Get location where it settled
     settlelong[j] <- thiscrab[1,"horizPos1"]
     settlelat[j] <- thiscrab[1,"horizPos2"]
     
     #Get how many crabs survived to settlement
     num[j] <- thiscrab[1,"number"]
     
     #Find where this crab started in the Zooea matrix
     startlong[j] <- filter(Zooea2017,id==thiscrab[1,'id'])$horizPos1
     startlat[j] <- filter(Zooea2017,id==thiscrab[1,'id'])$horizPos2
     
     #If the crab settled inside the grid track it, otherwise add it to mortality
     if((settlelong[j]>min(longGrid))&&(settlelong[j]<max(longGrid))&&(settlelat[j]>min(latGrid))&&(settlelat[j]<max(latGrid))){
      y <- min(last(which(startlong[j]>longGrid))+1,40)
  
      x <- first(which(startlat[j]>latGrid))
      y2 <- last(which(settlelong[j]>longGrid))+1
      x2 <- first(which(settlelat[j]>latGrid))

      probByLat[[yrind]][x,y,j] <- paste(x2,y2, sep=",")
    } else{
      recCumMort = recCumMort + (1/totspawn)
    }
  }
  
  #Across all of the crab that settled in a year, combine each into single grid cells based on where they started/ended up
  for(k in 1:length(yrSettled)){
    sumProb[[k]] <- probByLat[[k]][,,1]
    for(x in 1:(dim(probByLat[[k]])[3])){
      tmp <- probByLat[[k]][,,x][probByLat[[k]][,,x]!=""]
      if(length(tmp)>0){

        sumProb[[k]][probByLat[[k]][,,x]!=""]<-paste(sumProb[[k]][probByLat[[k]][,,x]!=""], probByLat[[k]][,,x][probByLat[[k]][,,x]!=""])
      }
    }
    
    #Transform grid cell full of strings into probabilities
    zooea1SettlementMap[[k]] <- array(dim=c(40,40,1600))
    for(l in 1:(dim(zooea1SettlementMap[[k]])[3])){
      if(sumProb[[k]][l]!=""){
        #If grid cell isn't emply, split into each entry
        transitionList<- unlist(strsplit(sumProb[[k]][l], " "))
        transitionList <- transitionList[transitionList!=""]
        map <- table(transitionList)
        coords<- as.numeric(unlist(strsplit(names(map),",")))
        if(length(coords)==2){
          #Put number into map if there's only one settlement
          zooea1SettlementMap[[k]][coords[1], coords[2],l] <- 1
        } else if(length(coords)>2){
          coords <- matrix(coords, ncol=2, byrow=T) 
          #If more than one settled in that grid, divide by the sum
          for(j in 1:nrow(coords)){
            zooea1SettlementMap[[k]][coords[j,1], coords[j,2],l] <- map[j]/sum(map)
          }
        }
      }
    }
  }
  
  #Write output
  for(l in 1:1600){
    write.table(zooea1SettlementMap[[k]][,,l],
                file=paste("Spawned",proj_years[i], "Settled",yrSettled[k],l,".csv", sep=""),
                col.names=as.character(round(longGrid,2)), 
                row.names=as.character(round(latGrid,2)), sep=",")
  }
  
  return(list("mort"=recCumMort, "time"=settletime))
  
}
