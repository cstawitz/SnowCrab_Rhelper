#Setup data frame to be plotted
mapsSetup<- function(df){
  df$posLongs <- df$horizPos1
  df$posLongs[df$horizPos1<0] <- ((df$horizPos1[df$horizPos1<0])+360)
  
  df$time <- as.POSIXct(df$time)
  df$year <- as.numeric(substr(df$time,1,4))
  return(df)
}

doGrid <- function(df){
  longs <-c(min(df$posLongs), max(df$posLongs))
  lats <- c(min(df$horizPos2), max(df$horizPos2))
  grid<-makeGrid(x=seq(longs[1],longs[2],by=.4),y=seq(lats[1],lats[2],by=.4), projection="LL")
  return(grid)
}


makeSizeByYear <- function(filename, df){
  pdf(filename)
  par(mar=c(0,2,0,2), mfrow=c(3,2))
  require(dplyr)
  brks <- calcBreaks(min(df$size), max(df$size))
  clrs <- brewer.pal(9,"PuRd")
  for(i in unique(df$year)){
    dframe <- filter(df,year==i)
    femEvent <- data.frame(X=dframe$posLongs, Y=dframe$horizPos2, EID = seq(1,nrow(dframe)),
                           Weight = dframe$weight, Size=dframe$size)
    require(RColorBrewer)
    events <- as.EventData(x=femEvent, projection="LL")
    locData<-findCells(events,grid)
    events$Z <- events$Size
    pdata <- combineEvents(events,locData, FUN=mean)

    pdata <- makeProps(pdata, brks, "col", clrs)
    plotMap(worldLLhigh, xlim=c(minLong,maxLong), ylim=c(minLat,maxLat))
    text(x=185,y=62, i)
    addPolys(grid, polyProps=pdata)
  }
  dev.off()
  system(paste("open", filename))
}

calcBreaks <- function(minimum, maximum){
  allSizes <- seq(minimum, maximum)
  interval <- round(length(allSizes)/10,0)
  colorbreaks <- rep(minimum, 9) + c(0,interval, interval*2, interval*3, interval*4, interval*5,interval*6, interval*7,
                                 interval*8)
  return(colorbreaks)
}

makeTimeToInstar<-function(filename, dataframe){
  pdf(filename)
  par(mar=c(0,2,0,2), mfrow=c(3,2))
  require(dplyr)
  for(i in unique(dataframe$year)){

    df <- filter(dataframe,year==i) %>% filter(ageInStage==0) %>% filter(instar==1)
    if(nrow(df)>0){
      malEvent <- data.frame(X=df$posLongs, Y=df$horizPos2, EID = seq(1,nrow(df)),
                             Weight = df$weight, Size=df$size,
                             Age = df$age)
      events <- as.EventData(x=malEvent, projection="LL")
      locData<-findCells(events,grid)
      events$Z <- events$Age
      pdata <- combineEvents(events,locData, FUN=mean)
      clrs <- brewer.pal(8,"PuBu")
 
      brks <- c(120,140,160,180,200,220,240,260,280)
      pdata <- makeProps(pdata, brks, "col", clrs)
      plotMap(worldLLhigh, xlim=c(minLong,maxLong), ylim=c(minLat,maxLat))
      text(x=185,y=62, i)
      addPolys(grid, polyProps=pdata)
    }
  }
  dev.off()
  system(paste("open", filename))
}

growthPerYear <- function(filename, dataframe){
  pdf(filename)
  par(mar=c(0,2,0,2), mfrow=c(2,2),mar=c(2,2,2,2))
  require(dplyr)
  yrs <- unique(dataframe$year)
  clrs <- brewer.pal(length(yrs),"Spectral")
  summaryTime<- bigresults %>% group_by(instar, year, ageInStage) %>% summarise(meanwt = mean(weight), meansize = mean(size))
  for(i in 1:9){
    eachCrab <- filter(summaryTime, instar==i)
    plot(NA, xlim=c(min(eachCrab$ageInStage),max(eachCrab$ageInStage)), ylim=c(min(eachCrab$meanwt),max(eachCrab$meanwt)))

    points(meanwt~ageInStage, data=eachCrab, col=paste(clrs[eachCrab$year-2005],"90", sep=''), pch=20)
    if(i %in% c(1,5,9)){legend("topleft",as.character(2006:2012), col=clrs, pch=19)}
  }

  dev.off()
  system(paste("open", filename))
}


sizePerYear <- function(filename, dataframe){
  pdf(filename)
  par(mar=c(0,2,0,2), mfrow=c(1,1),mar=c(2,2,2,2))
  require(dplyr)
  yrs <- unique(dataframe$year)[3:7]
  clrs <- brewer.pal(length(yrs),"Spectral")

  plot(NA, xlim=c(0,3500), ylim=c(3,45))
  timeByYear<- bigresults %>% group_by(year, instar) %>% summarise('timV'=mean(max(ageInInstar)), 'size'=mean(size))
  for(i in 1:length(yrs)){
  df <- filter(timeByYear,year==yrs[i], size<45)
  yval<-rep(df$size, df$timV)
  xval <- 1:sum(df$timV)
  lines(yval~xval, col=clrs[i], lwd=2)
  }
  abline(v=365)
  abline(v=730)
  abline(v=1095)
  abline(v=1460)
  abline(v=1825)
  abline(v=2190)
    legend("topleft",as.character(2008:2012), col=clrs, lwd=2)  
  dev.off()
  system(paste("open", filename))
}

makeTempPlot<- function(filename, dataframe){
  pdf(filename)
  par(mar=c(0,2,0,2), mfrow=c(3,2))
  require(dplyr)
  for(i in unique(dataframe$year)){
    
    df <- filter(dataframe,year==i)
    if(nrow(df)>0){
      malEvent <- data.frame(X=df$posLongs, Y=df$horizPos2, EID = seq(1,nrow(df)),
                             Temp = df$temperature)
      events <- as.EventData(x=malEvent, projection="LL")
      locData<-findCells(events,grid)
      events$Z <- events$Temp
      pdata <- combineEvents(events,locData, FUN=mean)
      clrs <- brewer.pal(9,"YlOrRd")

      brks <- c(-2,0,2,4,6,8,10,12, 14, 20)
      pdata <- makeProps(pdata, brks, "col", clrs)
      plotMap(worldLLhigh, xlim=c(minLong,maxLong), ylim=c(minLat,maxLat))
      text(x=185,y=62, i)
      addPolys(grid, polyProps=pdata)
    }
  }
  dev.off()
  system(paste("open", filename))
}

plotDepthByStage <- function(filename, dataframe){
  pdf(filename)
  par(mar=c(2,2,2,2), mfrow=c(3,2))
  lifestages <- unique(dataframe$typeName)
  for(i in 1:length(lifestages)){
    if(lifestages[i] %in% c("Zooea1", "Zooea2", "Megalopa")){
      maxD = 50
    } else{
      maxD = 500
    }
    df <- filter(dataframe, typeName==lifestages[i])
    plot(vertPos~age, data=df, ylim=c(0,maxD), col="gray")
    text(x=350, y=40, lifestages[i])
  }
  dev.off()
  system(paste("open", filename))
}

makeAllThePlots <- function(df, filenamesvector){
  
  
  #Make plot to check depth distribution across 
  plotDepthByStage(paste0(filenamesvector[1],".pdf"),df)
  df <- mapsSetup(df)
  grid <- doGrid(df)
  matureMales <- filter(df, typeName=="AdultMale")
  matureFemales <- filter(df, typeName=="AdultFemale")
  prefSizeMales <- filter(df, size>=101)
  immatureFemales <- filter(df, typeName %in% c("ImmatureFemale","AdolescentFemale"))
  immatureMales <- filter(df, typeName %in% c("ImmatureMale","AdolescentMale"))
  
  #Make mature male size by year
  makeSizeByYear(paste0(filenamesvector[2],".pdf"),matureMales)
  #Map preferred size males by year
  makeSizeByYear(paste0(filenamesvector[3],".pdf"),prefSizeMales)
  #Mature females by year
  makeSizeByYear(paste0(filenamesvector[4],".pdf"),matureFemales)
  #Male growth across years
  growthPerYear(paste0(filenamesvector[5],".pdf"),immatureMales)
  #Female growth across years
  growthPerYear(paste0(filenamesvector[6],".pdf"),immatureFemales)
  
}

makeStartDistByYear <- function(cpueData, years){
  returnmat <- array(0,dim=c(40,40,length(years)))
  latGrid <- seq(51.5,70, length.out=40)
  #This does not take into account differences in lines 
  longGrid <- seq(-179,-155,length.out=40)
  for(i in 1:length(years)){
    for(k in 1:39){
      for(j in 1:39){
        require(dplyr)
        df<- filter(cpueData, YEAR==years[i])
        latInd <- which((df$LATITUDE>latGrid[k])&(df$LATITUDE<latGrid[k+1]))
        longInd <- which((df$LONGITUDE>longGrid[j])&(df$LONGITUDE<longGrid[j+1]))
        if((length(latInd)>0)&(length(longInd)>0)){
          match.rows<- intersect(latInd,longInd)
          if(length(match.rows)==1){
            returnmat[k,j,i] <- df$numCPUE[match.rows]*ScaleUpFactor[k]
          } else{
            if(length(match.rows)>1){
              returnmat[k,j,i] <- mean(df$numCPUE[match.rows])*ScaleUpFactor[k]
            }
          }

        }
        
      }
    }
  }
  return(returnmat)
}


getLongGrid <- function(latGrid){
latToNM <- 69.172
NMLat <- latToNM*(max(latGrid)-min(latGrid))
longToNM <- function(x) return(cos(x)*latToNM*.8689)
NMtolong <- function(x) return(acos(x/(latToNM*.8689)))
deg2rad <- function(deg) {(deg * pi) / (180)}
longGrid <- matrix(ncol=80, nrow=length(latGrid))
nmPerLong <- rep(0, length(latGrid))
for(i in 1:length(latGrid)){
  #Calculate how many nm there are in each longitude
  nmPerLong[i] <- longToNM(deg2rad(latGrid[i]))
  #how many longitude units to get 600 nm
  longPer600 <- 600/nmPerLong[i]
  longGrid[i,] <- seq(-179,-179+longPer600,length.out=80)
}
ScaleUpFactor <- nmPerLong*69.172
return(longGrid)
}

loadResults <-function(wd, f){
setwd(paste0("C:\\Users\\Christine.Stawitz\\Documents\\SnowCrab_InputFiles\\",wd))
  col.call<- cols(
       mpg = col_double(),
       cyl = col_double()
    )
results <- read_csv_chunked("Results.csv", DataFrameCallback$new(f))
return(results)}