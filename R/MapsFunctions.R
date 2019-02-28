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
  for(i in unique(df$year)){
    dframe <- filter(df,year==i)
    femEvent <- data.frame(X=dframe$posLongs, Y=dframe$horizPos2, EID = seq(1,nrow(dframe)),
                           Weight = dframe$weight, Size=dframe$size)
    require(RColorBrewer)
    events <- as.EventData(x=femEvent, projection="LL")
    locData<-findCells(events,grid)
    events$Z <- events$Size
    pdata <- combineEvents(events,locData, FUN=mean)
    clrs <- brewer.pal(5,"YlGn")
    brks <- c(0,10,20,40,90,150)
    pdata <- makeProps(pdata, brks, "col", clrs)
    plotMap(worldLLhigh, xlim=c(minLong,maxLong), ylim=c(minLat,maxLat))
    text(x=185,y=62, i)
    addPolys(grid, polyProps=pdata)
  }
  dev.off()
  system(paste("open", filename))
}

makeTimeToInstar<-function(filename, dataframe){
  pdf(filename)
  par(mar=c(0,2,0,2), mfrow=c(3,2))
  require(dplyr)
  for(i in unique(dataframe$year)){

    df <- filter(dataframe,year==i) %>% filter(ageInStage==0)
    if(nrow(df)>0){
      malEvent <- data.frame(X=df$posLongs, Y=df$horizPos2, EID = seq(1,nrow(df)),
                             Weight = df$weight, Size=df$size,
                             Age = df$age)
      events <- as.EventData(x=malEvent, projection="LL")
      locData<-findCells(events,grid)
      events$Z <- events$Age
      pdata <- combineEvents(events,locData, FUN=mean)
      clrs <- brewer.pal(5,"YlGn")
      brks <- c(0,100,200,400,800)
      pdata <- makeProps(pdata, brks, "col", clrs)
      plotMap(worldLLhigh, xlim=c(minLong,maxLong), ylim=c(minLat,maxLat))
      text(x=185,y=62, i)
      addPolys(grid, polyProps=pdata)
    }
  }
  dev.off()
  system(paste("open", filename))
}
