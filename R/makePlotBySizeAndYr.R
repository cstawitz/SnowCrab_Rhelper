
makePlotBySizeAndYr <- function(theData, plotName){
  sizes <- sort(unique(theData$SIZE))
  years <- sort(unique(theData$YEAR))
  pdf(paste(plotName))
  bysizeyr <- vector("list")
  for(j in 1:length(years)){
    bysizeyr[[j]] <- vector("list")
    for(i in 1:length(sizes)){
      
      bysizeyr[[j]][[i]] <- filter(theData, SIZE==sizes[i]) %>% filter(YEAR==years[j])
      p<- ggplot(bysizeyr[[j]][[i]], aes(x=LONGITUDE, y=LATITUDE, size=numCPUE, fill=SIZE)) +
        geom_point() + ggtitle(paste(sizes[i], years[j])) + xlim(-180,-160) + ylim(54,63) +
        scale_size_continuous(breaks = c(1,100,1000,10000,100000,1000000), range=c(1,6))
      print(p)
    }
  }
  dev.off()
  return(bysizeyr)
}