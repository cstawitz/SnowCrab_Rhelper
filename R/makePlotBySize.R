makePlotBySize <- function(theData, plotName){
  sizes <- sort(unique(theData$SIZE))
  pdf(paste(plotName))
  for(i in 1:length(sizes)){
    
    p<- ggplot(filter(theData, SIZE==sizes[i]), aes(x=LONGITUDE, y=LATITUDE, size=numCPUE, fill=SIZE)) +
      geom_point() + ggtitle(sizes[i]) + xlim(-180,-160) + ylim(54,63) +
      scale_size_continuous(breaks = c(1,100,1000,10000,100000,1000000), range=c(1,6))
    print(p)
    
  }
  dev.off()
}
