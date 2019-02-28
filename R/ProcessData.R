df <- function(dataset){ 
  require(purrr)
  placedf <- dataset %>% 
  split(list(.$ageInStage, .$year)) %>% 
  map_dfr(~as.data.frame(t(as.matrix(quantile(.$weight, c(.025, .5,.975), na.rm=TRUE, names=F)))))
toPlot <- cbind(expand.grid(unique(dataset$ageInStage), unique(dataset$year)),placedf)
names(toPlot)<- c("age","year","lo","med","hi")
return(toPlot)
}