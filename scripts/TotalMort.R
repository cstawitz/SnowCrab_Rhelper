getCumMort <- function(temperature, CWmat, Ml, Mcrab){
  totalMortz1 <- exp(-Ml*max(0,(229.5/(temperature-.63))))
  totalMortz2 <- exp(-Ml*max((530.5)/temperature,0))
  totalMortm <- exp(-Ml*(417.3)/(temperature+2.24))
  totalMortcrab <- rep(0, length(CWmat))
  for(i in 1:length(CWmat)){
    predday <- (1573*(1
      -exp(.07*CWmat[i])))/(temperature+4.71*(1-exp(.14*CWmat[i])))
    totalMortcrab[i] <- exp(-Mcrab*predday)
  }
  return(totalMortz1*totalMortz2*totalMortm*sum(totalMortcrab))
}

femValues<-unlist(lapply(-1:10, getCumMort, c(3.19, 5.086, 7.34,10.03,12.52,
                 15.87,20.39,26.45,34.617,45.597,
                 60.367), .0161, .0014))
 
malValues <- unlist(lapply(-1:10, getCumMort, c(3.19, 5.086, 7.34,10.03,14.31,
                                   19.55,25.98,33.86,43.53,55.38,
                                   69.91), .0161, .000959))
 temperature <- 5
 