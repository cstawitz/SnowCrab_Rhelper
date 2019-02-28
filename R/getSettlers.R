
getSettlers <- function(tofind,immat_male){
  years <- unique(substr(tofind$time,1,4))
  settlers <- filter(immat_male, instar==1) %>% filter(ageInInstar==0)
  return(settlers)
}