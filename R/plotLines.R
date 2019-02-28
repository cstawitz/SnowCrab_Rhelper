
plotLines <- function(tofind, settlers, i, colores, ylev){
  ids_yr <- filter(tofind, substr(time,1,4)==years[i])$id
  settle_that_yr <- filter(settlers, id %in% ids_yr)
  format_time <- as.POSIXlt(settle_that_yr$time)
  settleTimeInt <- format_time$year + format_time$mon/12 + format_time$mday/31-min(format_time$year)
  quants<- quantile(settleTimeInt, c(.025,.5,.975))*12
  points(x=quants[2], y=ylev, pch=19, col=colores)
  segments(x0=quants[1], y0=ylev, x1=quants[3], y1=ylev, col=colores)
}
