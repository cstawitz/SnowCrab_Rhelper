#This script recreates the plots from the ECCWO presentation 


#Read in data for each LHS
snowcrabdir<- "F:\\SnowCrab_InputFiles"
#Read in data for each LHS
snowcrabdir<- "C:\\Users\\Christine.Stawitz\\Documents\\Snowcrab_Inputfiles"
#For laptop
snowcrabdir <- "C:\\Users\\chris\\Dropbox\\Postdoc\\Code"
immat_male<- readInLifeStage(snowcrabdir, "ImmatureMale", numAttributes=27,"Results.csv")
immat_male <- mapsSetup(immat_male)
immat_male2 <- mapsSetup(immat_male2)
loc45 <- sample_n(immat_male, 1000)
loc85 <- sample_n(immat_male2,1000)
with(loc45, plot(rev(posLongs), horizPos2, col=cols["RCP4.5"]))
with(loc85, points(rev(posLongs), horizPos2, col=cols["RCP8.5"]))
pdf("Latitudes.pdf")
boxplot(immat_male$horizPos2, immat_male2$horizPos2, 
        names=c("RCP 4.5", "RCP 8.5"), 
        col=cols[c("RCP4.5","RCP8.5")],
        ylab="Latitude", las=1)
dev.off()

summ <- immat_male %>%
  group_by(instar, year) %>%
  summarize(meantmp = mean(temperature),
            meanlat = mean(horizPos2),
            meandep = mean(vertPos))
summ2 <- immat_male2 %>%
  group_by(instar, year) %>%
  summarize(meantmp = mean(temperature),
            meanlat = mean(horizPos2),
            meandep = mean(vertPos))
which(summ$meantmp>summ2$meantmp)
plot(meantmp~instar, col=year,data=summ)
points(meantmp~instar, col=year, pch=3, data=summ2)
boxplot(immat_male$vertPos, immat_male2$vertPos)
names(immat_male)[24:27] <- c("shellthickness", "temperature", "salinity","pH")
pdf("ExperiencedTemps.pdf")
hist(immat_male$temperature, ylab="Frequency", 
     xlab="Temperature (d C)", 
     main="", col=paste(cols["RCP4.5"],"50",sep=""))
hist(immat_male2$temperature, col=paste(cols["RCP8.5"],"50",sep=""), add=T)
dev.off()

# immat_male$startTime <- as.POSIXct(immat_male$startTime)
immat_male$time <- as.POSIXct(immat_male$time)
hist(immat_male$time, breaks="months")
# round(immat_male$time-immat_male$startTime,0)==round(immat_male$age,0)
# par(mfrow=c(2,2))
# plot(size~age, immat_male)
require(RColorBrewer)
colorMat <- rep(brewer.pal(11,"Spectral"),2)
colorYear <- brewer.pal(11,"Spectral")
years <- unique(years)
immat_male2$year <- as.numeric(substr(immat_male2$time,1,4))
years <- unique(immat_male$year)




rcp45 <- df(immat_male2)

meantemp <- mean(immat_male$temperature)

plotOneCrab(9, meantemp)
cols <- c("RCP8.5"="#b94b75", "RCP4.5"="#50b47b", "DRPur"="#8177FF")
#melted <- reshape2::melt(toPlot, id.vars=c("age","year"))
require(ggplot2)
p <- ggplot(NULL) +
  geom_ribbon(data = filter(toPlot, year>2018),aes(x=age,ymin=lo, ymax=hi, group=as.factor(year),fill="RCP8.5"), alpha=.6) +
  geom_point(data = filter(toPlot, year>2018), aes(x=age, y=med, group=as.factor(year), color="RCP8.5"),alpha=.6) +
  scale_y_continuous("log(weight) [g]", trans="log", labels=round) +
  xlab("age (days)") +
  theme_classic() 

p + 
  geom_ribbon(data = filter(rcp45, year>2018), 
              aes(x=age,ymin=lo, ymax=hi, group=as.factor(year), fill="RCP4.5"),alpha=.6) +
  geom_point(data = filter(rcp45, year>2018), 
             aes(x=age, y=med, group=as.factor(year), color="RCP4.5"), alpha=.8) +
  scale_color_manual(name="Climate scenario", values=cols) +
  scale_fill_manual(name="Climate scenario", values=cols)


impd <- getIMPeriod(meantemp,siz)
wt<-rep(0,nrow(meang))
wt[1] <- aLW*siz[1]^bLW
for(i in 2:length(siz)){
  wt[i] <- wt[i-1]*exp(g[i]*impd[i])
}
meang <-data.frame(age=impd, w=siz^(g*impd))
meang$age <- 1:nrow(meang) %>% map_dbl(~sum(impd[1:.x]))
plot(x=meang$age, y=meang$w, log="y", type="l")



plot(weight~ageInStage, filter(immat_male, year==2021),
     col=1, ylim=c(0,10), cex=1.5)
legend("topleft",legend=years[2:3], col=colorYear[c(4,5)], pch=1, cex=2)
dev.off()
plot(weight~size, immat_male, col=colorMat[instar])
points(x=c(1,out[[2]]), y=predwt1, col="red", pch=19)
points(WEIGHT~WIDTH, filter(surveydat,SEX==1))
hist(immat_male$instar)
immat_male$number <- as.numeric(as.character(immat_male$number))
plot(number~ageInStage, immat_male)
m276 <- filter(immat_male, id==276)
m380 <- filter(immat_male, id==380)
plot(number~age, data=m276)
points(number~age, data=m380, col=2)