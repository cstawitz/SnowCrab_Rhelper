#Set relevant paths
surveydatadir <- "C:\\Users\\Christine.Stawitz\\Documents\\SurveyData"
surveydat <- read.csv(file.path(surveydatadir,"SnowCrab_HaulData.csv"))
surveysummarydir <- "C:\\Users\\Christine.Stawitz\\Dropbox\\Postdoc\\SurveyDataComparison"
snowcrabdir <- "C:\\Users\\Christine.Stawitz\\Documents\\SnowCrab_InputFiles"

#Read in female crab data
crabsummary <- vector("list")
yrs <- seq(2006,2011)
for(i in 1:length(yrs)){
  crabsummary[[i]]<-read.csv(file.path(surveysummarydir,paste("FemaleCrab", yrs[i], ".csv")))
}

##############################
# A bunch of scripts that look at DisMELS output and see if they make sense for each life stage

############################################
# immature females
immat_female<- readInLifeStage(snowcrabdir, "ImmatureFemale", numAttributes=27, "Results.csv")
names(immat_female)[24:27] <- c("shellthickness", "temperature", "salinity","pH")
#plot(size~age, immat_female)
require(RColorBrewer)

#Check growth curve for immature females
plot(weight~ageInStage, immat_female, col=colorMat[instar], pch=19)

#Check weight-at-length curve for immature females
plot(weight~size, immat_female, col=colorMat[instar], pch=19)
points(WEIGHT~WIDTH, filter(surveydat,SEX==2))
points(x=c(1,out[[2]]), y=predwt1, col="red", pch=19)
legend("bottomright",as.character(1:14),col=colorMat, pch=19)

#Check the number of females in each instar
hist(immat_female$instar)
immat_female$number <- as.numeric(as.character(immat_female$number))

#Check population growth over time
plot(number~age, immat_female)


require(dplyr)
group_by(immat_female,instar) %>%
  summarise(max(ageInInstar))

###############################################
# adolescent females
adol_female<- readInLifeStage(snowcrabdir, "AdolescentFemale", numAttributes=27)
names(adol_female)[24:27] <- c("shellthickness", "temperature", "salinity","pH")
plot(size~age, adol_female)
plot(weight~ageInStage, adol_female, col=colorMat[instar])
hist(adol_female$instar)
require(dplyr)
group_by(adol_female,instar) %>%
  summarise(max(ageInInstar))


##################################################
# adol male
adol_male<- readInLifeStage(snowcrabdir, "AdolescentMale", numAttributes=27)
names(adol_male)[24:27] <- c("shellthickness", "temperature", "salinity","pH")
adol_male$ageInStage <- adol_male$ageInStage+805
plot(weight~ageInStage, adol_male, col=colorMat[instar])

points(x=c(1,out[[2]]), y=predwt1, col="red", pch=19)
#################################################
#
zooea1 <- readInLifeStage(snowcrabdir, "Zooea1", numAttributes=24)
names(zooea1)
plot(weight~ageInStage,data=zooea1)
plot(number~age, data=zooea1)
zooea1$number

#################################################
#
zooea2 <- readInLifeStage(snowcrabdir, "Zooea2", numAttributes=24)
names(zooea2)
plot(weight~ageInStage, data=zooea2)
plot(number~ageInStage, data=zooea2)
#################################################
#
megalopa <- readInLifeStage(snowcrabdir, "Megalopa", numAttributes=24)
names(megalopa)
plot(number~ageInStage, data=megalopa)
megalopa$weight <- as.numeric(as.character(megalopa$weight))
plot(weight~ageInStage, data=megalopa)
 max(megalopa$weight)

 ######################################################
 adul_fem<- readInLifeStage(snowcrabdir, "AdultFemale", numAttributes=27, "Results (2).csv")
 names(adul_fem)[24:27] <- c("shellthickness", "temperature", "salinity","pH")
 
 immat_fem <- rbind(immat_female, adol_female)
 junJulyImm <- filter(immat_fem, as.numeric(substr(time,6,7))>6) %>%
   filter(as.numeric(substr(time,6,7))<8)
 
 ## Look at where females are in the summer
 junJulyFem <- filter(adul_fem, as.numeric(substr(time,6,7))>6) %>%
   filter(as.numeric(substr(time, 6,7))<8)
 junJulyFem2006 <- filter(junJulyFem, substr(time,1,4)=="2006")
 unique(junJulyFem2006$horizPos1)
 unique(junJulyFem2006$horizPos2)
 
 plot(weight~age, data=adul_fem)
 plot(pnorm(dat_pred,meanPre[5],sdPre[5])~dat_pred)
 
 ######################################################
 #Look at results between MIROC 4.5 and 8.5 scenarios
 adul_mal4.5 <- readInLifeStage(snowcrabdir, "AdultMale", numAttributes=27,"Results45.csv")
 adul_mal8.5 <- readInLifeStage(snowcrabdir, "AdultMale", numAttributes=27,"Results (2).csv")
  names(adul_mal)[24:27] <- c("shellthickness", "temperature", "salinity", "pH")
 adul_mal$ageInStage <- adul_mal$ageInStage + 1442
 plot(weight~ageInStage, adul_mal, col=colorMat[instar])
 points(x=c(1,out[[2]]), y=predwt1, col="red", pch=19)
 
 
########################################################
#Check maturity
 #These values taken from Burmeister and Sainte-Marie
 # Predicted mean and sd of post-molt chela width distribution at 5 degrees C
 sdmtemp <- standardDev(5, .20,.06,4,4.85)
 meanmtemp<- meanCW(5, .20,.06,4,4.85)
 



sdftemp <- standardDev(5, .16,.12,3.81,4.30)
meanftemp<- meanCW(5, .16,.12,3.81,4.30)

checkMaturity(adol_female, adul_fem, moltIncrementModF, sdftemp, meanftemp)
  

#Make plots!
install.packages("PBSmapping")
 require(PBSmapping)
data(worldLLhigh)
minLat <- 54
maxLat <- 65
minLong <- 164
maxLong <- 205
############################
#Plot meansize by area
if(length(immat_female$posLongs)==0){
  immat_female <- mapsSetup(immat_female)
}
grid <- doGrid(immat_female)
femEvent <- data.frame(X=immat_female$posLongs, Y=immat_female$horizPos2, EID = seq(1,nrow(immat_female)),
                       Weight = immat_female$weight, Size=immat_female$size)
events <- as.EventData(x=femEvent, projection="LL")
locData<-findCells(events,grid)
events$Z <- events$Size
pdata <- combineEvents(events,locData, FUN=mean)
clrs <- brewer.pal(5,"YlGn")
brks <- c(0,10,20,40,90,150)
pdata <- makeProps(pdata, brks, "col", clrs)
par(mar=c(2,2,2,2))
plotMap(worldLLhigh, xlim=c(minLong,maxLong), ylim=c(minLat,maxLat))
addPolys(grid, polyProps=pdata)

#By Year
if(length(immat_female$posLongs)==0){
  immat_female <- mapsSetup(immat_female)
}


##########
# Make plot for mean time to juvenile crab
if(length(immat_male$posLongs)==0){
  immat_male <- mapsSetup(immat_male)
}
grid <- doGrid(immat_male)
malEvent <- data.frame(X=immat_male$posLongs, Y=immat_male$horizPos2, EID = seq(1,nrow(immat_male)),
                       Weight = immat_male$weight, Size=immat_male$size,
                       Instar = immat_male$instar, Age = immat_male$age)
malEvent<- filter(malEvent, Instar==1)
grid<-makeGrid(x=seq(longs[1],longs[2],by=.5),y=seq(lats[1],lats[2],by=.5), projection="LL")
events <- as.EventData(x=malEvent, projection="LL")
locData<-findCells(events,grid)
events$Z <- events$Age
pdata <- combineEvents(events,locData, FUN=mean)
clrs <- brewer.pal(5,"YlGn")
brks <- c(0,50,100,200,400,800)
pdata <- makeProps(pdata, brks, "col", clrs)
par(mar=c(2,2,2,2))
plotMap(worldLLhigh, xlim=c(minLong,maxLong), ylim=c(minLat,maxLat))
addPolys(grid, polyProps=pdata)

#By Year
immat_male$time <- as.POSIXct(immat_male$time)
immat_male$year <- as.numeric(substr(immat_male$time,1,4))

#########Plot mean temperature experienced by crab in DisMELS by year
pdf("MeanTempByYear.pdf")
par(mar=c(0,2,0,2), mfrow=c(3,2))
require(dplyr)
for(i in unique(immat_male$year)){
  df <- filter(immat_male,year==i) 
  if(nrow(df)>0){
    malEvent <- data.frame(X=df$posLongs, Y=df$horizPos2, EID = seq(1,nrow(df)),
                           Temp = df$temperature)
    
    events <- as.EventData(x=malEvent, projection="LL")
    locData<-findCells(events,grid)
    events$Z <- events$Temp
    pdata <- combineEvents(events,locData, FUN=mean)
    clrs <- brewer.pal(6,"YlOrRd")
    brks <- c(-2,0,2,4,6,8,10)
    pdata <- makeProps(pdata, brks, "col", clrs)
    plotMap(worldLLhigh, xlim=c(minLong,maxLong), ylim=c(minLat,maxLat))
    text(x=185,y=55, i)
    addPolys(grid, polyProps=pdata)
  }
}
dev.off()
system("open MeanTempByYear.pdf")


#Read in all results
colNames <- c(names(resultsBig)[1:17], "ageInInstar", "number", "instar", "size", "weight", "shellcondition",
  names(resultsBig)[21:24])
resultsBig<- read.table("Results.csv", sep=",", skip=1,col.names = colNames, fill=TRUE)

#Get settled individuals
settled <- filter(resultsBig, typeName %in% c("AdultMale", "AdultFemale", "AdolescentFemale", "AdolescentMale",
                                              "ImmatureMale", "ImmatureFemale"))

#Get larvae
larvae <- filter(resultsBig, typeName %in% c("Zooea1", "Zooea2", "Megalopa"))
names(larvae) <- c(names(larvae)[-c(18,20,21,23)], "Blank1", "Blank2", "Blank3", "Blank4")
larvalcount <- larvae %>% group_by(typeName) %>% summarise(num = sum(number,na.rm=TRUE))
settcount <- settled %>% group_by(typeName) %>% summarise(num=sum(number))
countByLifeStage<-rbind(larvalcount, settcount)
countByLifeStage$typeName <- as.character(countByLifeStage$typeName)
sorted <- countByLifeStage[c(order(countByLifeStage$num)),]

#Plot number that settle by count
plot(sorted$num, log="y", axes=NA)
axis(2)
axis(1, labels=sorted$typeName, at=seq(1,9))


###Calculate centroid differences between output
centroidCalc(bigresults)
centroidCalc(adul_mal)
require(readr)
colspec <- cols(
  .default = col_integer(),
  typeName = col_character(),
  time = col_datetime(format = ""),
  horizPos1 = col_double(),
  horizPos2 = col_double(),
  vertPos = col_double(),
  number = col_double(),
  size = col_double(),
  weight = col_double(),
  temperature = col_double()
)


setwd(snowcrabdir)
f <- function(x, pos) subset(x, typeName=="AdultFemale")
bigresults <- read_csv_chunked(file="Results.csv", DataFrameCallback$new(f),col_types = colspec)
save(bigresults, file="AdultFem.RData")
bigresults <- mapsSetup(bigresults) 
centroidCalc(bigresults)

grid <- doGrid(bigresults)
makeSizeByYear("FemAdultSizeYr.pdf", bigresults)
makeTimeToInstar("FemAdultTimeTo1stInstar.pdf", bigresults)
rm(bigresults)

f <- function(x, pos) subset(x, typeName=="ImmatureMale")
bigresults <- read_csv_chunked(file="Results.csv", DataFrameCallback$new(f),col_types = colspec)

str(bigresults)
save(bigresults, file="ImmatMale.RData")
load("ImmatMale.RData")
bigresults <- mapsSetup(bigresults)
grid <- doGrid(bigresults)
makeTimeToInstar("ImmatMalTimeTo1stInstar.pdf", bigresults)

require(RColorBrewer)
growthPerYear("ImmatMalGrowth.pdf", bigresults)
sizePerYear("SizesByYear.pdf", bigresults)

####################################################
# Create plots of preferred size males by the fishery

f <- function(x, pos) subset(x, typeName=="AdolescentMale")
adol_male <- read_csv_chunked(file="Results.csv", DataFrameCallback$new(f),col_types = colspec)
f <- function(x, pos) subset(x, typeName=="AdultMale")
adult <- read_csv_chunked(file="Results.csv", DataFrameCallback$new(f),col_types = colspec)
oldermales <- rbind(adol_male, adult)
oldermales$month <- as.numeric(substr(oldermales$time,6,7))
bigmales <- filter(oldermales, size>100) %>% filter( month %in% c(10:12, 1:5))
bigmales <- mapsSetup(bigmales)
grid <- doGrid(bigmales)
makeSizeByYear("PreferredSizeMales.pdf", bigmales)

#########################################
# Make all the plots

allresults <- read.csv("Results.csv")
filenames <- c("DepthAcrossStage", "MatureMalesByYear", "PreferredSizeMalesByYear",
               "MatureFemalesByYear", "MaleGrowthAcrossYears","FemaleGrowthAcrossYears" )
if(length(zooea1_Start$posLongs)==0){
  zooea1_Start <- mapsSetup(zooea1_Start)
}
grid <- doGrid(zooea1_Start)
if(length(settled$posLongs)==0){
settled <- mapsSetup(settled)
}



############################################
# Make plots of settlement by year
########################################
settleEvent <- data.frame(X=settled$posLongs, Y=settled$horizPos2, EID = seq(1,nrow(settled)),
                       ID = settled$id, Time=settled$time, Num = settled$number)

years <- seq(2006,2011, length.out=6)
settleByYear <- vector("list")
clrs <- brewer.pal(5,"YlGn")
brks <- c(0,.1,.5,1,1.5,2)
pdf("SettleByYear.pdf")
for(i in 1:6){
  preTime <- as.POSIXct(strptime(paste(as.character(years[i]), "-09-01", sep=""),"%Y-%m-%d"))
  postTime <- as.POSIXct(strptime(paste(as.character(years[i]+1), "-04-01", sep=""), "%Y-%m-%d"))
settleByYear[[i]] <- filter(settleEvent, 
       (Time>preTime)&(postTime>Time))

events <- as.EventData(x=settleByYear[[i]], projection="LL")
locData<-findCells(events,grid)
events$Z <- events$Num
pdata <- combineEvents(events,locData, FUN=sum)

pdata <- makeProps(pdata, brks, "col", clrs)
par(mar=c(2,2,2,2))
plotMap(worldLLhigh, xlim=c(minLong,maxLong), ylim=c(minLat,maxLat))
text(x=185,y=62, paste("Settlement October",years[i], "- April", years[i]+1))
addPolys(grid, polyProps=pdata)
}
dev.off()



#Get Starting distributions from survey data
setwd("//akc0ss-n086/REFM_Users/Christine.Stawitz/My Documents/SurveyData")
data_csv <- read.csv("cpue.csv")
require(dplyr)
onlypos <- filter(data_csv, numCPUE>0)
immFemNew <- filter(data_csv, 
                    (SEX=="FEMALE")&(MATURITY=="IMMATURE")&(SHELL_CONDITION =="NEW_SHELL") & (numCPUE>0))
immFemOld <- filter(data_csv, 
                    (SEX=="FEMALE")&(MATURITY=="IMMATURE")&(SHELL_CONDITION =="OLD_SHELL")& (numCPUE>0))
matFemNew <- filter(data_csv, 
                    (SEX=="FEMALE")&(MATURITY=="MATURE")&(SHELL_CONDITION =="NEW_SHELL")& (numCPUE>0))
matFemOld <- filter(data_csv, 
                    (SEX=="FEMALE")&(MATURITY=="MATURE")&(SHELL_CONDITION =="OLD_SHELL")& (numCPUE>0))
MalNew <- filter(data_csv, 
                 (SEX=="MALE")&(SHELL_CONDITION =="NEW_SHELL"))
MalOld <- filter(data_csv, 
                 (SEX=="MALE")&(SHELL_CONDITION =="OLD_SHELL"))

par(mfrow=c(2,2))
plot(numCPUE~SIZE, data=immFemNew, main="Immature new shell Females", xlim=c(25,75))
plot(numCPUE~SIZE, data=immFemOld, main="Immature old shell Females", xlim=c(25,75))
plot(numCPUE~SIZE, data=matFemNew, main="Mature new shell Females", xlim=c(25,75))
plot(numCPUE~SIZE, data=matFemOld, main="Mature old shell Females", xlim=c(25,75))

years <- seq(1981,2017)
library(ggplot2)


makePlotBySize(immFemNew, "immNewFem.pdf")
makePlotBySize(immFemOld, "immOldFem.pdf")
makePlotBySize(matFemNew, "matFemNew.pdf")
makePlotBySize(matFemOld, "matFemOld.pdf")
makePlotBySize(rbind(immFemNew,immFemOld), "AllImmFem.pdf")

immfemnewyr <- makePlotBySizeAndYr(immFemNew, "immNewFemByYear.pdf")
immfemoldyr <-makePlotBySizeAndYr(immFemOld, "immOldFemYr.pdf")
matfemnewyr <-makePlotBySizeAndYr(matFemNew, "matFemNewYr.pdf")
matfemoldyr <-makePlotBySizeAndYr(matFemOld, "matFemOldYr.pdf")
allImmfemYr <-makePlotBySizeAndYr(rbind(immFemNew,immFemOld), "AllImmFemYr.pdf")






library(RColorBrewer)
my.pal <- c("#9e4a7f", "#53b548",
           "#8948c5",
           "#99b534",
           "#b271eb",
           "#e3aa2d",
           "#5c73e1",
           "#bba73a",
           "#a944a9",
           "#7abe77",
           "#dd4cb7",
           "#37814c",
           "#e83889",
           "#4cbdaf",
           "#d64f22",
           "#62a0d8",
           "#e27e2e",
           "#7955b2",
           "#cf9239",
           "#e27be1",
           "#6a792d",
           "#bb2a85",
           "#c2ab6d",
           "#516eb6",
           "#d65041",
           "#bd94dd",
           "#92692f",
           "#e46db0",
           "#a44d28",
           "#895c9c",
           "#e38d66",
           "#db87b1",
           "#cd374b",
           "#db7d7d",
           "#d54975",
           "#9c4556")
my.pal <- paste0(my.pal,"80")
par(mfrow=c(1,2))

##Initialize to full length of data
years <- unique(c(matFemNew$YEAR, matFemOld$YEAR))
sizes <- sort(unique(c(matFemNew$SIZE, matFemOld$SIZE)))
res <- getCentroidDifference(matFemNew, matFemOld, years, sizes)

pdf("AllYearsAllSizesNewShell.pdf")
plot(NA,xlim=c(-178,-164), ylim=c(56,62), main="New - old shell")
arrows(x0=res$long[,,1], y0=res$lat[,,1], x1=res$long[,,2], y1=res$lat[,,2], col=my.pal, 
       xlab="Long", ylab="Lat", lwd=(1:15)/2)
legend("topright", as.character(seq(1981,2016)),col=my.pal, lwd=1, ncol=2)
dev.off()

#By size
for(j in 1:length(sizes)){

res <- getCentroidDifference(matFemNew, matFemOld, years, sizes[j])
pdf(paste("AllYrsNewShell",sizes[j],".pdf"))
plot(NA,xlim=c(-178,-164), ylim=c(56,62), main="New - old shell")
arrows(x0=res$long[,1], y0=res$lat[,1], x1=res$long[,2], y1=res$lat[,2], col=my.pal, 
       xlab="Long", ylab="Lat", lwd=1)
legend("topright", as.character(seq(1981,2016)),col=my.pal, lwd=1, ncol=2)
dev.off()
}


par(mfrow=c(1,1))
resOld <- getCentroidDifference(matFemOld, matFemOld, years, sizes)
plot(NA,xlim=c(-178,-162), ylim=c(55,62), main="Old - old shell", xlab="Long", ylab="Lat")
arrows(x0=resOld$long[,,1], y0=resOld$lat[,,1], x1=resOld$long[,,2], y1=resOld$lat[,,2], col=my.pal,
       lwd = (1:15)/2)
legend("topright", as.character(seq(1981,2016)),col=my.pal, lwd=1)


#By size
for(j in 1:length(sizes)){
  
  res <- getCentroidDifference(matFemOld, matFemOld, years, sizes[j])
  pdf(paste("AllYrsOldShell",sizes[j],".pdf"))
  plot(NA,xlim=c(-178,-164), ylim=c(56,62), main="New - old shell")
  arrows(x0=res$long[,1], y0=res$lat[,1], x1=res$long[,2], y1=res$lat[,2], col=my.pal, 
         xlab="Long", ylab="Lat", lwd=1)
  legend("topright", as.character(seq(1981,2016)),col=my.pal, lwd=1, ncol=2)
  dev.off()
}



distance <- GetDist(resOld)
oldMagnitudeBySize <- apply(distance$dist,2,FUN=mean, na.rm=T)
oldMagnitudeByYear <- apply(distance$dist,1,FUN=mean, na.rm=T)
distanceNew <- GetDist(res)
newMagnitudeBySize <- apply(distanceNew$dist,2,FUN=mean, na.rm=T)
newMagnitudeByYear <- apply(distanceNew$dist,1,FUN=mean, na.rm=T)

library(plotrix)
distance$dist[is.na(distance$dist)] <- rep(0, length(distance$dist[is.na(distance$dist)]))
distanceNew$dist[is.na(distanceNew$dist)] <- rep(0, length(distanceNew$dist[is.na(distanceNew$dist)]))
distanceNew$angle[,15] <- rep(0, 36)
pdf(paste("RadialPrimipara.pdf"))
for(i in 1:(length(years)-1)){

  
#radial.plot(distance$dist[i,],radial.pos=distance$angle[i,],main=paste("Movement multipara",years[i], years[i+1]),lwd=3,line.col=my.pal[1:15],
 #           labels = seq(40,60, by=5), label.pos = distance$angle[i,], radial.lim=c(0,8.51))
radial.plot(distanceNew$dist[i,],radial.pos=distanceNew$angle[i,],main=paste("Movement primipara",years[i], years[i+1]),line.col=my.pal[1:15],
            lwd=3,
            labels = seq(40,60, by=5), label.pos = distanceNew$angle[i,], radial.lim=c(0,8.51))

}
dev.off()
mean(distanceNew$angle[distanceNew$angle!=0])
mean(distance$angle[distanceNew$angle!=0])
mean(distanceNew$dist[distanceNew$angle!=0])
mean(distance$dist[distanceNew$angle!=0])



# reset everything
par(oldpar)



yearsofdat <- unique(data_csv$YEAR)
cpueByLatImmFemNew <- makeStartDistByYear(immFemNew, yearsofdat)
  cpueByLatImmFemOld <- makeStartDistByYear(immFemOld, yearsofdat)
  cpueByLatMatFemNew <- makeStartDistByYear(matFemNew, yearsofdat)
  cpueByLatMatFemOld <-makeStartDistByYear(matFemOld, yearsofdat)
  cpueByLatMalNew <-  makeStartDistByYear(MalNew, yearsofdat)
  cpueByLatMalOld <- makeStartDistByYear(MalOld, yearsofdat)
  
  IFN <- IFO <- MFN <- MFO <- MN <- MO <- rep(0, length(yearsofdat))
  IFN <- cpueByLatImmFemNew  
  IFO <- cpueByLatImmFemNew  
  MFN <- cpueByLatImmFemNew  
  MFO <- cpueByLatImmFemNew  
  MN <- cpueByLatImmFemNew  
  MO <- cpueByLatImmFemNew  
for(i in 1:length(yearsofdat)){
IFN[,,i]<-  cpueByLatImmFemNew[,,i]/sum(cpueByLatImmFemNew[,,i])
sum(IFN[,,i]>0)
IFO[,,i]<-  cpueByLatImmFemOld[,,i]/sum(cpueByLatImmFemOld[,,i])

sum(IFO[,,i]>0)
MFN[,,i]<-  cpueByLatMatFemNew[,,i]/sum(cpueByLatMatFemNew[,,i])

sum(MFN[,,i]>0)
MFO[,,i]<-  cpueByLatMatFemOld[,,i]/sum(cpueByLatMatFemOld[,,i])

sum(MFO[,,i]>0)
MN[,,i]<-  cpueByLatMalNew[,,i]/sum(cpueByLatMalNew[,,i])

sum(MN[,,i]>0)
MO[,,i]<-  cpueByLatMalOld[,,i]/sum(cpueByLatMalOld[,,i])
sum(MO[,,i]>0)
}

skipMoltF <-  apply(IFO[,,2:length(yearsofdat)]>0, c(3), sum)/
  (apply(IFN[,,1:(length(yearsofdat)-1)]>0, c(3), sum) +apply(MFN[,,1:(length(yearsofdat)-1)]>0, c(3), sum))

skipMoltM <- apply(MO[,,2:length(yearsofdat)]>0, c(3), sum)/(apply(MN[,,1:(length(yearsofdat)-1)]>0, c(3), sum))
alltemp <- read.csv("TempYearAvg.csv")
temp<-filter(read.csv("TempYearAvg.csv"), (AKFIN_SURVEY_YEAR>1980)&(AKFIN_SURVEY_YEAR<2017))
plot(skipMoltM~temp$AKFIN_SURVEY_YEAR, ylim=c(-0.2,2.5))
plot(skipMoltF~temp$AKFIN_SURVEY_YEAR, col="red")
cor(temp$average, skipMoltM)

ma <- function(x,n=3){stats::filter(x,rep(1/n,n), sides=1)}
movavg<- ma(alltemp, n=3)
skippedF <- apply(IFO[,,2:length(yearsofdat)]>0, c(3), sum)
pmodel <- glm(skipMoltF~temp$average, family=binomial(link="logit"))
pmodel <- glm(skipMoltF~movavg[7:42,3], family=binomial(link="logit"))

lotsskipmolt <- which(skipMoltF>.05)
temp[lotsskipmolt+1,2]
summary(pmodel)
plot(movavg[,3]~movavg[,2])
plot(temp$average~temp$AKFIN_SURVEY_YEAR, type="l")

for(i in 1:yearsofdat){
  
  rbind()
}



settlers <- getSettlers(zooea1, immat_male)
tofind <- unique(zooea1$id)

forid <- unique(bigsettled$id)


years <- 2017:2020
all <- unique(tofind$id)
livelarvaeID <- all[which(all %in% forid)]

immat_male<- readInLifeStage(snowcrabdir, "ImmatureMale", numAttributes=27, "Results45.csv")
zooea1<- readInLifeStage(snowcrabdir, "Zooea1", numAttributes=24, "Results45.csv")
zooea12 <-readInLifeStage(snowcrabdir, "Zooea1", numAttributes=24, "Results (2).csv")

pdf("SettlementScenario.pdf")
plot(NA, xlim=c(8,24), ylim=c(0,5), axes=F, xlab="", ylab="")
axis(2, labels=c(years), at=1:4, las=1)
mtext("Birth year", 2, line=3)
axis(1, labels=c("Oct", "March", "Aug"), at=c(10,15,20))
mtext("Settlement month", 1, line=3)

for(i in 1:length(years)){
  tofind <- filter(zooea1, ageInStage==0)
  settlers <- getSettlers(tofind, immat_male)
  plotLines(tofind, settlers, i, cols["RCP4.5"], i)
}

for(i in 1:length(years)){
  tofind <- filter(zooea12, ageInStage==0)
  settlers <- getSettlers(tofind, immat_male2)
  plotLines(tofind, settlers, i,cols["RCP8.5"], i+.2)
}

dev.off()

i<-1

#Compare to batch runs
for(i in 1:length(years)){
  setwd(file.path("C:\\Users\\Christine.Stawitz\\Documents\\SnowCrab_InputFiles\\OldBatch",years[i]))
  immat_male<- readInLifeStage(getwd(), "ImmatureMale", numAttributes=27)
  zooea1<- readInLifeStage(getwd(), "Zooea1", numAttributes=24)
  tofind <- filter(zooea1, ageInStage==0)
  settlers <- getSettlers(tofind, immat_male)

  plotLines(tofind, settlers, i, "green", i+.4)
}
legend("bottom", c("Uniform egg distribution", "Using crab distribution"), pch=19, lty=1, col=c("purple", "green"), border=NA)


