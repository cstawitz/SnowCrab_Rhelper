

##############################################
# Crab bioenergetics
############################################
#Zoeae 1
#Yamamoto et al 2014
#paul et al - zoeae max consumption is 7 copepods per day
aIMperiodZ1 <- 229.5
bIMperiodZ1 <- 0.63
a <-getIMPeriodLarva(5,aIMperiodZ1,bIMperiodZ1)
#Respiration in terms of dry body weight for zoeae 1 & 2
R <- 1.198*w^(0.72)
#Zoeae 2 
#Yamamoto et al 2014
aIMperiodZ2 <- 530.5
bIMperiodZ2 <- 0
b <- getIMPeriodLarva(5,aIMperiodZ2,bIMperiodZ2)

#Megalope
#Consumption described by grams wet wt of Artemia nauplii consumed
#Yamamoto et al 2015
Consumption<- 1e4*CW-28077

#Yamamoto et al 2014
aIMperiod <- 417.3
bIMperiod <- -2.24 
c<- getIMPeriodLarva(5,aIMperiod,bIMperiod)



require(purrr)

#Proportion of max consumption
siz= c(3.19, 5.09, 7.34,10.03,14.31, 19.55, 25.98, 33.86, 43.53, 55.38,
       69.91, 87.72, 109.56)
daysPerStage <- unlist(getDays(5, siz))
ages <- 1:11 %>% map_dbl(~sum(daysPerStage[1:.x]))
mortFun <- getMort(temp = 1.83,CW = siz, plotTrue = F, 
                   z1Mort = .0161, AdultMort = .00063, longevity=6570)
curve(mortFun(x), from=0, to=6570, ylab="Daily mortality rate", xlab="Age (days)", las=1)
fullages <- c(0,ages,6570)
MidAges <- 2:13 %>% map_dbl(~mean(c(fullages[.x-1],fullages[.x])))

deltaw5 <- 2:length(predwt) %>% map_dbl(~predwt[.x]-predwt[.x-1])
iMCalc <- getIMPeriod(5,siz)
imActual <- c(iMCalc[1:8], rep(365,4))
g5 <- log(deltaw5)/imActual

mortFun(MidAges)
avgMaturity <- fullages[12]+365*2
avgTimesSpawning <- (6570-avgMaturity)/365


#Global parameters
FA=.89
UA=2.95
aR=exp(3.966)
bR=.722
cR=13.41
dR=5.36
QR=2.2
ACT=1.46
timestep=1
TmC = 14
dC = 9
QC = 5.5
aSDA = .019
bSDA = .17

#immature male parameters
ExA <- .001
ExB <- 2.243
aLWs <- 1.31
bLWs <- .946
aLW <- 2.8558E-4
bLW <- 3.081777
drywetRatio <- .136
perGramCal <- 3900
pC = .9

CW<- c(9.7, 12.9,
       17, 22.5, 29.7, 39.9, 51.5)
w0 <- (aLW*CW^bLW)
cons <- c(7,13,22,40,76,143,269)
dat <- data.frame(weight=w0, consumption=cons)
consumption <- lm(log(cons)~log(w0), data=dat)
plot(exp(fitted(consumption))~CW, type="l", log="xy")
points(cons~CW)
curve(20*((aLW*x^bLW))^.6, from=9, to=52, add=T, col="blue")



SAmolt <- data.frame(initCW = c(21.23, 22.2, 23.48,29.9,30.3,30.7,
                                44.2,44.7,64.7,67.6,67.9,74.5,79.9,
                                89.8,89.9,89.9,93.8,20),
                     postCW = c(26.41,28.1,28.27,39.9,40.3,40.5,
                                58.7,57.3,82.7,86,85.3,93.9,97.8,
                                110,112.1,112.3,117.6,26.3))
SAmoltFem <- data.frame(read.csv("StockAssessmentLengths.csv"))
moltIncrementModF <- lm(FemPost~FemPre, data=SAmoltFem)
moltIncrementMod1<- lm(postCW~initCW,data=SAmolt)

#Get parameters from fit to extracted data
aC <-exp(coef(consumption)[1])
bC <-coef(consumption)[2]



aCW <- moltIncrementMod1$coefficients[1]
betaCW <-moltIncrementMod1$coefficients[2]
#Molt increment params only work for larger crab, use values from Comeau et al. 1998 for smaller crab
aCWs <- 1.29
bCWs <- 1.19
   
temps <- 0:13
#Single crab across temperatures
#Weight of a male crab with a CW of 9.7 mm (using .136 ratio dry:wet wt)
g <- rep(0,length(siz))
for(i in 1:length(siz)){
  CW <- siz[i]
  g[i]<- plotOneCrab(CW=siz[i], temps=meantemp, mature=F)$growth
}


#################################################################################
############## Run projection

#Values to tweak to calibrate
#Ratio of dry:wet - mean .136, sd .13

#Calories per gram of somatic tissue - mean 3900, sd .4
perGramCal <- 3900
#Proportion of max consumption
siz= c(3.19, 5.09, 7.34,10.03,14.31, 19.55, 25.98, 33.86, 43.53, 55.38,
       69.91, 87.72, 109.56)
pC = c(1.33,
       1.1,
       1.02,
       1.02,
       0.97,
       0.95,
       0.94,
       0.94,
       0.96,
       1,
       1.05,
       1.11,
       1.19)

pCmat <- data.frame(CW=siz, propCons=pC) 
out<- runForYear(CWinit=3.19, temperature=5, duration=2000, FA=FA, UA=UA,aR=aR, bR=bR, cR=cR, dR=dR, eR=QR, 
           ACT=ACT, aSDA=aSDA, bSDA=bSDA, ExA=ExA, ExB=ExB,
           dt=timestep, mature=F, 
           RespFunc="Q10", cons.mat=pCmat, aC=aC, bC=bC, TmC=TmC, dC=dC , Qc=QC, perGramCal=perGramCal,
           CWslope1=betaCW, CWint=aCW, CWslopes=bCWs, CWints = aCWs)
pdf("RespTempCurve.pdf")
plot(Q10(aR, bR, cR, dR, QR, seq(0,10,by=.2),1)~
       seq(0,10,by=.2), type="l", 
     xlab="Temperature (Celsius)", 
     ylab="Respiration (microliter/hr)", las=1, 
     cex.axis=1.8, cex.lab=1.8, ylim=c(0,55))
dev.off()

pdf("ConsTempCurve.pdf")
lines(Q10(aC, bC, TmC, dC, QC, seq(0,10,by=.2),1)~
       seq(0,10,by=.2), type="l", 
     xlab="Temperature (Celsius)", 
     ylab="Respiration (microliter/hr)", las=1, 
     cex.axis=1.8, cex.lab=1.8)
dev.off()
predwt1<-aLW*out[[3]]^bLW

maleTestMetrics<- cbind(seq(1,14),c(0,out[[2]]), out[[3]], predwt1)
write.csv(maleTestMetrics,"MaleCheckValues.csv")

out[[1]][out[[2]]]
predwt1[-1]
plot(out[[1]], xlab="Days", ylab="Weight (g)", las=1, log="y")
#Caloric energy of a crab

predwt2 <- -10.154+3.453*out[[3]]
points(x=c(1,out[[2]]), y=predwt1, col="red", pch=19)
points(x=c(1,out[[2]]), y=predwt2, col="blue", pch=19)
abline(v=648)
legend("topleft", c("Model predicted weight", "Fitted PF equation predicted weight per instar",
                    "Fitted Hebert equation predicted weight per instar"),
       pch=c(1,19, 19), col=c(1,2, 4))


temps<-seq(1,8)
dur<-750
days <- 1:dur
wAtTemp <- lapply(X=temps, FUN=runForYear, CWinit=7, duration=dur, FA=FA, UA=UA,aR=aR, bR=bR, cR=cR, dR=dR, eR=QR, 
                  ACT=ACT, aSDA=aSDA, bSDA=bSDA, ExA=ExA, ExB=ExB,
                  dt=timestep, mature=F, 
                  RespFunc="Q10", p=pC, aC=aC, bC=bC, TmC=TmC, dC=dC , Qc=QC, perGramCal=perGramCal,
                  CWslope1=betaCW, CWint=aCW, CWslopes=bCWs, CWints = aCWs)
linecols <- rainbow(8)

plot(0, type="l", xlim=c(0,dur), ylim=c(0,50), ylab="weight (g)", xlab="days", las=1, cex.axis=1.5)
for(i in 1:8){
  lines(wAtTemp[[i]][[1]]~days, col=linecols[i])
  predwt1<-aLW*wAtTemp[[i]][[3]]^bLW/drywetRatio
  points(x=c(1,wAtTemp[[i]][[2]]), y=predwt1, col=linecols[i], pch=19)
}
legend("topleft", legend=temps, lty=1, col=linecols, cex=1.5)
#dev.off()

#Get biomasses:
init <- c(out[[1]][1],out[[1]][out[[2]]][1:(length(out[[2]])-1)])
final <-out[[1]][out[[2]]-1]
final/init


##########################################
## immature females!

#Fixed values:
ExA <- .000345
ExB <- 2.421
aLW <- .000065
bLW <- 2.869
aEggP=0.938
bEggP=-57.513
aEggM=1.117
bEggM=-52.761

#Calories per gram: mean 3700, sd 700
perGramCal <- 3700
#dry to wet ratio: mean .15 sd 1.4
drywetRatio <- .15
#proportion of max consumption
pC <- .91


#of weight-specific consumption Females
CW<- c(18.1, 23.9, 31.6, 41.7,53.6)
w0 <- aLW*CW^bLW
cons <- c(28,52,93,161,296)
dat <- data.frame(weight=w0, consumption=cons)
consumptionF<- lm(log(cons)~log(w0), data=dat)

aC <- exp(consumptionF$coefficients[1])
bC <- consumptionF$coefficients[2]

aCW <- moltIncrementModF$coefficients[2]
betaCW <-moltIncrementModF$coefficients[1]


plotOneCrab(CW=60, temps=-2:13, mature=T)

span <- 18*365
out<- runForYear(CWinit=3.19,temperature=5, duration=span, FA=FA, UA=UA,aR=aR, bR=bR, cR=cR, dR=dR, eR=QR, 
                 cons.mat = pCmat,
                 ACT=ACT, aSDA=aSDA, bSDA=bSDA, ExA=ExA, ExB=ExB,
                 dt=timestep, matureInit=F, 
                 RespFunc="Q10", p=pC, aC=aC, bC=bC, TmC=TmC, dC=dC , Qc=QC, perGramCal=perGramCal,
                 CWslope1=aCW, CWint=betaCW, CWints = aCWs, CWslopes=bCWs, aEggP=938, bEggP=-57513, aEggM=1117, bEggM=-52761, sex="F", primipInit=F)

#pdf("OrginialQ10Females.pdf")
plot(out[[1]], ylim=c(0,550), xlab="Days", ylab="Weight(g)", las=1, xlim=c(0,span))
#Caloric energy of a crab
predwt1<-aLW*out[[3]]^bLW/drywetRatio
points(x=c(1,out[[2]]), y=predwt, col="red", pch=19)
legend("topleft", c("Model predicted growth", "Predicted weight at each instar"), col=c(1,2),
       pch=c(1,19))
abline(v=out[[5]])
#dev.off()

plot(NA, xlim=c(0,1610), ylim=c(3,140), xlab="Days", ylab="Carapace width (mm)")
lines(maleRun[[3]]~c(0,maleRun[[2]]))
lines(femRun[[3]]~c(0,femRun[[2]]), lty=2)
points(maleRun[[3]]~c(0,maleRun[[2]]), col=colorMat[1], pch=111)
points(femRun[[3]]~c(0,femRun[[2]]), col=colorMat[11], pch=111)
legend("bottomright", c("Male", "Female"), pch=111, col=colorMat[c(1,11)])

temps<-seq(1,8)
dur<-900
days <- 1:dur
QC <- 5.5
wAtTemp <- runForYear(meantemp, CWinit=3, duration=dur, FA=FA, UA=UA,aR=aR, bR=bR, cR=cR, dR=dR, eR=QR, 
                  ACT=ACT, aSDA=aSDA, bSDA=bSDA, ExA=ExA, ExB=ExB,
                  dt=1, matureInit=F, 
                  RespFunc="Q10", p=pC, aC=aC, bC=bC, TmC=TmC, dC=dC , Qc=QC, perGramCal=perGramCal,
                  CWslope=aCW, CWint=betaCW, sex="F")
linecols <- rainbow(8)



#pdf("HigherQ10Female.pdf")

plot(0, type="l", xlim=c(0,dur), ylim=c(0,250), ylab="weight (g)", xlab="days", las=1, cex.axis=1.5)
for(i in 1:8){
  lines(wAtTemp[[i]][[1]]~days, col=linecols[i])
  predwt1<-aLW*wAtTemp[[i]][[3]]^bLW/drywetRatio
  points(x=c(1,wAtTemp[[i]][[2]]), y=predwt1, col=linecols[i], pch=19)
}
legend("topleft", legend=temps, lty=1, col=linecols, cex=1.5)
#dev.off()


#Weight of crab right before molt
postmolt<- out[[1]][out[[2]]-1]
#Starting weight just after molt
premolt<- c(out[[1]][1],out[[1]][out[[2]][1:(length(out[[2]])-1)]])
BiomassRatio <- postmolt/premolt
fit <- lm(BiomassRatio~out[[3]][1:6])
fitln <- lm(log(BiomassRatio)~log(out[[3]][1:6]))
plot(BiomassRatio~out[[3]][1:6], log="xy", xlab="CW (mm)", ylab="B1/B0")
lines(exp(fitted(fitln))~out[[3]][1:6])
coef(fitln)

CW <- out[[3]][1]
growthFun<- function(t, B){

  return(R)
} 


B<- rep(0,750)

CW <- out[[3]][1]
B[1] <- aLW*CW^bLW/drywetRatio
R <- coef(fitln)[1]+coef(fitln)[2]*log(CW)
D <- tMolt <- getIMPeriod(temp, CW)
k<-1
for(j in 2:750){
  B[j] <- exp(R/D)*B[j-1]
  if(j>=tMolt){
    D <- getIMPeriod(CW=out[[3]][k+1], temp=5)
    print(j)
    tMolt <- tMolt+D
    B[j] <- B[j]-(getExCosts(ExA, ExB, out[[3]][k], mature=F)/perGramCal)/drywetRatio
    k<-k+1
  }
}
plot(B)
points(out[[1]],col="red", pch=19)
points(x=c(1,out[[2]]), y=predwt,  col="blue")
