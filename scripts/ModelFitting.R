#Values at 3,5,7 deg c respectively
a <- c(28.76, 42.47, 152.44)
b <- c(0.7178, 0.8771, 0.5266)
W <- exp(seq(log(1.3),log(700),length.out=20))
plot(powerLaw(a[1],W,b[1])~W, log="xy", xlim=c(1,700), ylim=c(20,15000))
points(powerLaw(a[2],W,b[2])~W)
points(powerLaw(a[3],W,b[3])~W)

#Make fake data
crab_PaulFuji<- data.frame(Weight=rep(log(W),3), Respiration=log(c(powerLaw(a[1],W,b[1]), 
                                                                   powerLaw(a[2],W,b[2]), powerLaw(a[3],W,b[3]))),
                           Temp = rep(c(3,5,7),each=20))


#Read in real data
crab_GraphClick <- read.csv("TannerCrabGrowthAllTemp.csv", header=T)
#Plot real data
plot(Respiration~Weight, col=Temp, log="xy", data=crab_GraphClick, pch=19, xlim=c(.001,600))

RespMixtureModel <- function(alphas, betas, delta, betal, sigma, s){
  R0 <- alphas+betas*log(Weight)
  R1 <- alphas + delta*(betas-betal)+betal*log(Weight)
  Phi1 <- 1/(1+exp(-s*(Weight-delta)))
  Phi0 <- 1-Phi1

  logLik <- sum(-log(sigma)/2+log(exp(-((log(Respiration)-R0)/(2*sigma))^2)*Phi0+
                  exp(-((log(Respiration)-R1)/(2*sigma))^2)*Phi1))  
  return(-logLik)
}

getR<-function(alphas, betas, betal, delta, s){
R0 <- alphas+betas*log(Weight)
R1 <- alphas + delta*(betas-betal)+betal*log(Weight)
Phi1 <- 1/(1+exp(-s*(Weight-delta)))
Phi0 <- 1-Phi1
return(R0*Phi0+R1*Phi1)
}

loglik2 <- lapply(seq(1,10,length.out=1000), RespMixtureModel,  betas=.721, betal=.719, delta=10, sigma=1, s=1)
plot(unlist(loglik2)~seq(1,10,length.out=1000))

attach(crab_GraphClick)
dn<-seq(20.7,20.85, length.out=20)
for(i in 1:20){
fitMix<- mle2(RespMixtureModel, start=list(betas=.68, betal=.67, delta=20.8),
              data=list(crab_GraphClick), fixed=list(sigma=1, alphas=log(aR), s=10),
              method="SANN")
print(summary(fitMix))
}
detach(crab_GraphClick)

RespPar<- coef(fitMix)
r<- getR(RespPar[1], RespPar[2], RespPar[4], RespPar[3], RespPar[6])
plot(r~Weight, log="xy", type="l")
points(log(Respiration)~crab_GraphClick$Weight, col="blue")
lines(r~Weight)

#linear model with no interaction term
require(dplyr)
weightModel <- lm(log(Respiration)~log(Weight), data=crab_GraphClick)
summary(weightModel)
#log(a)=3.966; b=1.722
lines(exp(fitted(weightModel))~crab_GraphClick$Weight)
#linear model with no interaction term
linModel <- lm(log(Respiration)~log(Weight)+Temp, data=crab_GraphClick)
summary(linModel)
lines(exp(fitted(linModel)[1:19])~crab_GraphClick$Weight[1:19], pch=19)
lines(exp(fitted(linModel)[20:37])~crab_GraphClick$Weight[20:37], pch=19)
lines(exp(fitted(linModel)[38:57])~crab_GraphClick$Weight[38:57], pch=19)

#Add interaction term
intLinModel <- lm(log(Respiration)~log(Weight)+Temp+log(Weight)*Temp, data=crab_GraphClick)
summary(intLinModel)
lines(exp(fitted(intLinModel)[1:19])~crab_GraphClick$Weight[1:19], pch=19, col=6)
lines(exp(fitted(intLinModel)[20:37])~crab_GraphClick$Weight[20:37], pch=19, col=6)
lines(exp(fitted(intLinModel)[38:57])~crab_GraphClick$Weight[38:57], pch=19, col=6)

#Gutermuth Armstrong 1989
gaModel <- lm(log(Respiration)~I(Temp^2)+Temp+log(Weight) + log(Weight)*Temp, data=crab_GraphClick)
lines(exp(fitted(gaModel)[1:19])~crab_GraphClick$Weight[1:19], pch=19, col=3)
lines(exp(fitted(gaModel)[20:37])~crab_GraphClick$Weight[20:37], pch=19, col=5)
lines(exp(fitted(gaModel)[38:57])~crab_GraphClick$Weight[38:57], pch=19, col=7)
legend("topleft", col=c(3,5,7), legend=as.character(c(3,5,7)), pch=19)
summary(gaModel)

#Q10 rule
require(bbmle)
attach(crab_GraphClick)
cR<-13.32
dR<-8.04
fitQ <- mle2(fitQ10, start=list(cR=cR, dR=dR),
             data=list(crab_GraphClick), fixed=list(aR=exp(3.966), bR=.722, eR=2.2))
coef(fitQ)[3:4]


#Check predictions and plots across models
test <- data.frame(Temp=seq(1,20), Weight=seq(1,700))
predictedResp <- exp(predict(gaModel, newdata=test))
predict2 <- exp(predict(intLinModel, newdata=test))
predict3 <- exp(predict(linModel, newdata=test))
plot(predictedResp~test$Temp, xlim=c(.5,20), ylim=c(1,100000), 
     log="xy",ylab="Respiration (microL/hr)", col=2)
points(Respiration~Temp ,data=crab_GraphClick,  pch=19)
points(predict2~test$Temp, col=4, pch=19)
points(predict3~test$Temp, col=6, pch=19)
predict4 <- vector("list")
for(i in 1:12){
  predict4[[i]]<- Q10(exp(3.966), .722, 13.35, 8., 2.2, temp=i, weight=seq(1,700))
  points(predict4[[i]]~rep(i,700), col=7)
}

malExuv <- data.frame(CW=c(9.7, 12.9, 17.0, 22.5, 29.7, 39.9, 51.5), 
                      wt = c(0.3138616,  0.7556447,  1.7688743,  4.1961671,  9.8726727, 24.5227537, 53.8438290),
                      ex = c(2, 3, 4, 5, 7, 11, 18),
                      rest = c(6, 12, 20, 36, 68, 127, 240))


plot(ex/rest~CW, data=malExuv, type="l")
nonlog <- lm(ex/rest~CW, data=malExuv)
log <- lm(ex/rest~log(CW), data=malExuv)
loglog <- lm(log(ex/rest)~log(CW), data=malExuv)
points(exp(fitted(loglog))~malExuv$CW, col=2)
points(fitted(log)~malExuv$CW)
bin <- glm(cbind(ex,rest-ex)~CW, family=binomial(link="cloglog"), data=malExuv)
plot(residuals(bin)~fitted(bin))
points(fitted(bin)~malExuv$CW, col=3)

#With lines on real data
plot(Respiration~Temp, data=crab_GraphClick, xlim=c(1,13), pch=19, las=1,
     xlab="Temperature", ylab="Oxygen consumption (microL per crab per hr)", col="gray",
     axes=F, log="xy")
axis(1)
axis(2,las=1, line=-1)
curve(gutermuth_armstrong(aR=-.159, 
                          bR=-0.047, cR=2.002, dR=.974, eR=-1.713, Temp=x, weight=738.87), 
      from=1, to=20, add=T, lty=3, lwd=2)
curve(gutermuth_armstrong(aR=-.159, 
                          bR=-0.047, cR=2.002, dR=.974, eR=-1.713, Temp=x, weight=1.38), 
      from=1, to=20, add=T, lty=3, lwd=2)
curve(Q10(exp(3.966), .722, 13, 6.9, 16.19, temp=x, weight=738.87), from=1, to=20, add=T, col="Blue", lty=2)
curve(Q10(exp(3.966), .722, 13, 6.9, 16.19, temp=x, weight=1.38), from=1, to=20, add=T, col="Blue", lty=2)
legend("topright",c("Gutermuth & Armstrong 1989", "Kitchell et al. 1977", "Paul & Fuji 1989 observations"), lty=c(3,2,NA), pch=c(NA, NA, 19),
       col=c("black","blue","gray"))

lines(Q103~seq(1,700))
lines(Q105~seq(1,700))
lines(Q107~seq(1,700))
AIC(linModel, intLinModel, gaModel, fitQ)



#Fit data to molt increment to come up with predicted CW~f(init CW)
moltIncrementMale<-data.frame(initCW =c(7.9,11.6, 16.4, 21.4), 
                              postCW = c(11.6, 16.4,21.4,29.6))

SAmolt <- data.frame(initCW = c(21.23, 22.2, 23.48,29.9,30.3,30.7,
                                44.2,44.7,64.7,67.6,67.9,74.5,79.9,
                                89.8,89.9,89.9,93.8,20),
                     postCW = c(26.41,28.1,28.27,39.9,40.3,40.5,
                                58.7,57.3,82.7,86,85.3,93.9,97.8,
                                110,112.1,112.3,117.6,26.3))
SAmoltFem <- data.frame(read.csv("StockAssessmentLengths.csv"))
moltIncrementModF <- lm(FemPost~FemPre, data=SAmoltFem)
moltIncrementMod1<- lm(postCW~initCW,data=SAmolt)



plot(postCW~initCW, data=SAmolt)
curve(-5.23+1.52*x,from=20, to=32.11, add=T)
curve(-5.23+32.11*(1.52-1.15)+1.15*x,from=32.11, to=95, add=T)

-5.23+1.52*(CW)

#Fit consumed calories predicted from Paul & Fuji to get slope and asymptote for immature males
#of weight-specific consumption
CW<- c(9.7, 12.9,
       17, 22.5, 29.7, 39.9, 51.5)
w0 <- (aLW*CW^bLW)/drywetRatio
cons <- c(7,13,22,40,76,143,269)
dat <- data.frame(weight=w0, consumption=cons)
consumption <- lm(log(cons)~log(w0), data=dat)
plot(exp(fitted(consumption))~CW, type="l", log="xy")
points(cons~CW)
curve(20*((aLW*x^bLW)/drywetRatio)^.6, from=9, to=52, add=T, col="blue")



#Fit consumed calories predicted from Paul & Fuji to get slope and asymptote for immature males
#of weight-specific consumption Females
CW<- c(18.1, 23.9, 31.6, 41.7,53.6)
aLW <- 5.284922E-4
bLW <- 2.895
w0 <- aLW*CW^bLW
cons <- c(28,52,93,161,296)
dat <- data.frame(weight=w0, consumption=cons)
consumptionF<- lm(log(cons)~log(w0), data=dat)





####################
# Mature female model fitting

paulFujiEggs <- data.frame(num=c(rep(169300,5),
                                 rep(195000,38),
                                 rep(250000,58),
                                 rep(281000,9)),
                           wt=c(rep(7.11,5),
                                rep(8.19,38),
                                rep(10.15,58),
                                rep(11.8,9)),
                           cal=c(rep(44.1,5),
                                 rep(50.7,38),
                                 rep(62.9,58),
                                 rep(73.1,9)))

meanwt <- mean(paulFujiEggs$wt/paulFujiEggs$num)
meancal <- mean(paulFujiEggs$cal/paulFujiEggs$num)
meancal

#P(length at maturity) based on stock assessment
lengths <- seq(2.5, 132.5, 5)
probs <- c(rep(0,5), 0.00982844, 0.02801290, 0.07984330, 0.22516600, 0.55343100, 0.85447400, 0.85267100,
           0.89130300, rep(1.00000000, 14))
plot(probs~lengths, lwd=3, type="l", xlab="CW(mm)", ylab="p(terminal molt)", las=1)
lines(probs~predict(moltIncrementModF, newdata=data.frame(FemPre=lengths)), lwd=3)

randLengths <-runif(10000,min(lengths),103)
postLengths<- predict(moltIncrementModF, newdata=data.frame(FemPre=randLengths))
hist(randLengths)
mature <- unlist(lapply(randLengths, function (x) return(rbinom(1,1,maturityProb("F", x)))))
CW<- NULL
for(i in 1:length(mature)){
  if(mature[i]){
    preLength <- c(preLength,randLengths[i])
    CW <- c(CW, predict(moltIncrementModF, newdata=data.frame(FemPre=randLengths[i])))
  }
}
CW<- as.numeric(CW)
p <- hist(CW, freq=F)
probByCw <- glm(mature~randLengths, family=binomial)
summary(probByCw)
backTransformed <- get_logistic_transform(predict(probByCw), hi=1,lo=0)

plot(backTransformed~randLengths, xlab="CW(mm)", ylab="p(terminal molt)", las=1)
points(backTransformed~postLengths, col="blue")


#Use quantiles to compute median and standard deviation of the std normal

sdftemp <- standardDev(-1:8, .16,.12,3.81,4.30)
meanftemp<- meanCW(-1:8, .16,.12,3.81,4.30)

sdPre <- (sdftemp/coefficients(moltIncrementModF)[2])
meanPre <- (meanftemp-coefficients(moltIncrementModF)[1])/coefficients(moltIncrementModF)[2]
i<-5
par(mfrow=c(3,3))
for(i in 1:9){
  dat_pred<-rnorm(10000,meanftemp[i],sdftemp[i])
  transformed <- (dat_pred-coefficients(moltIncrementModF)[1])/coefficients(moltIncrementModF)[2]
  hist(transformed)
  abline(v=meanPre[i])
}


plot(pnorm(dat_pred,meanPre[5],sdPre[5])~dat_pred)
df <- rep(seq(-1,8), each=100)
CW <- rep(seq(20,80, length.out=100),10)
pdf("FemaleMoltProb.pdf")
plot(probs~lengths, lwd=3, type="l", xlab="CW(mm)", ylab="p(terminal molt)", las=1)
lines(probs~predict(moltIncrementModF, newdata=data.frame(FemPre=lengths)), lwd=3, col="blue")
points(x=tempq, y=rep(.5,10), col=linecol, pch=2, lwd=3, cex=.5)
temp99<-(getQuant(-1:8, .12,4.30)+.968)/1.345
points(x=temp99, y=rep(.99,10), col=linecol,pch=3, lwd=3, cex=.5)
temp1<-(getQuant(-1:8, .16,3.81)+.968)/1.345
points(x=temp1, y=rep(.01,10), col=linecol,pch=4, lwd=3, cex=.5)
for(i in 1:10){
 lines(sort(pnorm(CW,meanPre[i],sdPre[i]))~sort(CW), pch=1, col=linecol[i], lwd=1, lty=2)
}


legend("bottomright",c("Temperature-independent premolt CW", "Temperature-independent postmolt CW",
                   "Predicted Burmeister 1% quantile pre CW at temperatures -1:8",
                   "Predicted Burmeister median pre CW at temperatures -1:8",
                   "Predicted Burmeister 99% quantile pre CW at temperatures -1:8",
                   "Predicted temperature-dependent CW at temperatures -1:8"),col=c("black", "blue", rep("red",4)),
       pch=c(NA,NA,2,3,4,NA), lwd=c(3,3,3,3,3,1), lty=c(1,1,rep(NA,3),2), cex=.6)
legend("topleft", as.character(-1:8), col=linecols, lty=1, title="Temperature (degrees C)", cex=.6)
dev.off()


###################################
# Males
sdmtemp <- standardDev(-1:8, .20,.06,4,4.85)
meanmtemp<- meanCW(-1:8, .20,.06,4,4.85)

sdPre <- (sdmtemp/coefficients(moltIncrementMod1)[2])
meanPre <- (meanmtemp-coefficients(moltIncrementMod1)[1])/coefficients(moltIncrementMod1)[2]
CW <- rep(seq(20,80, length.out=100),10)
pnorm(CW,sdPre,meanPre)

temp99<-(getQuant(-1:8, .06,4.85)-coefficients(moltIncrementMod1)[1])/coefficients(moltIncrementMod1)[2]
temp1<-(getQuant(-1:8, .20,4)-coefficients(moltIncrementMod1)[1])/coefficients(moltIncrementMod1)[2]
tempq<-(getQuant(-1:8, .16,4.46)-coefficients(moltIncrementMod1)[1])/coefficients(moltIncrementMod1)[2]
probs <- c(rep(0,5), 0.00606599, 0.01042500, 0.01791640, 0.03073110, 0.05169140, 0.08207390, 0.11872800,
           0.15459900, 0.17998200, 0.19536400, 0.20094100, 0.19821500, 0.20424600, 0.24323400,
           0.35519700, 0.62285300, 0.98887200, rep(1.00000000, 5))
pdf("MaleMoltProb.pdf")
plot(probs~lengths, lwd=3, type="l", xlab="CW(mm)", ylab="p(terminal molt)", las=1)
lines(probs~predict(moltIncrementMod1, newdata=data.frame(initCW=lengths)), lwd=3, col="blue")
points(x=tempq, y=rep(.5,10), col=linecol, pch=2, lwd=3, cex=.5)
points(x=temp99, y=rep(.99,10), col=linecol,pch=3, lwd=3, cex=.5)
points(x=temp1, y=rep(.01,10), col=linecol,pch=4, lwd=3, cex=.5)
CW <- seq(40,120)
for(i in 1:10){
  lines(sort(pnorm(CW,meanPre[i],sdPre[i]))~sort(CW), pch=1,  lwd=1, lty=2)
}
legend("topleft",c("Temperature-independent premolt CW", "Temperature-independent postmolt CW",
                   "Predicted Burmeister 1% quantile pre CW at temperatures -1:8",
                   "Predicted Burmeister median pre CW at temperatures -1:8",
                   "Predicted Burmeister 99% quantile pre CW at temperatures -1:8",
                   "Predicted temperature-dependent CW at temperatures -1:8"),col=c("black", "blue", rep("red",4)),
       pch=c(NA,NA,2,3,4,NA), lwd=c(3,3,3,3,3,1), lty=c(1,1,rep(NA,3),2), cex=.6)
legend("left", as.character(-1:8), col=linecols, lty=1, title="Temperature (degrees C)", cex=.6)
dev.off()

pBySize<-lapply(X=0:10, FUN=maturityProbSA, sex="M", size=1:200)
pBySizeF<-lapply(X=0:10, FUN=maturityProbSA, sex="F", size=1:200)
par(cex=1.5)
plot(NA, xlim=c(30,140), ylim=c(0,1), las=1, ylab="P(terminal molt)", xlab="Pre-molt carapace width (mm)")
colorMat <- brewer.pal(length(pBySize),"Spectral")
for(i in 1:length(pBySize)){

lines(pBySize[[i]]~seq(1,200), col=colorMat[length(pBySize)-i+1])
lines(pBySizeF[[i]]~seq(1,200), lty=2, col=colorMat[length(pBySize)-i+1])
}
legend("bottomright", title="Temperature", legend=as.character(seq(0,10)), col=rev(colorMat), lty=1)
legend(y=.12, x=120, title="Sex", legend=c("Female", "Male"), lty=c(2,1))


#Find consumption multiplier for males
siz= c(3.19, 5.09, 7.34,10.03,14.31, 19.55, 25.98, 33.86, 43.53, 55.38,
       69.91, 87.72, 109.56)
pC = c(1.29, 1.11,  1.03, 1.04, .99,  .95,  .95,  0.95,  .9, 0.9,
       0.675, 1, 1)

ConsumptionBYSizeM <- data.frame(CW=siz, propCons=pC) 
pCmat<- findConsumptionMultiplier(startCW=3.19, temp=5, initDur = 2000,
                                  cons.mat=ConsumptionBYSizeM, 
                FA=FA, UA=UA,aR=aR, bR=bR, cR=cR, dR=dR, eR=QR, 
                 ACT=ACT, aSDA=aSDA, bSDA=bSDA, ExA=ExA, ExB=ExB,
                 dt=timestep, mature=F, 
                 RespFunc="Q10", aC=aC, bC=bC, TmC=TmC, dC=dC , Qc=QC, perGramCal=perGramCal,
                 CWslope1=betaCW, CWint=aCW, CWslopes=bCWs, CWints = aCWs)
write.csv(pCmat, "pCByCarapaceWidthM.csv")
pred.pC <- approxfun(x=pCmat$CW, y=pCmat$propCons)
curve(pred.pC(x), from=3, to=100, ylim=c(.5,1.5))
points(propCons~CW, data=pCmat)

#FInd consumption multiplier for females
siz= c(out[[3]]
pC = rep(1,length(out[[3]]))

ConsumptionBYSizeF <- data.frame(CW=siz, propCons=pC) 
pCmat<- findConsumptionMultiplier(startCW=80.23, temp=5, initDur = 1500,
                                  cons.mat=ConsumptionBYSizeF[12,], 
                                  FA=FA, UA=UA,aR=aR, bR=bR, cR=cR, dR=dR, eR=QR, 
                                  ACT=ACT, aSDA=aSDA, bSDA=bSDA, ExA=ExA, ExB=ExB,
                                  dt=timestep, matureInit=F, sex="F",
                                  RespFunc="Q10", aC=aC, bC=bC, TmC=TmC, dC=dC , Qc=QC, perGramCal=perGramCal,
                                  CWslope1=aCW, CWint=betaCW, CWslopes=bCWs, CWints = aCWs,
                                  primipInit=F, aEggP=938, bEggP=-57513, aEggM=1117, bEggM=-52761)
write.csv(pCmat, "pCByCarapaceWidthF.csv")
pred.pC <- approxfun(x=pCmat$CW, y=pCmat$propCons)
curve(pred.pC(x), from=3, to=100, ylim=c(.5,1.5))
points(propCons~CW, data=pCmat)
