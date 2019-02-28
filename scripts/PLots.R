pdf("TempCosts50mm.pdf")
plotOneCrab(CW=50, temps=0:13, mature=F)
dev.off()

acrossSize <- plotOneCrab(CW=10:64, temps=0:13, mature=F)
pdf("CostsIntake7deg.pdf")
plot(acrossSize$costs[,8]~seq(10,64), xlab="CW (mm)", 
     ylab="Calories/day")
lines(acrossSize$intake[,8]~seq(10,64))
legend("bottomright",c("Costs", "Consumption"), 
       pch=c(1,NA), lty=c(NA,1), cex=.8)
dev.off()



#Check intermolt period
SteMarie <- c(6,12,17,22,28,39.6,51.6,68.4,80.4)*(365/12)

plot(SteMarie, pch=19, col="red", ylab="Intermolt period (days)", las=1, xlab="Molt #")
points(out[[2]], pch=19)
points(out[[2]], pch=19, col="blue")
legend("topleft", c("Sainte-Marie observations", "Predicted at mean bottom temperature for 50m (0.75)", 
                    "Predicted at mean bottom temperature for 100m (-0.1)"),
       pch=19, col=c(2,1,4), cex=2)



CW <- 11.6
temps<-seq(1,8)
dur<-750
days <- 1:dur
wAtTemp <- lapply(X=temps, runForYear, CWinit=CW, duration=dur, FA=FA, UA=UA,aR=aR, bR=bR, cR=cR, dR=dR, eR=QR, 
                  ACT=ACT, aSDA=aSDA, bSDA=bSDA, ExA=ExA, ExB=ExB,
                  dt=timestep, mature=F, 
                  RespFunc="Q10", p=pC, aC=aC, bC=bC, TmC=TmC, dC=dC , Qc=QC, perGramCal=perGramCal,
                  CWslope1=betaCW, CWint=aCW)
linecols <- rainbow(12)

#pdf("GrowthTempMaleJuveniles.pdf")
plot(0, type="l", xlim=c(0,dur), ylim=c(0,200), ylab="weight (g)", xlab="days", las=1, cex.axis=1.5)
for(i in 1:8){
  lines(wAtTemp[[i]][[1]]~days, col=linecols[i])
  predwt1<-aLW*wAtTemp[[i]][[3]]^bLW/drywetRatio
  points(x=c(1,wAtTemp[[i]][[2]]), y=predwt1, col=linecols[i], pch=19)
}
legend("topleft", legend=temps, lty=1, col=linecols, cex=1.5, pch=19)
#dev.off()