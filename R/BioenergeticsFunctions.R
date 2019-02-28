#'@description Function to get the numbers of days to molt under different temperatures and chela widths
#'@param temp ambient temperature
#'@param CWs vector of Chela widths
getDays <- function(temp, CWs){
  
  #Paramters for zoea1 from Yamamoto et al
  aIMperiodZ1 <- 229.5
bIMperiodZ1 <- 0.63
a <-getIMPeriodLarva(temp,aIMperiodZ1,bIMperiodZ1)

#Zoeae 2 
#Yamamoto et al 2014
aIMperiodZ2 <- 530.5
bIMperiodZ2 <- 0
b <- getIMPeriodLarva(temp,aIMperiodZ2,bIMperiodZ2)

#Yamamoto et al 2014 - megalopae
aIMperiod <- 417.3
bIMperiod <- -2.24 
c<- getIMPeriodLarva(temp,aIMperiod,bIMperiod)

#Get intermold period of immature crab
d <- lapply(X=CWs, FUN=getIMPeriod,temp=temp)
return(c(a,b,c,d))
}


#' @description function that extrapolates between adult and zoea1 mortality to get mortality rate at each life stage, assuming exponential decay of mortality rate with time
#' @param temperature - ambient mean temperature
#' @param CWs - vector, chela widths to look at
#' @param plotTrue - boolean, whether or not to plot mortality curve
#' #' @param z1Mort - scalar double, zoea1 mortality rate
#' @param AdultMort - scalar double, adult mortality rate
#' @param longevity - scalar int, how long adults live (days)
getMort <- function(temp, CWs, plotTrue, z1Mort, AdultMort, longevity){
  require(purrr)
  intervals <-unlist(getDays(temp, CWs))
  
  if(plotTrue){
    plot(intervals, col=c(1,2,3,rep(4,length(unlist(d)))), ylab="Days in stage")
  }
  
  #Map ages in each stage to total age (in days) at the end of each stage
  ages <- 1:length(intervals) %>% map_dbl(~sum(intervals[1:.x]))
   
  #Assume exponential decay of mortality

    #Calculate slope
   slope <- -(log(z1Mort)-log(AdultMort))/(mean(c(0,ages[1]))-mean(c(ages[8],longevity)))
    alpha <- -(log(z1Mort)-slope*mean(c(0,ages[1])))

  
  return(function(x) exp(-(alpha+slope*x)))
}



powerLaw <- function(a, W, b){
  return(a*W^b) 
}


#'Function to calculate growth rate based on bioenergetics
#'@param aC coefficient of weight power law relationship with max consumption
#'@param bC Exponent of weight power law relationship with max consumption
#'@param TmC max temperature at which consumption occurs
#'@param dC difference between TmC and temperature at which consumption is maximized
#'@param c1c consumption coefficient
#'@param FA proportion of consumed energy that is absorbed (less egestion)
#'@param UA amount of excretion per g of crab
#'@param aR coefficient of resting metabolism power law with weight
#'@param bR exponent of resting metabolism power law with weight
#'@param ACT activity multiplier for energetic budget
#'@param CW initial carapace width
#'@param dt growth interval
#'@param perGramCal calories per gram of wet weight
#'@param eggE energy diverted towards egg production by females
#'@param cR temperature of max respiration
#'@param dR difference between temperature at max respiration at temperature above which respiration stops
#'@param eR Q10 coefficient for respiration
#'@param aSDA - coefficient for weight-specific costs of specific dynamic action
#'@param bSDA - coefficient for weight-specific exponent of specific dynamic action
#'@param temp - ambient temperature
#'@param w0 - minimum weight at first instar
#'@param mature boolean, if the individual is mature
#'@param respFunc string, which respiration curve to use; one of Gutermuth_Armstrong, NoInteraction, Interaction, and Q10
#'@param p double - percentage of maximum consumption that is realized; this is a tuning parameter
#'@return growth over time dt
# Currently males only - needs extra for females
bioenergeticGrowthRate <- function(FA, UA, aR, bR, cR, dR, eR, ACT, aSDA, bSDA,
                                   dt, temp, w0, mature, RespFunc, 
                                   p, aC, bC, TmC, dC, Qc, perGramCal, eggE){

  # #Calculate consumption
   maxC <- aC*w0^bC
   tempC <- exp(calcF(temp, Tm=TmC, d=dC, Q=Qc))
   C <- p*maxC*tempC
   
  # #For now, use temperature-dependent function from Foyle et al 1989 that doesn't use weight
  # #in kcal/(kg)
  #C <- (2.2*w0*exp(-((temp-5.2)^2)/30.7))
  
  #Resting metabolic costs as a function of weight and temperature
  #This is returned in units microliter/ind/hr so
  #Convert to calories per day by multiplying .00463 cal per microliter o2 * 24 hours/day
  Rest <- switch(RespFunc,
                 "Gutermuth_Armstrong"=gutermuth_armstrong(aR,bR,cR,dR,eR,  Temp=temp, weight=w0),
                 "NoInteraction"=no_interaction(aR, bR, cR, Temp=temp, weight=w0),
                 "Interaction"=interaction(aR, bR, cR, dR, Temp=temp, weight=w0),
                 "Q10"=Q10(aR, bR, cR, dR, eR,  temp=temp, weight=w0))*.00463*24
  Abs <- C*FA
  #Total respiration including activity multiplier and specific dynamic action
  #units here are in kcal/hr
  SDA <- 24*aSDA*exp(bSDA*temp)/1000
  Resp <- (Rest*ACT)+SDA*(C-Abs)

  
  #Energy lost to excretion
  # UA*w0 is in units mg ammonia per day, 5.94 cal/mg of ammonia
  A <- (UA*w0*.00594)
  moltEnergy <- exp(.9876)*CW^(-.9281)
  #Growth rate is what's left over from absorbed energy after excretion, exuviae, and respiration 
  #in Calories, so daily growth in g is gr 
  #Absorbed energy
  g <- Abs-Resp-A-moltEnergy-eggE
  gramGrowth <- (g/perGramCal)/drywetRatio
  return(list(Rest=Rest, Resp=Resp, A=A, intake=C, growth=g, gramGrowth=gramGrowth, eggs=eggE))
}

#Function to get intermolt period for immature crab
getIMPeriod <- function(temp, CW){
  K <- 1573*(1-exp(-.07*CW))
  alpha <- -4.71*(1-exp(-.14*CW))
  D <- K/(temp-alpha)
  return(D)
}

#Larval intermolt period ()
getIMPeriodLarva <- function(temp, a, b){
  D <- a/(temp-b)
  return(D)
}

#'@param Tm - max temp for respiration
#'@param T0 - temp of max resp
#'@param a - respiration coefficient
calcF <- function(Temp,Tm, d, Q){
  v <- (Tm-Temp)/(d)
  w <- log(Q)*(d)
  y <- log(Q)*(d+2)
  x <- ((w*(1+sqrt(1+40/y)))^2)/400
  f <- x*log(v) + x*(1-v)
  return(f)
}

#'@param aR coefficient on temperature squared term
#'@param bR coefficient on temp:log(w) interaction
#'@param cR coefficient on temp term
#'@param dR coefficient on log weight
#'@param eR intercept
gutermuth_armstrong <- function(aR, bR, cR, dR, eR, Temp, weight){
  return(exp(aR*(Temp^2)+bR*Temp*log(weight)+cR*Temp+dR*log(weight)+eR))
}

#'@param aR coefficient on temperature effect
#'@param bR coefficient on log weight effect
#'@param cR intercept
no_interaction <- function(aR, bR, cR, Temp, weight){
  return(exp(aR*(Temp)+bR*log(weight)+cR))
}

#'@param aR coefficient on temperature respiration effect
#'@param bR coefficient on log weight respiration effect
#'@param cR coefficient on temperature:respiration interaction
#'@param dR intercept
interaction <- function(aR, bR, cR, dR, Temp, weight){
  return(exp(aR*(Temp)+bR*log(weight)+cR*Temp*log(weight)+dR))
}

#' @param aR coefficient of weight specific max respiration
#' @param bR exponent of weight specific max respiration
#' @param cR max temperature for respiration
#' @param dR temperature at which respiration is maximized
#' @param eR respiration coefficient
Q10 <- function(aR, bR, cR, dR, eR, temp, weight){
  maxR = log(aR)+bR*log(weight)
  
  r = maxR+calcF(temp,cR,dR,eR)
  return(exp(r))
}

#Function to fit the Q10 relationship using log-likelihood - need to fix quite a few parameters to get it to converge
fitQ10 <- function(sigma, aR, bR, cR, dR, eR, temp, weight, resp){
  pred.resp <- function(aR, bR, cR, dR, eR, temp, weight){
    #cR_trans <- get_logistic_transform(cR, 8, 20)
    #dR_trans <- get_logistic_transform(dR, 0.5, 7)
    maxR = log(aR)+bR*log(weight)  
    r = maxR+calcF(temp,cR,dR,eR)
    return(r)
  }

  predicted <- pred.resp(aR, bR, cR, dR, eR, crab_GraphClick$Temp, crab_GraphClick$Weight)
  logLik <- sum(-(predicted-log(crab_GraphClick$Respiration))^2)
            
  return(-logLik)
}

inv_logistic <- function(par, lo, hi){
  par <- -log(hi-par)+log(par-lo)
  return(par)
}

get_logistic_transform <- function(par, lo, hi){
  par <- lo + (hi - lo)/(1+exp(-par))
  return(par)
}

#' @param temperature - ambient temperature
#' @param duration - number of days to grow the crab
#' @param CWinit - starting carapace width in mm
#' @param mature - is crab mature?
#' @param CWslope - slope of post:pre CW line
#' @param CWint - intercept of post:pre CW line
#' @param ExA - coefficient of exuvia weight relationship
#' @param ExB - exponent OR intercept (if mature) of exuvia:weight relationship
#' @return list of crab weight, days on which molts occurred, and post molt CW on molting days
runForYear <- function(temperature, duration=365, CWinit, CWslope1, matureInit = F, primipInit=F,
                       CWint, CWints, CWslopes, ExA, ExB, sex="m", aEggP, bEggP, aEggM, bEggM, 
                       cons.mat, ...){
  #Create objects
  weight<-rep(0,duration)
  moltTime<-NULL
  eggWt <- 0
  maturityT<-NULL
  
  #Initialize starting values
  CWmat <- CWinit
  CW <- CWinit
  mature <- matureInit
  primip <- primipInit
  weight[1] <- (aLW*CW^bLW)
  tMolt <- getIMPeriod(temperature,CW)
  eggEnergy <- 0

  #Run for a year
  for(t in 2:duration){
 
    #if crab is mature, no exuvial costs
    if(!mature){
      moltInc <- getIMPeriod(temperature,CW)
      moltE <- getExCosts(ExA, ExB, CW, mature)
      moltEnergy <- moltE/moltInc
    } else{
      moltEnergy<-0
      
      #If mature female, calculate egg energy
      if(sex=="F"){
        if(primip){
          eggEnergy <- aEggP*CW+bEggP

        } else{
          eggEnergy <- aEggM*CW+bEggM
        }
      }
    }
    
    pred.function <- approxfun(x=cons.mat$CW, y=cons.mat$propCons)
    pC <- pred.function(CW)
    if(length(pC)==0){
      return(paste(CW, " carapace width is not valid."))
      }
    #Calculate daily growth
    bioE <- bioenergeticGrowthRate(FA=FA, UA=UA,aR=aR, 
                                   bR=bR, cR=cR, dR=dR, eR=QR, 
                                   ACT=ACT, aSDA=aSDA, bSDA=bSDA, 
                                   dt=timestep, temp=temperature, w0=weight[t-1], mature=mature,
                                   moltEnergy=moltEnergy, RespFunc="Q10",
                                   p=pC, aC=aC, bC=bC, TmC=TmC, dC=dC, Qc=QC, perGramCal=perGramCal,
                                   eggE = eggEnergy/365)
    if(!mature){
      if(t>tMolt){

      #Length-based probability of terminal molt
      pMaturing <- maturityProb(sex, CW, temperature)
      if(pMaturing>0){
        die.roll <- rbinom(1,1,pMaturing)
        if(die.roll==1){
          maturityT <-t
          moltE <- getExCosts(ExA, ExB, CW, T)
          moltEnergy <- moltE/moltInc
          if(sex=="F"){
            #if this is the molt to maturity, the female is primiparous
            primip <- T
          }
        }
      }
      
      
        #Append molt timing to vector, calculate new CW, get new molt timing, and subtract ex weight from crab weight
        moltTime<-c(moltTime, t)
        if(CW>9){
        CW <- CWint+CW*CWslope1
        } else{
          CW <- CWints+CW*CWslopes
        }

        moltInc <- getIMPeriod(temperature,CW)
        moltE <- getExCosts(ExA, ExB, CW, T)/moltInc
        CWmat <- c(CWmat, CW)
        tMolt <- t+moltInc
        #weight[t] <- (aLW*CW^bLW)/drywetRatio
        weight[t] <- bioE$gramGrowth +
          weight[t-1] - (moltE)/perGramCal
        
      } else{
        weight[t] <-  bioE$gramGrowth+weight[t-1]   
      }
      } else{
        if(primip){primip<-F}
        #If the crab is mature check if it is time to spawn
        if((t%%365==0)&mature&(sex=="F")){
          eggs <- (eggEnergy/perGramCal)/drywetRatio
          weight[t] <-  bioE$gramGrowth+weight[t-1]-eggs
          eggWt <- eggWt+eggs
        } else{
          weight[t] <-  bioE$gramGrowth+weight[t-1]
        }

    }

  }
  return(list(weight,moltTime, CWmat, eggWt, maturityT))
}

#' @description function to get exuvial costs of a crab of chela width CW
#' @param mature boolean
#' @param ExA coefficient of chela with
#' @param ExB exponent of chela width
getExCosts <- function(ExA, ExB, CW, mature){
  if(mature){
    Ex <- ExA*CW-ExB
  }else{
    Ex <- ExA*CW^ExB
  }
  #Convert to cal from kcal
  return(Ex*1000)
}


#'@description Function to plot the weight trajectory of a single crab at one chela width and multiple temps
plotOneCrab <-function(CW, temps, mature=F){
  w0 <- (aLW*CW^bLW)
  if(length(CW)==1){
  pred<-data.frame(costs=rep(0, length(temps)), intake=rep(0, length(temps)), 
                   growth=rep(0, length(temps)))
  } else{
    costs<-matrix(nrow=length(CW), ncol=length(temps))
    intake<-matrix(nrow=length(CW), ncol=length(temps))
  }
  for(i in 1:length(temps)){
    Ex<-getExCosts(ExA,ExB, CW, mature)

    growth <- bioenergeticGrowthRate(FA=FA, UA=UA,aR=aR, 
                                     bR=bR, cR=cR, dR=dR, eR=QR, 
                                     ACT=ACT, aSDA=aSDA, bSDA=bSDA, 
                                     dt=timestep, temp=temps[i], w0=w0, mature=F,
                                     RespFunc="Q10",
                                     p=pC, aC=aC, bC=bC, TmC=TmC, dC=dC, Qc=QC, perGramCal=perGramCal, 0)
    if(length(CW)==1){
    pred$costs[i] <- (growth$Resp+growth$A)
    pred$intake[i] <- growth$intake
    pred$growth[i] <- growth$growth
    pred$eggs[i] <- growth$eggs
    } else{
      costs[,i] <- growth$Resp+growth$A
      intake[,i] <- growth$intake
      eggs[,i] <- growth$eggs
    }
  }
  if(length(CW)==1){
  par(oma=c(2,2,2,2), cex=1.2)
  plot(pred$costs/w0~temps, pch=19, ylim=c(0,max(pred$intake/w0)),
       xlab="Temp (deg Celsius)", ylab="Calories/day", las=1,cex.axis=1.2)
  lines(pred$intake/w0~temps)
  legend("bottomleft",
         c("Costs","Consumption"), pch=c(19,NA), lty=c(NA,1),
         cex=1.2)
  return(pred)
  } else{
    plot(costs[1,]~temps, ylim=c(0,max(intake)))
    lines(intake[1,]~temps)
    for(i in 2:nrow(costs)){
      points(costs[i,]~temps, col=i)
      lines(intake[i,]~temps, col=i)
    }
    return(list(costs=costs,intake=intake))
  }
}

#'Calculates the probability of maturing during this molt based on CW (size)
#'Maturity probabilities taken from stock assessment
maturityProbSA <- function(sex, size){
  lengths <- seq(2.5, 200, 5)
  if(sex=="F"){
    probs <- c(rep(0,5), 0.00982844, 0.02801290, 0.07984330, 0.22516600, 0.55343100, 0.85447400, 0.85267100,
               0.89130300, rep(1.00000000, 27))
  } else{
    probs <- c(rep(0,5), 0.00606599, 0.01042500, 0.01791640, 0.03073110, 0.05169140, 0.08207390, 0.11872800,
               0.15459900, 0.17998200, 0.19536400, 0.20094100, 0.19821500, 0.20424600, 0.24323400,
               0.35519700, 0.62285300, 0.98887200, rep(1.00000000, 18))
  }
  #linearly interpolate to get probability at length
  pred.function <- approxfun(x=lengths, y=probs)
  return(pred.function(size))
}


#' @description Function to get the standard deviation of the distribution of maturity sizes according to temperature from the 1st and 99th quantile
#' quantiles taken from Burmeister & Sainte-Marie study
#' @param temp ambient temperature
#' @param a0 coefficient of post molt CW at maturity at 1st quantile
#' @param b0 exponent of post molt CW at maturity at 1st quantile
#' @param a1 coefficient of post molt CW at maturity at 99th quantile
#' @param b1 exponent of post molt CW at maturity at 99th quantile
standardDev <- function(temp, a0, a1, b0, b1){
  hi <- exp(a1*log(temp+2)+b1)
  onequant <- exp(a0*log(temp+2)+b0)
  return((hi-onequant)/(qnorm(.99)-qnorm(.01)))
}

#' @description Function to get the mean of the distribution of post-molt maturity sizes according to temperature from the 1st and 99th quantile
#' quantiles taken from Burmeister & Sainte-Marie study
#' @param temp ambient temperature
#' @param a0 coefficient of post molt CW at maturity at 1st quantile
#' @param b0 exponent of post molt CW at maturity at 1st quantile
#' @param a1 coefficient of post molt CW at maturity at 99th quantile
#' @param b1 exponent of post molt CW at maturity at 99th quantile
meanCW <- function(temp, a0, a1, b0, b1){
  hi <- exp(a1*log(temp+2)+b1)
  onequant <- exp(a0*log(temp+2)+b0)
  return((onequant*(qnorm(.99))-hi*(qnorm(.01)))/(qnorm(.99)-qnorm(.01)))
}


#'@description Function to use the normal distributions defined by functions above to calculate maturity probabilities
maturityProb <- function(sex, size, temp){
  lengths <- seq(2.5, 200, 5)
  
  #Females
  if(sex=="F"){
    #Get parameters of normal distribution describing post-molt cw
    sdVal <- standardDev(temp, .16,.12,3.81,4.30)
    meanVal <- meanCW(temp, .16,.12,3.81,4.30)
    
    #Parameters to translate the post-molt CW to pre-molt CW
    intercept <- -0.968
    slope <- 1.345
  } else{
    sdVal <- standardDev(temp, .20,.06,4,4.85)
    meanVal <- meanCW(temp, .20,.06,4,4.85)
    intercept <- 2.012
    slope <- 1.226
  }
  sdPre <- (sdVal/slope)
  meanPre <- (meanVal-intercept)/slope
  #linearly interpolate to get probability at length
  return(pnorm(size,meanPre,sdPre))
}

getQuant <- function(temp, slope, int){
  return(exp(slope*log(temp+2)+int))
}


#' @description function to calibrate p, the proportion of realized maximum consumption from the bioenergetics equations
findConsumptionMultiplier <- function(startCW, temp, initDur, cons.mat,...){
  #First, find how many iterations you have to run to get the size classes
  firstRun <- runForYear(CWinit=startCW, temperature=temp, duration=initDur, 
             cons.mat=cons.mat, ...)

  searchInc <- seq(1.3,.7, by=-.01)
  
  
  for(i in 1:(nrow(cons.mat))){
    #if(i==10){browser()}
    minDiff <- 1e500
    for(j in 1:length(searchInc)){
      cons.mat[i,2] <- searchInc[j]
      eachRun <- runForYear(CWinit=cons.mat[i,1], temperature=temp, duration=(firstRun[[2]][i])+1, 
                 cons.mat=cons.mat, ...)
      predWts <- (aLW*eachRun[[3]]^bLW)
      diffQ <- abs(eachRun[[1]][eachRun[[2]]]-predWts[2])
      if(!is.na(diffQ)){
      if(diffQ<minDiff){
        minDiff<-diffQ
        index <- j
      }
      }
      
    }
    cons.mat[i,2] <- searchInc[index]
    print(cons.mat)
  }
  return(cons.mat)
}


#'@description - Function to automate writing the bioenergetic parameters to the XML input to DisMELS
#'Not sure I ever finished this one!
getParsAndMakeXML <- function(LifeStage){
  #Global parameters
  FA=.89
  UA=2.95
  aR=exp(3.966)
  bR=.722
  cR=13.32
  dR=5.28
  QR=2.2
  ACT=1.46
  timestep=1
  TmC = 14
  dC = 9
  QC = 2.2
  aSDA = .019
  bSDA = .17
  
  
  if(LifeStage=="ImmatureMale"){
    #immature male parameters
    ExA <- .001
    ExB <- 2.243
    aLWs <- 1.31
    bLWs <- .946
    aLW <- .00005
    bLW <- 2.903
    drywetRatio <- .136
    perGramCal <- 3900
    pC = 1.03
    aCWs <- 1.29
    bCWs <- 1.19
    
    #Get parameters from fit to extracted data
    aC <-exp(coef(consumption)[1])
    bC <-coef(consumption)[2]
    
    aCW <- moltIncrementMod1$coefficients[1]
    betaCW <-moltIncrementMod1$coefficients[2]

      
    #Run once to get carapace widths
    carapaceWidth <- runForYear(CWinit=3.19, temperature=5, duration=2000, FA=FA, UA=UA,aR=aR, bR=bR, cR=cR, dR=dR, eR=QR, 
                           ACT=ACT, aSDA=aSDA, bSDA=bSDA, ExA=ExA, ExB=ExB,
                           dt=timestep, mature=F, 
                           RespFunc="Q10", cons.mat=pCmat, aC=aC, bC=bC, TmC=TmC, dC=dC , Qc=QC, perGramCal=perGramCal,
                           CWslope1=betaCW, CWint=aCW, CWslopes=bCWs, CWints = aCWs)
    siz <- carapaceWidth[[3]]
    pC = rep(1.1,length(siz))
    
    ConsumptionBYSizeM <- data.frame(CW=siz, propCons=pC) 
    pCmat<- findConsumptionMultiplier(startCW=siz[1], temp=5, initDur = 2000,
                                      cons.mat=ConsumptionBYSizeM, 
                                      FA=FA, UA=UA,aR=aR, bR=bR, cR=cR, dR=dR, eR=QR, 
                                      ACT=ACT, aSDA=aSDA, bSDA=bSDA, ExA=ExA, ExB=ExB,
                                      dt=timestep, mature=F, 
                                      RespFunc="Q10", aC=aC, bC=bC, TmC=TmC, dC=dC , Qc=QC, perGramCal=perGramCal,
                                      CWslope1=betaCW, CWint=aCW, CWslopes=bCWs, CWints = aCWs)
    write.csv(pCmat, "pCByCarapaceWidthM.csv")
    
  } else if(LifeStage=="ImmatureFemale"){
    
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
    
    aC <- exp(consumptionF$coefficients[1])
    bC <- consumptionF$coefficients[2]
    
    aCW <- moltIncrementModF$coefficients[2]
    betaCW <-moltIncrementModF$coefficients[1]
    
    #Run once to get carapace widths
    carapaceWidth <- runForYear(CWinit=3.19, temperature=5, duration=2000, FA=FA, UA=UA,aR=aR, bR=bR, cR=cR, dR=dR, eR=QR, 
                                ACT=ACT, aSDA=aSDA, bSDA=bSDA, ExA=ExA, ExB=ExB,
                                dt=timestep, mature=F, 
                                RespFunc="Q10", cons.mat=pCmat, aC=aC, bC=bC, TmC=TmC, dC=dC , Qc=QC, perGramCal=perGramCal,
                                CWslope1=betaCW, CWint=aCW, CWslopes=bCWs, CWints = aCWs, aEggP=938, bEggP=-57513, aEggM=1117, 
                                bEggM=-52761, sex="F", primipInit=F)
    siz <- carapaceWidth[[3]]
    pC = rep(1.1,length(siz))
    
    ConsumptionBYSizeF <- data.frame(CW=siz, propCons=pC) 
    pCmat<- findConsumptionMultiplier(startCW=siz[1], temp=5, initDur = 2000,
                                      cons.mat=ConsumptionBYSizeM, 
                                      FA=FA, UA=UA,aR=aR, bR=bR, cR=cR, dR=dR, eR=QR, 
                                      ACT=ACT, aSDA=aSDA, bSDA=bSDA, ExA=ExA, ExB=ExB,
                                      dt=timestep, mature=F, 
                                      RespFunc="Q10", aC=aC, bC=bC, TmC=TmC, dC=dC , Qc=QC, perGramCal=perGramCal,
                                      CWslope1=betaCW, CWint=aCW, CWslopes=bCWs, CWints = aCWs)
    write.csv(pCmat, "pCByCarapaceWidthF.csv")
    
    
  } else if(LifeStage=="AdolescentMale"){
    
  }
                            
  paramList <- list(aC = aC, bC =bC, cmT = TmC, 
                    coT = TmC-dC, c1c = QC, ACT = ACT,
                    aR = aR, bR = bR, rmT = cR, roT = cR-dR,
                    c1r = QR, FA = 1-FA, aSDA = aSDA, bSDA = bSDA,
                    UA = UA, sigRt=0, wRat = drywetRatio, calPerGram=perGramCal)
  descr <- c("linear coefficient of weight-dependent max consumption",
                "exponent coefficient of weight-dependent max consumption",
                "max temperature at which consumption occurs",
                "temperature at which consumption is maximized",
                "consumption coefficient",
                "respiration activity multiplier",
                "linear coefficient of weight-dependent respiration",
                "exponent coefficient of weight-dependent respiration",
                "max temperature at which respiration occurs",
                "temperature at which respiration is maximized",
                "respiration coefficient",
                "fraction of consumption lost to egestion",
                "coefficient of fraction of assimilated energy lost to SDA",
                "exponent of fraction of assimilated energy lost to SDA",
                "excretion fraction",
                "daily cost of exuviae",
                "std. dev. in linear growth rate",
                "dry to wet weight ratio of crab",
                "calories per gram of crab tissue")
  openString <- "<void method=\"setParameterValue\">\n<string>"
  bwNameValue <- "</string>\n<double>"
  bwvalueDescrip <- "</double>\n</void>\n<void method=\"setParameterDescription\">\n<string>"
  bwDescrip <- "</string>\n<string>"
  endStr <- " </string></void>\n"
  XMLString <- ""
  for(i in 1:length(paramList)){
   XMLString <- paste(XMLString, openString,
                      names(paramList)[i],
                      bwNameValue,
                      round(paramList[[i]],2),
                      bwvalueDescrip,
                      names(paramList)[i],
                      bwDescrip, 
                      descr[i],
                      endStr, sep="")
 }                          
  writeLines(XMLString, "BioEGrowthParsFem.xml") 
  
  paramList <- list(a=aCW, b=betaCW, mat=F)
  descr <- c("intercept of molt increment", "exponent of molt increment", "if crab is mature")
  openString <- "<void method=\"setParameterValue\">\n<string>"
  bwNameValue <- "</string>\n<double>"
  bwvalueDescrip <- "</double>\n</void>\n<void method=\"setParameterDescription\">\n<string>"
  bwDescrip <- "</string>\n<string>"
  endStr <- " </string></void>\n"
  XMLString <- ""
  for(i in 1:length(paramList)){
    XMLString <- paste(XMLString, openString,
                       names(paramList)[i],
                       bwNameValue,
                       round(paramList[[i]],2),
                       bwvalueDescrip,
                       names(paramList)[i],
                       bwDescrip, 
                       descr[i],
                       endStr, sep="")
  }                          
  writeLines(XMLString, "BioEMoltInc.xml") 
  
  
  return(pCmat)
}