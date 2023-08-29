# This Script is meant to run the the PPE over data collected by Tension and Anderson (unknown year)
# The goal of this script is to run the PPE the data and attempt to estimate 
# Cognitive Stages of learning 
#--------------------------------------------------------------------------------------------------------------------------------
#Load packages
library(dplyr)
library(R2jags)
library(Hmisc)
#--------------------------------------------------------------------------------------------------------------------------------
#Load Functions Need to Run this code
{
  
  Parameter_Format=function(parameter, distribution){
    if(distribution=="beta" ){ x = betaABfromMeanSD(    mean(parameter), sd(parameter)) }
    if(distribution=="gamma"){ x = gammaShRaFromMeanSD( mean(parameter), sd(parameter)) }
    return(paste(as.character(x[1]), as.character(x[2]) ) )
  }
  
  betaABfromMeanSD = function( mean , sd ) {
    if ( mean <=0 | mean >= 1) stop("must have 0 < mean < 1")
    if ( sd <= 0 ) stop("sd must be > 0")
    kappa = mean*(1-mean)/sd^2 - 1
    if ( kappa <= 0 ) stop("invalid combination of mean and sd")
    a = mean * kappa
    b = ( 1.0 - mean ) * kappa
    return( list( a=a , b=b ) )
  }
  
  getmode <- function(v) {
    
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
  }
  
  gammaShRaFromMeanSD = function( mean , sd ) {
    if ( mean <=0 ) stop("mean must be > 0")
    if ( sd <=0 ) stop("sd must be > 0")
    shape = mean^2/sd^2
    rate = mean/sd^2
    return( list( shape=shape , rate=rate ) )
  }
  HDIofMCMC = function( sampleVec , credMass=0.95 ) {
    # Computes highest density interval from a sample of representative values,
    #   estimated as shortest credible interval.
    # Arguments:
    #   sampleVec
    #     is a vector of representative values from a probability distribution.
    #   credMass
    #     is a scalar between 0 and 1, indicating the mass within the credible
    #     interval that is to be estimated.
    # Value:
    #   HDIlim is a vector containing the limits of the HDI
    sortedPts = sort( sampleVec )
    ciIdxInc = ceiling( credMass * length( sortedPts ) )
    nCIs = length( sortedPts ) - ciIdxInc
    ciWidth = rep( 0 , nCIs )
    for ( i in 1:nCIs ) {
      ciWidth[ i ] = sortedPts[ i + ciIdxInc ] - sortedPts[ i ]
    }
    HDImin = sortedPts[ which.min( ciWidth ) ]
    HDImax = sortedPts[ which.min( ciWidth ) + ciIdxInc ]
    HDIlim = c( HDImin , HDImax )
    return( HDIlim )
  }
}
#--------------------------------------------------------------------------------------------------------------------------------
#Load Data from Caitlen's  study
setwd( paste(dirname(rstudioapi::getSourceEditorContext()$path), sep="") ) 
#Data = read.csv("/Users/michaelcollins/Documents/Predictive Performance Equation/Catilin Three Stage Model/Data/TenisonSpacingData_PPE_Variables_v2.csv")
Data = read.csv(paste(dirname(rstudioapi::getSourceEditorContext()$path), "/TenisonSpacingData_PPE_Stage_Specific_Variables_v2.csv", sep="") )
Data = Data[Data$Subject.ID %in% unique(Data$Subject.ID)[1:20],]
#--------------------------------------------------------------------------------------------------------------------------------
#Create an ID varaible with a combination of the Subject's ID and the Problem type that they were given
Data$ID = paste(Data$Subject.ID, Data$Problem.ID)
#------------------------------------------------------------------
Behavior_File_Name = "PPE_Cal1_Pred3_20_Subjects_Kappa_50"
Parameter_File_Name ="PPE_Cal1_Pred3_Parameters"
Cal_Day = 1
#------------------------------------------------------------------
for(i in 1 : length(unique(Data$ID)) ){ 
  
  i=1
  print(paste(i, " out of ", length(unique(Data$ID))))
  Subject = Data[Data$ID == unique(Data$ID)[i],]
  #Subject = Subject[Subject$Instances >= 5,]
#--------------------------------------------------------------------------------
#This sets all of the parameters for the simulations
#Num_Cal                    - Number of events the model will use for a likelihood.
#Num_Sub                    - Number of Subjects per a simulation
#Perf                       - Performance measure
#N                          - Task exposures
#T_var                      - Model Time
#CumAveInvLogTime_diff_lag1 - Decay variable
Num_Cal                    =  length(Subject$Instances[Subject$TrialTimeDay < 3])
    #
Num_Sub                    = 1
Num_Rounds                 = length(Subject$N)
Perf                       = 1 - ( Subject$Latency[Subject$TrialTimeDay < 3] / max(Data$Latency) )
  #1 - ( Subject$Latency / max(Data$Latency) )
  #
Perf                       = ifelse(Perf == 0, .001, Perf)
N                          = Subject$N
T_var                      = Subject$T
CumAveInvLogTime_diff_lag1 = Subject$CumAveInvLogTime_diff_lag1
#--------------------------------------------------------------------------------
data = c("N", "T_var", "CumAveInvLogTime_diff_lag1", "Num_Sub", "Perf", "Num_Rounds", "Num_Cal")
#--------------------------------------------------------------------------------
#Set Parameters 
parameters <- c("b", "m", "a", "tau", "k_par", "Prediction",
                "b_Prior", "m_Prior", "a_Prior", "tau_Prior", "k_par_Prior", "Prediction_Prior")
myinits <-	list( 
  list(b = 0.15, m = .15, a = .5, tau = .5, k_par = 900, Prediction=rep(1,Num_Rounds) ), 
  list(b = 0.15, m = .15, a = .5, tau = .5, k_par = 900, Prediction=rep(1,Num_Rounds) ),
  list(b = 0.15, m = .15, a = .5, tau = .5, k_par = 900, Prediction=rep(1,Num_Rounds) ) 
  ) # Initial group assignment
#----------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------
#Jags Model
setwd( paste(dirname(rstudioapi::getSourceEditorContext()$path), sep="") ) 
samples <- jags(data, inits=myinits, parameters,
                model.file ="PPE_Standard_Model_Jags.txt", n.chains=3, n.iter=2000, 
                n.burnin=1000, n.thin = 1, DIC = T )
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
#Obtain all of the informaiton relevent from the model
{
b          = samples$BUGSoutput$sims.list$b
m          = samples$BUGSoutput$sims.list$m
a          = samples$BUGSoutput$sims.list$a
tau        = samples$BUGSoutput$sims.list$tau
Prediction = samples$BUGSoutput$sims.list$Prediction
k_par      = samples$BUGSoutput$sims.list$k_par
#----------------------------------------------------
b_Prior          = samples$BUGSoutput$sims.list$b_Prior
m_Prior          = samples$BUGSoutput$sims.list$m_Prior
a_Prior          = samples$BUGSoutput$sims.list$a_Prior
tau_Prior        = samples$BUGSoutput$sims.list$tau_Prior
Prediction_Prior = samples$BUGSoutput$sims.list$Prediction_Prior
k_par_Prior      = samples$BUGSoutput$sims.list$k_par_Prior
#---------------------------------------------------------------
Prediction_Mean           = (1 - apply(Prediction, 2, mean))            * max(Data$Latency)
Prediction_HDI            = (1 - apply(Prediction, 2, HDIofMCMC))       * max(Data$Latency)
Prior_Prediction_Mean     = (1 - apply(Prediction_Prior, 2, mean))      * max(Data$Latency)
Prior_Prediction_HDI      = (1 - apply(Prediction_Prior, 2, HDIofMCMC)) * max(Data$Latency)

#----------------------------------------
Subject$Prediction = Prediction_Mean
Subject$Prior_Prediction = Prior_Prediction_Mean
Subject$Prior_High_HDI   = Prior_Prediction_HDI[1,]
Subject$Prior_Low_HDI    = Prior_Prediction_HDI[2,]
Subject$High_HDI   = Prediction_HDI[1,]
Subject$Low_HDI    = Prediction_HDI[2,]
Subject$Cal_Point  = Num_Cal
}
#---------------------------------------------------------------------------------------------------------------------------
Parameters=cbind(unique(Subject$Subject.ID), 
                 unique(Subject$Problem.ID),
                 unique(Subject$Spacing),
                 unique(Subject$Base),
                 Parameter_Format(b,           "beta"),
                 Parameter_Format(m,           "beta"), 
                 Parameter_Format(a,           "beta"), 
                 Parameter_Format(tau,         "beta"), 
                 Parameter_Format(b_Prior,     "beta"),  
                 Parameter_Format(m_Prior,     "beta"), 
                 Parameter_Format(a_Prior,     "beta"),  
                 Parameter_Format(tau_Prior,   "beta"), 
                 Parameter_Format(k_par,       "gamma"),
                 Parameter_Format(k_par_Prior, "gamma")
                 )

colnames(Parameters) = c("Subject.ID", "Problem.ID", "Spacing", "Base",
                         "b_Post",  "m_Post",  "a_Post",  "tau_Post", 
                         "b_Prior", "m_Prior", "a_Prior", "tau_Prior",
                         "k_par", "k_par_Prior")
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
setwd("/Users/michaelcollins/Documents/Predictive Performance Equation/Catilin Three Stage Model/Results")
if(i == 1) { write.table(Subject, Behavior_File_Name , col.names = T)    } else { write.table(Subject, Behavior_File_Name , col.names=F, append = TRUE) }
#if(i == 1) { write.table(Parameters, Parameter_File_Name, col.names = T) } else { write.table(Parameters, Parameter_File_Name, col.names=F, append = TRUE) }
}
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
par(mfrow=c(2,2)   )
plot( density(b),   col=1, xlim=c(0,.25)   )
lines(density(b_Prior), col=1, lty=3 )

plot( density(m),   col=2, xlim=c(0,.25)   )
lines(density(m_Prior), col=2, lty=3 )

plot( density(a),   col=3, xlim=c(0,1)     )
lines(density(a_Prior), col=3, lty=3 )

plot( density(tau), col=4, xlim=c(0,1)     )  
lines(density(tau_Prior), col=4, lty=3     )
#------------------------------------------
par(mfrow=c(1,1))
plot( density(k_par), col=1, xlim=c(100,1000)   )  
lines(density(k_par_Prior), col=1, lty=3 )
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
par( mfrow = c( 1, 2 ) )
{
plot( Subject$Latency[Subject$TrialTimeDay < 3]/1000,
     type="b", lwd=2,
     xlab="Within Session Instance", ylab="Response Time", 
     main="Day 1",
     ylim=c(0, (max(c(Subject$Latency, Prediction_Mean)) / 1000) ))
lines(Prediction_Mean[Subject$TrialTimeDay < 3]/1000, col=2, lwd=2, lty=3, type="b")
errbar( 1:length(Prediction_Mean[Subject$TrialTimeDay < 3]),    
       Prediction_Mean[Subject$TrialTimeDay < 3]/1000, 
        yplus  = Prediction_HDI[1,1:length(Prediction_Mean[Subject$TrialTimeDay < 3])]/1000,
        yminus = Prediction_HDI[2,1:length(Prediction_Mean[Subject$TrialTimeDay < 3])]/1000, 
        add=T, col=2, errbar.col = 2)
# lines(Prior_Prediction_Mean[Subject$TrialTimeDay < 3]/1000, col=4, lwd=2, lty=3, type="b")
# errbar( 1:length(Prior_Prediction_Mean[Subject$TrialTimeDay < 3]),    
#         Prior_Prediction_Mean[Subject$TrialTimeDay < 3]/1000, 
#         yplus  = Prior_Prediction_HDI[1,1:length(Prior_Prediction_Mean[Subject$TrialTimeDay < 3])]/1000,
#         yminus = Prior_Prediction_HDI[2,1:length(Prior_Prediction_Mean[Subject$TrialTimeDay < 3])]/1000, 
#         add=T, col=4, errbar.col = 4)
#----------------------------------------------------------------------
plot(Subject$Latency[Subject$TrialTimeDay > 3]/1000, type="b", lwd=2,
     xlab="Within Session Instance", ylab="Response Time", 
     main="Day 3",
     ylim=c(0, (max(c(Subject$Latency, Prediction_Mean))/1000 ) ))
#lines(Prediction_Mean[Subject$TrialTimeDay > 1]/1000, col=2, lwd=2, lty=3, type="b")
lines(Prediction_Mean[Subject$TrialTimeDay > 3]/1000, col=2, lwd=2, lty=3, type="b")
# errbar(1:length(Prediction_Mean[Subject$TrialTimeDay > 3]),    
#        Prediction_Mean[Subject$TrialTimeDay > 1]/1000, 
#        yplus  = Prediction_HDI[1, length(Prediction_Mean[Subject$TrialTimeDay < 3])+1:length(Prediction_Mean[Subject$TrialTimeDay > 3]) ]/1000,
#        yminus = Prediction_HDI[2, length(Prediction_Mean[Subject$TrialTimeDay < 3])+1:length(Prediction_Mean[Subject$TrialTimeDay > 3]) ]/1000, 
#        add=T, col=2, errbar.col = 2 )
errbar(1:length(Prediction_Mean[Subject$TrialTimeDay > 3]),    
       Prediction_Mean[Subject$TrialTimeDay > 3]/1000, 
       yplus  = Prediction_HDI[1, Subject$TrialTimeDay > 3 ]/1000,
       yminus = Prediction_HDI[2, Subject$TrialTimeDay > 3 ]/1000, 
       add=T, col=2, errbar.col = 2 )
# lines(Prior_Prediction_Mean[Subject$TrialTimeDay > 3]/1000, col=4, lwd=2, lty=3, type="b")
# errbar( 1:length(Prior_Prediction_Mean[Subject$TrialTimeDay > 3]),    
#         Prior_Prediction_Mean[Subject$TrialTimeDay > 3]/1000, 
#         yplus  = Prior_Prediction_HDI[1,1:length(Prior_Prediction_Mean[Subject$TrialTimeDay > 3])]/1000,
#         yminus = Prior_Prediction_HDI[2,1:length(Prior_Prediction_Mean[Subject$TrialTimeDay > 3])]/1000, 
#         add=T, col=4, errbar.col = 4)
# #----------------------------------------------------------------------
}
#============================================================================================
mean(ifelse(Subject$Latency >= Subject$Low_HDI & Subject$Latency <= Subject$High_HDI, 1,0))
mean(ifelse(Subject$Latency[Subject$TrialTimeDay < 3] >= Subject$Low_HDI[Subject$TrialTimeDay < 3] & Subject$Latency[Subject$TrialTimeDay < 3] <= Subject$High_HDI[Subject$TrialTimeDay < 3], 1,0))
mean(ifelse(Subject$Latency[Subject$TrialTimeDay > 3] >= Subject$Low_HDI[Subject$TrialTimeDay > 3] & Subject$Latency[Subject$TrialTimeDay > 3] <= Subject$High_HDI[Subject$TrialTimeDay > 3], 1,0))
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------


par(mfrow=c(1,1))
plot(0,0)
legend(-1,1, 
       legend=c("Human Data", "Calibration", "Prediction"), 
       col=c(1,2,4), lwd=3, lty=c(1,3,3), pch=c(20,8,8),
       horiz=T)
#---------------------------------------------
plot(0,0)
legend(-1,1, 
       legend=c("Tenison's Change Points", "PPE Change Points"), 
       col=c(1,4), lwd=3, lty=c(1,3), pch=c(20,20),
       horiz=F)
