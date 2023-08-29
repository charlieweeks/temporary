# This Script is meant to run the the PPE over data collected by Tension and Anderson (unknown year)
# The goal of this script is to run the PPE the data and attempt to estimate 
# Cognitive Stages of learning 
#--------------------------------------------------------------------------------------------------------------------------------
#Load packages
library(dplyr)
library(R2jags)
library(Hmisc)
#--------------------------------------------------------------------------------------------------------------------------------
#Load Functions
{
  gammaShRaFromMeanSD = function( mean , sd ) {
    if ( mean <=0 ) stop("mean must be > 0")
    if ( sd <=0 ) stop("sd must be > 0")
    shape = mean^2/sd^2
    rate = mean/sd^2
    return( list( shape=shape , rate=rate ) )
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
  Parameter_Format=function(Name,Num_Sub, Data){
    
    Data$Parameter = Name
    for(i in 1:Num_Sub){
      assign("x", Data[ , c(i,5)]  )
      x$Participant=i
      colnames(x)=c("Post"," Parameter", "Subject")
      if(i ==1 ) { y = x }  else { y = rbind(y, x) }
      if(i == Num_Sub) { y$Participant =  rep( unique( MEG_Group$Subject_ID[MEG_Group$Session == unique(MEG_Group$Session)[Num_Groups]] ), each = length(x$Post) ) }  
    }  
    
    return( y )
  }
  
  Extract_Information=function(samples){
    
    b           = samples$BUGSoutput$sims.list$b
    m           = samples$BUGSoutput$sims.list$m
    a           = samples$BUGSoutput$sims.list$a
    tau         = samples$BUGSoutput$sims.list$tau
    tauPrime    = samples$BUGSoutput$sims.list$tauPrime
    w           = samples$BUGSoutput$sims.list$w
    w_rank      = samples$BUGSoutput$sims.list$w_rank
    wTmp        = samples$BUGSoutput$sims.list$wTmp
    Prediction  = samples$BUGSoutput$sims.list$Prediction
    OS_Prediction  = samples$BUGSoutput$sims.lis$OS_Prediction
    wTmp_round          = apply(wTmp, 2, median)
    w_round_mean        = apply(   w, 2, mean)
    w_round_median      = apply(   w, 2, median)
    w_round_mode        = apply(   w, 2, getmode)
    #------------------------------------------------------------------    
    Weight_Set_5 = length(w[w==5])/length(w)
    Weight_Set_4 = length(w[w==4])/length(w)
    Weight_Set_3 = length(w[w==3])/length(w)
    Weight_Set_2 = length(w[w==2])/length(w)
    Weight_Set_1 = length(w[w==1])/length(w)
    #------------------------------------------------------------------  
    #tauPrime = 
    #-----------------------------------------------------------------------  
    PPE     = ((1 - apply(Prediction, 2, mean) ) * max(Data$ComputeDuration)) 
    PPE_HDI = ((1 - apply(Prediction, 2, HDIofMCMC) ) * max(Data$ComputeDuration)) 
    #------------------------------------------------------------------------------------------------
    OS_PPE_1      = ((1 - apply(OS_Prediction[,,1], 2, mean) ) * max(Data$ComputeDuration) ) 
    OS_PPE_2      = ((1 - apply(OS_Prediction[,,2], 2, mean) ) * max(Data$ComputeDuration) ) 
    OS_PPE_3      = ((1 - apply(OS_Prediction[,,3], 2, mean) ) * max(Data$ComputeDuration) ) 
    OS_PPE_4      = ((1 - apply(OS_Prediction[,,4], 2, mean) ) * max(Data$ComputeDuration) ) 
    OS_PPE_5      = ((1 - apply(OS_Prediction[,,5], 2, mean) ) * max(Data$ComputeDuration) ) 
    #------------------------------------------------------------------------------------------------
    HDI_OS_PPE_1  = ((1 - apply(OS_Prediction[,,1], 2, HDIofMCMC ) ) * max(Data$ComputeDuration) ) 
    HDI_OS_PPE_2  = ((1 - apply(OS_Prediction[,,2], 2, HDIofMCMC ) ) * max(Data$ComputeDuration) ) 
    HDI_OS_PPE_3  = ((1 - apply(OS_Prediction[,,3], 2, HDIofMCMC ) ) * max(Data$ComputeDuration ) ) 
    HDI_OS_PPE_4  = ((1 - apply(OS_Prediction[,,4], 2, HDIofMCMC ) ) * max(Data$ComputeDuration ) ) 
    HDI_OS_PPE_5  = ((1 - apply(OS_Prediction[,,5], 2, HDIofMCMC ) ) * max(Data$ComputeDuration ) ) 
    #------------------------------------------------------------------------------------------------
    Overall_Average        = ((1 - apply( OS_Prediction[,,], 2, mean) ) * max(Data$ComputeDuration) ) 
    Overall_Average_HDI    = ((1 - apply( OS_Prediction[,,], 2, HDIofMCMC) ) * max(Data$ComputeDuration) )
    #_________________________________________________________________________________________________________     
    Weighted_Average  =  (1 - apply(
      rbind(
        apply( OS_Prediction[,,5], 2, sample, Weight_Set_5 * length(OS_Prediction[,1,1]) ),
        apply( OS_Prediction[,,4], 2, sample, Weight_Set_4 * length(OS_Prediction[,1,1]) ),
        apply( OS_Prediction[,,3], 2, sample, Weight_Set_3 * length(OS_Prediction[,1,1]) ),
        apply( OS_Prediction[,,2], 2, sample, Weight_Set_2 * length(OS_Prediction[,1,1]) ),
        apply( OS_Prediction[,,1], 2, sample, Weight_Set_1 * length(OS_Prediction[,1,1]) )
      ), 2, mean ) ) * max(Data$ComputeDuration)
    Weighted_Average_HDI  =  (1 - apply(
      rbind(
        apply( OS_Prediction[,,5], 2, sample, Weight_Set_5 * length(OS_Prediction[,1,1]) ),
        apply( OS_Prediction[,,4], 2, sample, Weight_Set_4 * length(OS_Prediction[,1,1]) ),
        apply( OS_Prediction[,,3], 2, sample, Weight_Set_3 * length(OS_Prediction[,1,1]) ),
        apply( OS_Prediction[,,2], 2, sample, Weight_Set_2 * length(OS_Prediction[,1,1]) ),
        apply( OS_Prediction[,,1], 2, sample, Weight_Set_1 * length(OS_Prediction[,1,1]) )
      ), 2, HDIofMCMC ) ) * max(Data$ComputeDuration)
    #_________________________________________________________________________________________________________     
    Weighted_Median  =  (1 - apply(
      rbind(
        apply( OS_Prediction[,,5], 2, sample, (length(w_round_median[w_round_median == 5])/length(w_round_median))  * length(OS_Prediction[,1,1]) ),
        apply( OS_Prediction[,,4], 2, sample, (length(w_round_median[w_round_median == 4])/length(w_round_median))  * length(OS_Prediction[,1,1]) ),
        apply( OS_Prediction[,,3], 2, sample, (length(w_round_median[w_round_median == 3])/length(w_round_median))  * length(OS_Prediction[,1,1]) ),
        apply( OS_Prediction[,,2], 2, sample, (length(w_round_median[w_round_median == 2])/length(w_round_median))  * length(OS_Prediction[,1,1]) ),
        apply( OS_Prediction[,,1], 2, sample, (length(w_round_median[w_round_median == 1])/length(w_round_median))  * length(OS_Prediction[,1,1]) )
      ), 2,mean ) ) * max(Data$ComputeDuration)
    
    Weighted_Median_HDI  =  (1 - apply(
      rbind(
        apply( OS_Prediction[,,5], 2, sample, (length(w_round_median[w_round_median == 5])/length(w_round_median))  * length(OS_Prediction[,1,1])  ),
        apply( OS_Prediction[,,4], 2, sample, (length(w_round_median[w_round_median == 4])/length(w_round_median))  * length(OS_Prediction[,1,1]) ),
        apply( OS_Prediction[,,3], 2, sample, (length(w_round_median[w_round_median == 3])/length(w_round_median))  * length(OS_Prediction[,1,1]) ),
        apply( OS_Prediction[,,2], 2, sample, (length(w_round_median[w_round_median == 2])/length(w_round_median))  * length(OS_Prediction[,1,1]) ),
        apply( OS_Prediction[,,1], 2, sample, (length(w_round_median[w_round_median == 1])/length(w_round_median))  * length(OS_Prediction[,1,1]) )
      ), 2, HDIofMCMC ) ) * max(Data$ComputeDuration)
    #_________________________________________________________________________________________________________     
    Weighted_Mode  =  (1 - apply(
      rbind(
        apply( OS_Prediction[,,5], 2, sample, (length(w_round_mode[w_round_mode == 5])/length(w_round_mode))  * length(OS_Prediction[,1,1]) ),
        apply( OS_Prediction[,,4], 2, sample, (length(w_round_mode[w_round_mode == 4])/length(w_round_mode))  * length(OS_Prediction[,1,1]) ),
        apply( OS_Prediction[,,3], 2, sample, (length(w_round_mode[w_round_mode == 3])/length(w_round_mode))  * length(OS_Prediction[,1,1]) ),
        apply( OS_Prediction[,,2], 2, sample, (length(w_round_mode[w_round_mode == 2])/length(w_round_mode))  * length(OS_Prediction[,1,1]) ),
        apply( OS_Prediction[,,1], 2, sample, (length(w_round_mode[w_round_mode == 1])/length(w_round_mode))  * length(OS_Prediction[,1,1]) )
      ), 2,mean ) ) * max(Data$ComputeDuration)
    
    Weighted_Mode_HDI  =  (1 - apply(
      rbind(
        apply( OS_Prediction[,,5], 2, sample, (length(w_round_mode[w_round_mode == 5])/length(w_round_mode))  * length(OS_Prediction[,1,1]) ),
        apply( OS_Prediction[,,4], 2, sample, (length(w_round_mode[w_round_mode == 4])/length(w_round_mode))  * length(OS_Prediction[,1,1]) ),
        apply( OS_Prediction[,,3], 2, sample, (length(w_round_mode[w_round_mode == 3])/length(w_round_mode))  * length(OS_Prediction[,1,1]) ),
        apply( OS_Prediction[,,2], 2, sample, (length(w_round_mode[w_round_mode == 2])/length(w_round_mode))  * length(OS_Prediction[,1,1]) ),
        apply( OS_Prediction[,,1], 2, sample, (length(w_round_mode[w_round_mode == 1])/length(w_round_mode))  * length(OS_Prediction[,1,1]) )
      ), 2, HDIofMCMC ) ) * max(Data$ComputeDuration)
    #_________________________________________________________________________________________________________     
    #--------------------------------------------------
    #--------------------------------------------------
    
    
    
    #--------------------------------------------------  
    #--------------------------------------------------      
    Predictions = data.frame(
      PPE_1 = as.numeric( c(PPE, OS_PPE_1) ),
      PPE_2 = as.numeric( c(PPE, OS_PPE_2) ),
      PPE_3 = as.numeric( c(PPE, OS_PPE_3) ),
      PPE_4 = as.numeric( c(PPE, OS_PPE_4) ),
      PPE_5 = as.numeric( c(PPE, OS_PPE_5) ),
      #-------------------------------------
      H_HDI_PPE_1 = as.numeric( c(PPE_HDI[1,], HDI_OS_PPE_1[1,]) ),
      H_HDI_PPE_2 = as.numeric( c(PPE_HDI[1,], HDI_OS_PPE_2[1,]) ),
      H_HDI_PPE_3 = as.numeric( c(PPE_HDI[1,], HDI_OS_PPE_3[1,]) ),
      H_HDI_PPE_4 = as.numeric( c(PPE_HDI[1,], HDI_OS_PPE_4[1,]) ),
      H_HDI_PPE_5 = as.numeric( c(PPE_HDI[1,], HDI_OS_PPE_5[1,]) ),
      #----------------------------------
      L_HDI_PPE_1 = as.numeric( c(PPE_HDI[2,], HDI_OS_PPE_1[2,]) ),
      L_HDI_PPE_2 = as.numeric( c(PPE_HDI[2,], HDI_OS_PPE_2[2,]) ),
      L_HDI_PPE_3 = as.numeric( c(PPE_HDI[2,], HDI_OS_PPE_3[2,]) ),
      L_HDI_PPE_4 = as.numeric( c(PPE_HDI[2,], HDI_OS_PPE_4[2,]) ),
      L_HDI_PPE_5 = as.numeric( c(PPE_HDI[2,], HDI_OS_PPE_5[2,]) ),
      #-------------------------------------
      Overall_Average  = as.numeric( c(PPE, Overall_Average )  ),
      Weighted_Average = as.numeric( c(PPE, Weighted_Average)  ),
      Weighted_Median  = as.numeric( c(PPE, Weighted_Median )  ),
      Weighted_Mode    = as.numeric( c(PPE, Weighted_Mode   )  ),
      #---------------------------------------------------
      H_HDI_Overall_Average  = as.numeric( c(PPE_HDI[1,], Overall_Average_HDI[1,]  )   ),
      H_HDI_Weighted_Average = as.numeric( c(PPE_HDI[1,], Weighted_Average_HDI[1,] )   ),
      H_HDI_Weighted_Median  = as.numeric( c(PPE_HDI[1,], Weighted_Median_HDI[1,]  )   ),
      H_HDI_Weighted_Mode    = as.numeric( c(PPE_HDI[1,], Weighted_Mode_HDI[1,]    )   ),
      #---------------------------------------------------   
      L_HDI_Overall_Average  = as.numeric( c(PPE_HDI[2,] , Overall_Average_HDI[2,] )  ),
      L_HDI_Weighted_Average = as.numeric( c(PPE_HDI[2,] , Weighted_Average_HDI[2,]) ),
      L_HDI_Weighted_Median  = as.numeric( c(PPE_HDI[2,] , Weighted_Median_HDI[2,] ) ),
      L_HDI_Weighted_Mode    = as.numeric( c(PPE_HDI[2,] , Weighted_Mode_HDI[2,]   ) ),
      #---------------------------------------------------
      Cal_Point   =     as.numeric(Num_Cal),
      Tau_Prime_1 = as.numeric(getmode(tauPrime[,1] )  ),
      Tau_Prime_2 = as.numeric(getmode(tauPrime[,2] )  ),
      Tau_Prime_3 = as.numeric(getmode(tauPrime[,3] )  ),
      Tau_Prime_4 = as.numeric(getmode(tauPrime[,4] )  ),
      w_round_mean   = as.numeric( c(w_round_mean, rep(0, Num_Rounds-Num_Cal  ) )  ),
      w_round_median = as.numeric( c(w_round_median,rep(0, Num_Rounds-Num_Cal ) )  ),
      w_round_mode   = as.numeric( c(w_round_mode,  rep(0, Num_Rounds-Num_Cal ) )  ),
      Weight_Set_5   = as.numeric( Weight_Set_5 ),
      Weight_Set_4   = as.numeric( Weight_Set_4 ),
      Weight_Set_3   = as.numeric( Weight_Set_3 ),
      Weight_Set_2   = as.numeric( Weight_Set_2 ),
      Weight_Set_1   = as.numeric( Weight_Set_1 )
      #-------------------------------------
    )
    return(Predictions)
  } 
}
#--------------------------------------------------------------------------------------------------------------------------------
#Load Data
setwd( paste(dirname(rstudioapi::getSourceEditorContext()$path), sep="") ) 
Data = read.csv("Example_Data.csv")
#--------------------------------------------------------------------------------------------------------------------------------
#Create a unique variable for each participant
Data$ID = paste(Data$Participant, Data$ID)

#Creat a name for the output file
Behavior_File_Name = "PPE_Change_Detection_2_Subject_Full_Estimation"
#Data = Data[Data$Participant %in% unique(Data$Participant)[2], ]
Data  = Data[Data$Trained %in% c(4,  8, 16, 32), ]

# Now you can run run the following loop over the data.
# KEEP in mind this model does take ahwile to run, so if you're just trying to get it started then maybe don't run it over the full data right away
#------------------------------------------------------------------------------------------------------------------------------------
for(i in 1 : length( unique(Data$ID) ) ) {
  i=1
  print("####################################################################################")
  print(paste(i, " out of ", length(unique(Data$ID)) ))
  print("####################################################################################")
#--------------------------------------------------------------------------------
  #Select data for a specific participant
  Subject       = Data[Data$ID == unique(Data$ID)[i],]
  Subject$Experiment_Instance = 1 : length(Subject$Instances)
  #--------------------------------------------------------------------------------
  #Load the data that you are interested in
  Num_Sub                        = 1 #Number of items the model is being run over
  #------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  #Get the learninge curve for each of the thee curves
  #------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  {
  Num_Cal                    =  ifelse(length(Subject$Instances) == 45, 25, 12 )            # Determine the number of trials that the model will calibrate
  Num_Rounds                 =  length(Subject$Instances)                                   # Number of total Rounds
  Perf                       =  1 - ( Subject$ComputeDuration / max(Data$ComputeDuration) ) # Normalize the Performance Variable
  Perf                       =  ifelse(Perf == 0, .001, Perf)                               # Modify the lower bound
  #Perf_2                     = Perf                                
  #Perf_3                     = Perf
  #Perf_4                     = Perf
  #Perf_5                     = Perf
  N                          =  Subject$Instances - 1                                       # N Variable for PPE Model
  T_var                      =  Subject$T                                                   # T Variable 
  CumAveInvLogTime_diff_lag1 =  Subject$CumAveInvLogTime_diff_lag1                          # CumAveInvLogTime_diff_lag1 variable for decay variable
  gamma                      =  4                                                           # Determine the number of possible change points
  Spike_Slabe_Prior          = c( Num_Cal / (2 * Num_Cal - 1), rep(1 / (2 * Num_Cal - 1), Num_Cal - 1) ) # Defind spike and slab prior 
  #--------------------------------------------------------------------
  #Here we tell JAGS what Variables it needs to run the model
  data = c("Num_Cal", "Num_Rounds", "Perf", "N", "T_var",
           "CumAveInvLogTime_diff_lag1", "gamma", 
           "Spike_Slabe_Prior")
  #, "Perf_2", "Perf_3", "Perf_4"
  #--------------------------------------------------------------------------------
  #Set Parameters to be saved from JAGS
  parameters <- c("b", "m", "a", "tau", "k_par", "tauPrime", 
                  "w", "wTmp", "Prediction", "pred", "OS_Prediction") 
  #Determine initial parameters values
  myinits <-	list(
    list( b = rep(.01, gamma+1), m = rep(.15, gamma+1), a = rep(.15, gamma+1), tau = rep(.5, gamma+1), k_par=rep(999, gamma+1), tauPrime = c(1,5,10,25) ) ,
    list( b = rep(.03, gamma+1), m = rep(.15, gamma+1), a = rep(.15, gamma+1), tau = rep(.5, gamma+1), k_par=rep(999, gamma+1), tauPrime = c(1,5,10,25) ) ,
    list( b = rep(.02, gamma+1), m = rep(.15, gamma+1), a = rep(.15, gamma+1), tau = rep(.5, gamma+1), k_par=rep(999, gamma+1), tauPrime = c(1,5,10,25) )
  )
  # Initial group assignment
  #----------------------------------------------------------------------------------------------------
  #----------------------------------------------------------------------------------------------------
  #Jags Model
  #--------------------------------------------------------------------------------
  setwd( paste(dirname(rstudioapi::getSourceEditorContext()$path), sep="") ) 
  timestamp()
  samples <- jags(data, inits=myinits, parameters,
                  model.file ="PPE_Simple_Change_Detection_Cal_Pred.txt",
                  n.chains = 3, 
                  n.iter   = 6000, 
                  n.burnin = 2000, 
                  n.thin   = 1,      DIC=T )
  timestamp()
  #---------------------------------------------------------------------------------------
  #---------------------------------------------------------------------------------------
  #Extract output from first Model
  Predictions = Extract_Information(samples)
  #---------------------------------------------------------------------------------------
  #---------------------------------------------------------------------------------------
  Subject = cbind(Subject, Predictions)

  }
  #------------------------------------------------------------------------------------------------------------
  #State Estimation Model 
  #----------------------
  Subject$Stage_Marker = Subject$w_round_median
  for(j in 1 : length( unique( Subject$Stage_Marker ) ) ){
    Set = Subject[Subject$Stage_Marker == unique(Subject$Stage_Marker)[j], ]
    #------------------------------------------------------------------------------------
    #================================================================================================================
    DV         = Set$ComputeDuration
    N          = length(DV)
    #================================================================================================================
    data       =  c("DV", "N")
    parameters =  c("Stage", "RT_Prediction",
                    "Cog_mean" , "Ass_mean" , "Pro_mean",
                    "Cog_SD"   , "Ass_SD"   , "Pro_SD",
                    "Cog_Prior", "Ass_Prior", "Pro_Prior",
                    "Cog_Shape", "Ass_Shape", "Pro_Shape",
                    "Cog_Rate" , "Ass_Rate" , "Pro_Rate"
    )
    myinits     =  list(
      list(Stage = c(1), Cog_mean=8469.779,  Ass_mean = 3146.525,  Pro_mean =  1334.677 ),
      list(Stage = c(2), Cog_mean=8469.779,  Ass_mean = 3146.525,  Pro_mean =  1334.677 ),
      list(Stage = c(3), Cog_mean=8469.779,  Ass_mean = 3146.525,  Pro_mean =  1334.677 )
    )
    #--------------------------------------------------------------------------------------------------------------------
    setwd( paste(dirname(rstudioapi::getSourceEditorContext()$path), sep="") ) 
    Stage_Est <- jags(data, inits=myinits, parameters,
                    model.file ="Simple_Stage_Estimation_Model.txt",
                    n.chains = 3, 
                    n.iter   = 9000, 
                    n.burnin = 2000, 
                    n.thin   = 1,      DIC=T )
    #--------------------------------------------------------------------------------------------------------------------
    #---------------------------------------------------------------------
    Stage             = Stage_Est$BUGSoutput$sims.list$Stage
    Estamated_Stage   = getmode(Stage)
    #----------------------------------------------------------
    #----------------------------------------------------------
    Set$Estimated_Stage = Estamated_Stage
    Set$Stage_1_Freq    = length(Stage[Stage == 1]) 
    Set$Stage_2_Freq    = length(Stage[Stage == 2]) 
    Set$Stage_3_Freq    = length(Stage[Stage == 3]) 
    #----------------------------------------------------------
    #----------------------------------------------------------
    #------------------------------------------------------------------------------------------------
    if(j == 1){ Subject_Estimation = Set } else { Subject_Estimation = rbind(Subject_Estimation, Set) }
  }
  #------------------------------------------------------------------------------------------------------------
  #Extract information from Stage Estimation Model
  {
  Subject = Subject_Estimation  
  #------------------------------------------------------------------------------------------------------------
  #Subject$Estimated_Stage[Subject$Day == 1] 
  #Subject$w_round_median [Subject$Day == 1] 
  #------------------------------------------------------------
  Fourth_Weight = length(Subject$w_round_median [Subject$Day == 1 & Subject$w_round_median == 4]) / length(Subject$w_round_median [Subject$Day == 1 & Subject$w_round_median >= 4 ] ) 
   Fifth_Weight = length(Subject$w_round_median [Subject$Day == 1 & Subject$w_round_median == 5]) / length(Subject$w_round_median [Subject$Day == 1 & Subject$w_round_median >= 4 ] ) 
  #------------------------------------------------------------
  D_Weight = length(Subject$Estimated_Stage[Subject$Day == 1 & Subject$Estimated_Stage == 1]) / length(Subject$Estimated_Stage[Subject$Day == 1])
  A_Weight = length(Subject$Estimated_Stage[Subject$Day == 1 & Subject$Estimated_Stage == 2]) / length(Subject$Estimated_Stage[Subject$Day == 1])
  P_Weight = length(Subject$Estimated_Stage[Subject$Day == 1 & Subject$Estimated_Stage == 3]) / length(Subject$Estimated_Stage[Subject$Day == 1])
  #------------------------------------------------------------
  Declarative = unique( Subject$w_round_median[Subject$Estimated_Stage == 1  & Subject$Day == 1 ] )
  Associaitve = unique( Subject$w_round_median[Subject$Estimated_Stage == 2  & Subject$Day == 1 ] )
  Procedural  = unique( Subject$w_round_median[Subject$Estimated_Stage == 3  & Subject$Day == 1 ] )
  #-----------------------------------------------
  Sample_Extraction = function(Mark, samples){
  x= rbind(
  if( 1 %in% Mark  ) { samples$BUGSoutput$sims.list$OS_Prediction[,,1] },
  if( 2 %in% Mark  ) { samples$BUGSoutput$sims.list$OS_Prediction[,,2] },
  if( 3 %in% Mark  ) { samples$BUGSoutput$sims.list$OS_Prediction[,,3] },
  if( 4 %in% Mark  ) { samples$BUGSoutput$sims.list$OS_Prediction[,,4] },
  if( 5 %in% Mark  ) { samples$BUGSoutput$sims.list$OS_Prediction[,,5] }
  )
  #--------------  
  return(x)
  #--------------  
  }
  #------------------------------------------------------------------------------------------------------  
  Declarative_Samples =  if( length( Declarative ) != 0 ) { (1 - Sample_Extraction(Declarative,  samples) ) * max( Data$ComputeDuration) } else { (1 - Sample_Extraction(1,  samples) ) * 0  }
  Associative_Samples =  if( length( Associaitve ) != 0 ) { (1 - Sample_Extraction(Associaitve,  samples) ) * max( Data$ComputeDuration) } else { (1 - Sample_Extraction(1,  samples) ) * 0  } 
  Procedural_Samples  =  if( length( Procedural  ) != 0 ) { (1 - Sample_Extraction(Procedural ,  samples) ) * max( Data$ComputeDuration) } else { (1 - Sample_Extraction(1,  samples) ) * 0  }
  #-----------------
  Fourth_Sample        =  if( length(  Fourth_Weight  ) != 0 ) { (1 - Sample_Extraction(4 ,  samples) ) * max( Data$ComputeDuration) } else { (1 - Sample_Extraction(1,  samples) ) * 0  }
  Fifth_Sample         =  if( length(  Fifth_Weight  ) != 0 )  { (1 - Sample_Extraction(5 ,  samples) ) * max( Data$ComputeDuration) } else { (1 - Sample_Extraction(1,  samples) ) * 0  }
  #------------------------------------------------------------------------------------------------------
  Declarative_Prediction   = if(  length(Declarative_Samples) != 0 ) {apply(Declarative_Samples, 2, mean)       }  else { rep(0, length(Subject$Instances[Subject$Day== 2]) )   }
  Associative_Prediction   = if(  length(Associative_Samples) != 0 ) {apply(Associative_Samples, 2, mean)       }  else { rep(0, length(Subject$Instances[Subject$Day== 2]) )   }
  Procedural_Prediction    = if(  length(Procedural_Samples)  != 0 ) {apply( Procedural_Samples, 2, mean)       }  else { rep(0, length(Subject$Instances[Subject$Day== 2]) )   }
  Declarative_HDI          = if(  length(Declarative_Samples) != 0 ) {apply(Declarative_Samples, 2, HDIofMCMC ) }  else { rep(0, length(Subject$Instances[Subject$Day== 2]) )   }
  Associative_HDI          = if(  length(Associative_Samples) != 0 ) {apply(Associative_Samples, 2, HDIofMCMC ) }  else { rep(0, length(Subject$Instances[Subject$Day== 2]) )   }
  Procedural_HDI           = if(  length(Procedural_Samples)  != 0 ) {apply( Procedural_Samples, 2, HDIofMCMC ) }  else { rep(0, length(Subject$Instances[Subject$Day== 2]) )   }
  #----------------------------------------------------------------------------------------------
  N_Weighted = apply( rbind(
    rbind(
      
      Fourth_Sample * Fourth_Weight, 
      Fifth_Sample * Fifth_Weight) 
  
    ), 2, mean )
  
  N_Weighted = apply( rbind(
    rbind(
      Fourth_Sample[ 1 : length(Fourth_Sample[,1]) * Fourth_Weight ,],
      Fifth_Sample[ 1 : length(Fifth_Sample[,1]) * Fifth_Weight ,]
      ) 
  ), 2, mean )
  
  Fourth_Sample[ 1 : length(Fourth_Sample[,1]) * Fourth_Weight ,]
   Fifth_Sample[ 1 : length(Fifth_Sample[,1]) * Fifth_Weight ,]
  
  N_Weighted_HDI = apply( rbind(
     rbind(
      Fourth_Sample[ 1 : length(Fourth_Sample[,1]) * Fourth_Weight ,],
      Fifth_Sample[ 1 : length(Fifth_Sample[,1]) * Fifth_Weight ,]
     )  
  ), 2, HDIofMCMC )
  
  Stage_Weighted = apply( rbind(
     rbind(
     
     Declarative_Samples * D_Weight, 
     Associative_Samples * A_Weight, 
     Procedural_Samples  * P_Weight
     
     ) 
    ), 2, mean )
  
  Stage_Weighted_HDI = apply( rbind(
    rbind(
      Declarative_Samples * D_Weight, 
      Associative_Samples * A_Weight, 
      Procedural_Samples  * P_Weight) 
  ), 2, HDIofMCMC )
  
  AP_Back_Weighted = apply( rbind(
    rbind(
      Associative_Samples *   if((A_Weight/(A_Weight + P_Weight) ) == "NaN") {0} else {(A_Weight/(A_Weight + P_Weight) )} , 
      Procedural_Samples  *   if((P_Weight/(A_Weight + P_Weight) ) == "NaN") {0} else {(P_Weight/(A_Weight + P_Weight) )} ) 
      ), 2, mean )
  
  AP_Back_Weighted_HDI = apply( rbind(
    rbind(
      Associative_Samples *   if((A_Weight/(A_Weight + P_Weight) ) == "NaN") {0} else { (A_Weight/(A_Weight + P_Weight) )}, 
      Procedural_Samples  *   if((P_Weight/(A_Weight + P_Weight) ) == "NaN") {0} else { (P_Weight/(A_Weight + P_Weight) )} ) 
      ), 2, HDIofMCMC )
  
  
  DA_Back_Weighted = apply( rbind(
    rbind(
      Declarative_Samples *   if((D_Weight/(D_Weight + A_Weight) ) == "NaN") {0} else { (D_Weight/(D_Weight + A_Weight) )} ,
      Associative_Samples *   if((D_Weight/(D_Weight + A_Weight) ) == "NaN") {0} else {(A_Weight/(D_Weight + A_Weight) )} ) 
  ), 2, mean )
  
  DA_Back_Weighted_HDI = apply( rbind(
    rbind(
      Declarative_Samples *   if((D_Weight/(D_Weight + A_Weight) ) == "NaN") {0} else {  (D_Weight/(D_Weight + A_Weight) ) } ,
      Associative_Samples *   if((D_Weight/(D_Weight + A_Weight) ) == "NaN") {0} else { (A_Weight/(D_Weight + A_Weight ) ) }) 
  ), 2, HDIofMCMC )
  #------------------------------------------------------------------------------------------------------------
  Stage_Predictions =as.data.frame( cbind(
      c(Subject$Weighted_Median[Subject$Day == 1],        N_Weighted),
      c(Subject$L_HDI_Weighted_Median[Subject$Day == 1] , N_Weighted_HDI[1, ]),
      c(Subject$H_HDI_Weighted_Median[Subject$Day == 1] , N_Weighted_HDI[2, ]),
    
      c(Subject$Weighted_Median[Subject$Day == 1],        Declarative_Prediction),
      c(Subject$L_HDI_Weighted_Median[Subject$Day == 1] , Declarative_HDI[1, ]),
      c(Subject$H_HDI_Weighted_Median[Subject$Day == 1] , Declarative_HDI[2, ]),
    
      c(Subject$Weighted_Median[Subject$Day == 1],        Associative_Prediction),
      c(Subject$L_HDI_Weighted_Median[Subject$Day == 1] , Associative_HDI[1, ]),
      c(Subject$H_HDI_Weighted_Median[Subject$Day == 1] , Associative_HDI[2, ]),
    
      c(Subject$Weighted_Median[Subject$Day == 1],        Procedural_Prediction) ,
      c(Subject$L_HDI_Weighted_Median[Subject$Day == 1] , Procedural_HDI[1, ]),
      c(Subject$H_HDI_Weighted_Median[Subject$Day == 1] , Procedural_HDI[2, ]),
    
      c(Subject$Weighted_Median[Subject$Day == 1],        Stage_Weighted          ),  
      c(Subject$L_HDI_Weighted_Median[Subject$Day == 1] , Stage_Weighted_HDI[1, ] ),
      c( Subject$H_HDI_Weighted_Median[Subject$Day == 1] , Stage_Weighted_HDI[2, ]),
  
      c(Subject$Weighted_Median[Subject$Day == 1],        AP_Back_Weighted),
      c(Subject$L_HDI_Weighted_Median[Subject$Day == 1] , AP_Back_Weighted_HDI[1, ]),
      c(Subject$H_HDI_Weighted_Median[Subject$Day == 1] , AP_Back_Weighted_HDI[2, ]),
  
      c(Subject$Weighted_Median[Subject$Day == 1]       , DA_Back_Weighted),
      c(Subject$L_HDI_Weighted_Median[Subject$Day == 1] , DA_Back_Weighted_HDI[1, ] ),
      c(Subject$H_HDI_Weighted_Median[Subject$Day == 1] , DA_Back_Weighted_HDI[2, ] )
  ))
  
  colnames(Stage_Predictions)=c("N_Stage"  , "N_Stage_L_HDI"  ,  "N_Stage_H_HDI"  ,
                                "D_Stage"  , "D_Stage_L_HDI"  ,  "D_Stage_H_HDI"  ,
                                "A_Stage"  , "A_Stage_L_HDI"  ,  "A_Stage_H_HDI"  ,
                                "P_Stage"  , "P_Stage_L_HDI"  ,  "P_Stage_H_HDI"  ,
                                "All_Stage", "All_Stage_L_HDI",  "All_Stage_H_HDI",
                                "AP_Stage" , "AP_Stage_L_HDI" ,  "AP_Stage_H_HDI" ,
                                "DA_Stage" , "DA_Stage_L_HDI" ,  "DA_Stage_H_HDI"  )
  
Subject = cbind(Subject, Stage_Predictions)
  }
  #------------------------------------------------------------------------------------------------------------
  #------------------------------------------------------------------------------------------------------------
  #Out the data 
  setwd( paste(dirname(rstudioapi::getSourceEditorContext()$path), sep="") ) 
  if(i == 1) { write.table(Subject,        Behavior_File_Name , col.names = T) } else { write.table(Subject, Behavior_File_Name , col.names=F, append = TRUE) }
  #if(i == 1) {write.csv2(Parameter_File, paste(Behavior_File_Name, "_Parameters"))} else { write.table(Parameter_File, paste(Behavior_File_Name, "_Parameters"), col.names=F, append = TRUE) }
}
#------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------

