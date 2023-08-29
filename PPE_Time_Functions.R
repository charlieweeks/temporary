#--------------------------------------------------------------------------------------------------------------------
# This function takes in data from a particular dataset and creats the time variables needed 
# to run the Predictive Performance Equation (PPE).
# Data needs to be in a preset format. 
#--------------------------------------------------------------------------------------------------------------------
#Data Structure
#             user_id : uniquely identifies a user/subject/participant in a dataset
#         stimulus_id : uniquely identifies a stimulus/item in a dataset
#               time : timestamp associated with each observation in a dataset
#       performance : performance metric associated with each observation (accuracy, percent correct, etc.); must be in range [0, 1]
#      [condition] : experimental condition; if applicable
#  [response_time] : time between stimulus onset and first recorded key press; in seconds.
#       [stimulus] : the stimulus associated with the stimulus_id; e.g., the Japanese-English word-pair that was studied
#            [...] : any additional, dataset-specific columns can be retained
#--------------------------------------------------------------------------------------------------------------------
TrialTime_Sec_Diff                 = function(TrialTimeSec){
    if( length(TrialTimeSec) == 1 ){
      TrialTimeSecDiff = NA
    } else {
      TrialTimeSecDiff = TrialTimeSec - c(NA, TrialTimeSec[1 : (length(TrialTimeSec) - 1) ])
    }
    return( TrialTimeSecDiff )
}

TrialTime_Sec_sum_Diff_Function    = function( TrialTimeSec_sumDiff, TrialTimeSec){
  
    for (i in 1:length( TrialTimeSec ) ) {TrialTimeSec_sumDiff[i] = sum(( TrialTimeSec[i] -  TrialTimeSec[1:i]))}
  
    return(TrialTimeSec_sumDiff)
}


TrialTimeSec_sumDiffneg75_Function = function(TrialTimeSec_sumDiffneg75, TrialTimeSec ){
    
    for (i in 1 : length( TrialTimeSec) ) {
      TrialTimeSec_sumDiffneg75[i] = sum( ( (   ifelse(  TrialTimeSec[i] -  TrialTimeSec[1:i-1] == 0, .01,   TrialTimeSec[i] -  TrialTimeSec[1:i-1] ) )^-0.75 ) ) 
    }
  
    return(TrialTimeSec_sumDiffneg75)
  }

T_Function                         = function(T, TrialTimeSec, TrialTimeSec_sumDiffneg75){  
    for (i in 1:length( TrialTimeSec )){
      T[i] = sum( 
        (  
          ( (  ifelse( TrialTimeSec[i] -  TrialTimeSec[1:i-1] == 0, .01,  TrialTimeSec[i] -  TrialTimeSec[1:i-1] )  ) ^ -0.75 ) /  TrialTimeSec_sumDiffneg75[i] ) * ( ( ifelse( TrialTimeSec[i] -  TrialTimeSec[1:i-1] == 0,.01,  TrialTimeSec[i] -  TrialTimeSec[1:i-1] ) ) ) )
      
    }
    return(T)
}

#time = TrialTime_Sec_Diff(t)
#TrialTimeSec_sumDiff      = rep(0, length(t))
#TrialTime_Sec_sum_Diff    = TrialTime_Sec_sum_Diff_Function(TrialTimeSec_sumDiff, t)
#TrialTimeSec_sumDiffneg75 = rep(0, length(t))
#TrialTimeSec_sumDiffneg75 = TrialTimeSec_sumDiffneg75_Function(TrialTimeSec_sumDiffneg75, t)
#T = rep(0, length(t))
#T = T_Function(T, t, TrialTimeSec_sumDiffneg75)
#plot(T)
#--------------------------------------------------------------------------------------------------------------------
PPE_Time_Variables = function(Data_set){
  require(dplyr)
  
  #----------------------------------------------------
  #----------------------------------------------------  
  a   = group_by(Data_set, ID)
  b   = mutate(a, 
               Instances = 1:length(ID),
               TrialTimeSec = time - min(time),
               TrialTimeSec_Diff = TrialTime_Sec_Diff(TrialTimeSec),
               TrialTimeSec_sumDiff = NA,
               TrialTimeSec_sumDiffneg75 = NA,
               T = NA,
               TrialTimeSec_sumDiff          = TrialTime_Sec_sum_Diff_Function(TrialTimeSec_sumDiff, TrialTimeSec),
               TrialTimeSec_sumDiffneg75     = TrialTimeSec_sumDiffneg75_Function(TrialTimeSec_sumDiffneg75, TrialTimeSec),
               T                             = T_Function (T, TrialTimeSec, TrialTimeSec_sumDiffneg75),
               T = ifelse(T==0,1, T),
               N = Instances-1,
               c = .1,
               s = .1,
               InvLogTime_diff = 1 / (log(TrialTimeSec_Diff + exp(1))),
               #Changed Bout to Trial
               InvLogTime_diff            = ifelse(Instances == 1, 0, InvLogTime_diff),
               CumAveInvLogTime_diff      = ifelse(Instances != 1, cumsum(InvLogTime_diff)/(seq_along(InvLogTime_diff)), 0),
               CumAveInvLogTime_diff_lag1 = c(0, CumAveInvLogTime_diff[ -length(CumAveInvLogTime_diff) ])
  )
  
  d1a = b

  return(d1a)
  
  }
#--------------------------------------------------------------------------------------------------------------------


