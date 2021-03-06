setwd("D:/Sand_Lance/Habitat_Model_Data")
SL_TimeLag<-read.csv("SL_Time_Lag_Calanus_Pred_Data_200612.csv", header=TRUE)
#get rid of anything after 2008 for sand lance 
SL_TimeLag<-SL_TimeLag[1:32,]

#load in the processed copepod data
library(plyr)
CSI<-read.csv("NES_CSI_190930.csv", header=TRUE)
CSI$Year<-CSI$row.names
#setting the years to match up with sand lance (eg make zoops from 1977 match sl from 1980 for 3 yr lag)
CSI$Three_Year_Lag<-CSI$Year+3

#set a four year lag for spring survey 
CSI$Four_Year_Lag<-CSI$Year+4
#now calculate rolling means for everything 
Zoops_CSI<-CSI[,3:14]
library(zoo)
Rolling_means<-NULL
for (i in 1:ncol(Zoops_CSI)){ 
  #take a 2 yr "rolling mean" (really midpoint)
  Rolling_means[[i]]<-rollmean(Zoops_CSI[,i],2)
}

Rolling_Cope_means<-as.data.frame(Rolling_means, col.names=names(Zoops_CSI))

#assign them their "lagged" years, 3-4 yrs for zoops
Rolling_Cope_means$Lagged_year<-seq(1981, 2018, by=1)

SL_Spring_Correlations<-merge(SL_TimeLag, Rolling_Cope_means, by.x=c("Year"), by.y=c("Lagged_year"))

#parsing dataframe for rolling means
Herring_mack<-SL_TimeLag[,names(SL_TimeLag) %in% c("Herr_Ind", "Year")]
Herring_mack$Three_Year_Lag<-Herring_mack$Year+3

#now do the rolling mean for 2-3 year lag to get at larval predation 
Herring_Roll_mean<-rollmean(Herring_mack$Herr_Ind, 2)
Three_Year_Lag<-seq(1980, 2010, by=1)
#attach years to the rolling means 
Herring_Mack_Roll<-cbind(Herring_Roll_mean,Three_Year_Lag)

#add these to the bottom trawl with the three year lag 
SL_Spring_Full<-merge(SL_Spring_Correlations, Herring_Mack_Roll, by.x=c("Year"), by.y=c("Three_Year_Lag"))

library(ecodata)#super useful package from NEFSC
WSW<-slopewater[slopewater$Var=="WSW proportion ne channel",]
SL_Slopewater<-merge(SL_Spring_Full, WSW, by.x=c("Year"), by.y=c("Time"))
#chop off the NAs
SL_Slopewater<-SL_Slopewater[1:25,]
SL_Slopewater$PropWSW<-SL_Slopewater$Value/100#make it a proportion

slopewater_glm<-lm(SL_SP_Anom~GOM_Spring_Cal+Herring_Roll_mean+PropWSW, data=SL_Slopewater, na.action = "na.fail")

summary(slopewater_glm)

#######predictions##########

#now for projections
stepbetaprediction<- function (Time){
  L_tau<-ifelse(Time<=2009, 0, 1)
  Mu<-exp(-18.8799732099464+0.00981681192228245*Time+0.260410052373209*(Time-2009)*L_tau)/(1+exp(-18.8799732099464+0.00981681192228245*Time+0.260410052373209*(Time-2009)*L_tau))
  return(Mu)}
#this is the changepoint beta regression with the fitted parameter values from Matlab

#subset warm slope water to be that before 2008 
WSW_old<-WSW[WSW$Time<=2008,]

#get the estimates from the beta regression
Time<-seq(from=2020, to=2100, by=1)
Fitted_Mu<-stepbetaprediction(Time)-0.0000001
Mu_Pre_2008<-mean(WSW_old$Value)/100
Cope_Mean<-mean(Rolling_Cope_means$GOM_Spring_Cal)
Cope_Var<-sd(Rolling_Cope_means$GOM_Spring_Cal)

#to get the final values for each RCP scenario 
End_Cope_mean_8.5=log(0.5)
End_Cope_mean_4.5=log(0.68)

#make it a constant, linearly decreasing mean
Increment_8.5=(End_Cope_mean_8.5-Cope_Mean)/81
Increment_4.5=(End_Cope_mean_4.5-Cope_Mean)/81

#make this the change in value with time
Seqs<-seq(1,81, by=1)
Calanus_Time_4.5<-Cope_Mean+Increment_4.5*Seqs
Calanus_Time_8.5<-Cope_Mean+Increment_8.5*Seqs

#herring in low abundance from 1973 to 1988
Low_Herring<-c(-1.2726,-0.3906 ,-0.5673,-0.2997,-1.6673,-1.0049,-0.7824,-1.0498,-1.2079,-1.2498,-1.1129,-0.5398,-1.407,-0.867,-0.6402,-0.395)
Low_Herring_mean<-mean(Low_Herring)
Low_Herring_var<-sd(Low_Herring)
#set average to a mean of 0 and standard deviation of 1 
Avg_herring<-0
Avg_Herring_var<-1

#now preallocate variable for each scenario
N=81
SL_Prediction_Plaus_4.5<-data.frame(PropWSW=double(N),GOM_Spring_Cal=double(N),Herring_Roll_mean=double(N))
SL_Prediction_Plaus_8.5<-data.frame(PropWSW=double(N),GOM_Spring_Cal=double(N),Herring_Roll_mean=double(N))

SL_Prediction_Red_WSW_4.5<-data.frame(PropWSW=double(N),GOM_Spring_Cal=double(N),Herring_Roll_mean=double(N))
SL_Prediction_Red_WSW_8.5<-data.frame(PropWSW=double(N),GOM_Spring_Cal=double(N),Herring_Roll_mean=double(N))

SL_Prediction_Avg_herr_4.5<-data.frame(PropWSW=double(N),GOM_Spring_Cal=double(N),Herring_Roll_mean=double(N))
SL_Prediction_Avg_herr_8.5<-data.frame(PropWSW=double(N),GOM_Spring_Cal=double(N),Herring_Roll_mean=double(N))

SL_Prediction_Avg_herr_WSW_4.5<-data.frame(PropWSW=double(N),GOM_Spring_Cal=double(N),Herring_Roll_mean=double(N))
SL_Prediction_Avg_herr_WSW_8.5<-data.frame(PropWSW=double(N),GOM_Spring_Cal=double(N),Herring_Roll_mean=double(N))

#preallocate fitted values for each scenario
SL_Predicted_Values_Plaus_4.5<-data.frame(matrix("",ncol=1000,nrow=81))
SL_Predicted_Values_Plaus_8.5<-data.frame(matrix("",ncol=1000,nrow=81))

#preallocate values to fill in for the time series
Plausible_4.5_Calanus<-data.frame(matrix("",ncol=1000,nrow=81))
Plausible_4.5_Herring<-data.frame(matrix("",ncol=1000,nrow=81))
Plausible_4.5_WSW



#same set of preallocation, reduced WSW
SL_Predicted_Values_Red_WSW_4.5<-data.frame(matrix("",ncol=1000,nrow=81))
SL_Predicted_Values_Red_WSW_8.5<-data.frame(matrix("",ncol=1000,nrow=81))

#avg herr
SL_Predicted_Values_Avg_Herr_4.5<-data.frame(matrix("",ncol=1000,nrow=81))
SL_Predicted_Values_Avg_Herr_8.5<-data.frame(matrix("",ncol=1000,nrow=81))

#optimistic
SL_Predicted_Values_Avg_Herr_WSW_4.5<-data.frame(matrix("",ncol=1000,nrow=81))
SL_Predicted_Values_Avg_Herr_WSW_8.5<-data.frame(matrix("",ncol=1000,nrow=81))

#create a matrix to track the postive or negative anomalies for sand lance & herring
SL_Pos_neg_Plaus_4.5<-as.matrix(data.frame(matrix("",ncol=1000,nrow=81)))
SL_Pos_neg_Plaus_8.5<-as.matrix(data.frame(matrix("",ncol=1000,nrow=81)))
SL_Pos_neg_Plaus_4.5_herring<-as.matrix(data.frame(matrix("",ncol=1000,nrow=81)))
SL_Pos_neg_Plaus_8.5_herring<-as.matrix(data.frame(matrix("",ncol=1000,nrow=81)))

#for reduced wsw
SL_Pos_neg_Red_WSW_4.5<-as.matrix(data.frame(matrix("",ncol=1000,nrow=81)))
SL_Pos_neg_Red_WSW_8.5<-as.matrix(data.frame(matrix("",ncol=1000,nrow=81)))
SL_Pos_neg_Red_WSW_4.5_herring<-as.matrix(data.frame(matrix("",ncol=1000,nrow=81)))
SL_Pos_neg_Red_WSW_8.5_herring<-as.matrix(data.frame(matrix("",ncol=1000,nrow=81)))

#for average herring
SL_Pos_neg_Avg_Herr_4.5<-as.matrix(data.frame(matrix("",ncol=1000,nrow=81)))
SL_Pos_neg_Avg_Herr_8.5<-as.matrix(data.frame(matrix("",ncol=1000,nrow=81)))
SL_Pos_neg_Avg_Herr_4.5_herring<-as.matrix(data.frame(matrix("",ncol=1000,nrow=81)))
SL_Pos_neg_Avg_Herr_8.5_herring<-as.matrix(data.frame(matrix("",ncol=1000,nrow=81)))

#for optimistic 
SL_Pos_neg_Avg_Herr_WSW_4.5<-as.matrix(data.frame(matrix("",ncol=1000,nrow=81)))
SL_Pos_neg_Avg_Herr_WSW_8.5<-as.matrix(data.frame(matrix("",ncol=1000,nrow=81)))
SL_Pos_neg_Avg_Herr_WSW_4.5_herring<-as.matrix(data.frame(matrix("",ncol=1000,nrow=81)))
SL_Pos_neg_Avg_Herr_WSW_8.5_herring<-as.matrix(data.frame(matrix("",ncol=1000,nrow=81)))


for (i in 1:1000){
  #find alpha and beta from mean* dispersion parameter (comes from fitted regression in Matlab)
  Alpha_high<-Fitted_Mu*8.81315270996855
  Beta_high<-(1-Fitted_Mu)*8.81315270996855
  
  #random beta variable following this pattern
  SL_Prediction_Plaus_4.5$PropWSW<-rbeta(81, Alpha_high, Beta_high)
  SL_Prediction_Plaus_8.5$PropWSW<-rbeta(81, Alpha_high, Beta_high)
  
  #random decreasing calanus 
  SL_Prediction_Plaus_4.5$GOM_Spring_Cal<-rnorm(81,mean=Calanus_Time_4.5, sd=Cope_Var)
  SL_Prediction_Plaus_8.5$GOM_Spring_Cal<-rnorm(81,mean=Calanus_Time_8.5, sd=Cope_Var)
  
  #random herring
  SL_Prediction_Plaus_4.5$Herring_Roll_mean<-rnorm(81,mean=Low_Herring_mean, sd=Low_Herring_var)
  SL_Prediction_Plaus_8.5$Herring_Roll_mean<-rnorm(81,mean=Low_Herring_mean, sd=Low_Herring_var)
  
  #predicted values
  SL_Predicted_Values_Plaus_4.5[,i]<-predict(slopewater_glm, SL_Prediction_Plaus_4.5)
  SL_Predicted_Values_Plaus_8.5[,i]<-predict(slopewater_glm, SL_Prediction_Plaus_8.5)
  
  #reduced wsw
  #finding alpha and beta based on mean*dispersion parameter 
  Alpha_Red<-Mu_Pre_2008*7.97349252723927
  Beta_Red<-(1-Mu_Pre_2008)*7.97349252723927
  
  SL_Prediction_Red_WSW_4.5$PropWSW<-rbeta(81, Alpha_Red, Beta_Red)
  SL_Prediction_Red_WSW_8.5$PropWSW<-rbeta(81, Alpha_Red, Beta_Red)
  
  SL_Prediction_Red_WSW_4.5$GOM_Spring_Cal<-rnorm(81,mean=Calanus_Time_4.5, sd=Cope_Var)
  SL_Prediction_Red_WSW_8.5$GOM_Spring_Cal<-rnorm(81,mean=Calanus_Time_8.5, sd=Cope_Var)
  
  SL_Prediction_Red_WSW_4.5$Herring_Roll_mean<-rnorm(81,mean=Low_Herring_mean, sd=Low_Herring_var)
  SL_Prediction_Red_WSW_8.5$Herring_Roll_mean<-rnorm(81,mean=Low_Herring_mean, sd=Low_Herring_var)
  
  SL_Predicted_Values_Red_WSW_4.5[,i]<-predict(slopewater_glm, SL_Prediction_Red_WSW_4.5)
  SL_Predicted_Values_Red_WSW_8.5[,i]<-predict(slopewater_glm, SL_Prediction_Red_WSW_8.5)
  
  #repeat for average herring scenario
  SL_Prediction_Avg_herr_4.5$PropWSW<-rbeta(81, Alpha_high, Beta_high)
  SL_Prediction_Avg_herr_8.5$PropWSW<-rbeta(81, Alpha_high, Beta_high)
  
  SL_Prediction_Avg_herr_4.5$GOM_Spring_Cal<-rnorm(81,mean=Calanus_Time_4.5, sd=Cope_Var)
  SL_Prediction_Avg_herr_8.5$GOM_Spring_Cal<-rnorm(81,mean=Calanus_Time_8.5, sd=Cope_Var)
  
  SL_Prediction_Avg_herr_4.5$Herring_Roll_mean<-rnorm(81,mean=Avg_herring, sd=Avg_Herring_var)
  SL_Prediction_Avg_herr_8.5$Herring_Roll_mean<-rnorm(81,mean=Avg_herring, sd=Avg_Herring_var)
  
  SL_Predicted_Values_Avg_Herr_4.5[,i]<-predict(slopewater_glm, SL_Prediction_Avg_herr_4.5)
  SL_Predicted_Values_Avg_Herr_8.5[,i]<-predict(slopewater_glm, SL_Prediction_Avg_herr_8.5)
  
  #avg herring and avg wsw 
  #average herring, what we now call the optimistic scenario
  Alpha_Red<-Mu_Pre_2008*7.97349252723927
  Beta_Red<-(1-Mu_Pre_2008)*7.97349252723927
  
  SL_Prediction_Avg_herr_WSW_4.5$PropWSW<-rbeta(81, Alpha_Red, Beta_Red)
  SL_Prediction_Avg_herr_WSW_8.5$PropWSW<-rbeta(81, Alpha_Red, Beta_Red)
  
  SL_Prediction_Avg_herr_WSW_4.5$GOM_Spring_Cal<-rnorm(81,mean=Calanus_Time_4.5, sd=Cope_Var)
  SL_Prediction_Avg_herr_WSW_8.5$GOM_Spring_Cal<-rnorm(81,mean=Calanus_Time_8.5, sd=Cope_Var)
  
  SL_Prediction_Avg_herr_WSW_4.5$Herring_Roll_mean<-rnorm(81,mean=Avg_herring, sd=Avg_Herring_var)
  SL_Prediction_Avg_herr_WSW_8.5$Herring_Roll_mean<-rnorm(81,mean=Avg_herring, sd=Avg_Herring_var)
  
  SL_Predicted_Values_Avg_Herr_WSW_4.5[,i]<-predict(slopewater_glm, SL_Prediction_Avg_herr_WSW_4.5)
  SL_Predicted_Values_Avg_Herr_WSW_8.5[,i]<-predict(slopewater_glm, SL_Prediction_Avg_herr_WSW_8.5)
  
  for (k in 1:81){
    #'flag' if negative anomaly
    SL_Pos_neg_Plaus_4.5[k,i]<-ifelse(SL_Predicted_Values_Plaus_4.5[k,i]<=1.45,1,0)
    SL_Pos_neg_Plaus_8.5[k,i]<-ifelse(SL_Predicted_Values_Plaus_8.5[k,i]<=1.45,1,0)
    #same for herring
    SL_Pos_neg_Plaus_4.5_herring[k,i]<-ifelse(SL_Prediction_Plaus_4.5$Herring_Roll_mean[k]<=0,1,0)
    SL_Pos_neg_Plaus_8.5_herring[k,i]<-ifelse(SL_Prediction_Plaus_8.5$Herring_Roll_mean[k]<=0,1,0)
    
    SL_Pos_neg_Red_WSW_4.5[k,i]<-ifelse(SL_Predicted_Values_Red_WSW_4.5[k,i]<=1.45,1,0)
    SL_Pos_neg_Red_WSW_8.5[k,i]<-ifelse(SL_Predicted_Values_Red_WSW_8.5[k,i]<=1.45,1,0)
    
    SL_Pos_neg_Red_WSW_4.5_herring[k,i]<-ifelse(SL_Prediction_Red_WSW_4.5$Herring_Roll_mean[k]<=0,1,0)
    SL_Pos_neg_Red_WSW_8.5_herring[k,i]<-ifelse(SL_Prediction_Red_WSW_8.5$Herring_Roll_mean[k]<=0,1,0)
    
    
    SL_Pos_neg_Avg_Herr_4.5[k,i]<-ifelse(SL_Predicted_Values_Avg_Herr_4.5[k,i]<=1.45,1,0)
    SL_Pos_neg_Avg_Herr_8.5[k,i]<-ifelse(SL_Predicted_Values_Avg_Herr_8.5[k,i]<=1.45,1,0)
    
    SL_Pos_neg_Avg_Herr_4.5_herring[k,i]<-ifelse(SL_Prediction_Avg_herr_4.5$Herring_Roll_mean[k]<=0,1,0)
    SL_Pos_neg_Avg_Herr_8.5_herring[k,i]<-ifelse(SL_Prediction_Avg_herr_8.5$Herring_Roll_mean[k]<=0,1,0)
    
    
    SL_Pos_neg_Avg_Herr_WSW_4.5[k,i]<-ifelse(SL_Predicted_Values_Avg_Herr_WSW_4.5[k,i]<=1.45,1,0)
    SL_Pos_neg_Avg_Herr_WSW_8.5[k,i]<-ifelse(SL_Predicted_Values_Avg_Herr_WSW_8.5[k,i]<=1.45,1,0)
    
    SL_Pos_neg_Avg_Herr_WSW_4.5_herring[k,i]<-ifelse(SL_Prediction_Avg_herr_WSW_4.5$Herring_Roll_mean[k]<=0,1,0)
    SL_Pos_neg_Avg_Herr_WSW_8.5_herring[k,i]<-ifelse(SL_Prediction_Avg_herr_WSW_8.5$Herring_Roll_mean[k]<=0,1,0)
    
  }
}



SL_Pos_Plaus_4.5<-as.data.frame(SL_Pos_neg_Plaus_4.5)
SL_Pos_Plaus_4.5_herring<-as.data.frame(SL_Pos_neg_Plaus_4.5_herring)

SL_Pos_Plaus_8.5<-as.data.frame(SL_Pos_neg_Plaus_8.5)
SL_Pos_Plaus_8.5_herring<-as.data.frame(SL_Pos_neg_Plaus_8.5_herring)

SL_Pos_Red_WSW_4.5<-as.data.frame(SL_Pos_neg_Red_WSW_4.5)
SL_Pos_Red_WSW_4.5_herring<-as.data.frame(SL_Pos_neg_Red_WSW_4.5_herring)

SL_Pos_Red_WSW_8.5<-as.data.frame(SL_Pos_neg_Red_WSW_8.5)
SL_Pos_Red_WSW_8.5_herring<-as.data.frame(SL_Pos_neg_Red_WSW_8.5_herring)


SL_Pos_Avg_Herr_4.5<-as.data.frame(SL_Pos_neg_Avg_Herr_4.5)
SL_Pos_Avg_Herr_4.5_herring<-as.data.frame(SL_Pos_neg_Avg_Herr_4.5_herring)

SL_Pos_Avg_Herr_8.5<-as.data.frame(SL_Pos_neg_Avg_Herr_8.5)
SL_Pos_Avg_Herr_8.5_herring<-as.data.frame(SL_Pos_neg_Avg_Herr_8.5_herring)

SL_Pos_Avg_Herr_WSW_4.5<-as.data.frame(SL_Pos_neg_Avg_Herr_WSW_4.5)
SL_Pos_Avg_Herr_WSW_4.5_herring<-as.data.frame(SL_Pos_neg_Avg_Herr_WSW_4.5_herring)

SL_Pos_Avg_Herr_WSW_8.5<-as.data.frame(SL_Pos_neg_Avg_Herr_WSW_8.5)
SL_Pos_Avg_Herr_WSW_8.5_herring<-as.data.frame(SL_Pos_neg_Avg_Herr_WSW_8.5_herring)


for (i in 1:1000){
  SL_Pos_Plaus_4.5[,i]<-as.numeric(as.character(SL_Pos_Plaus_4.5[,i]))
  SL_Pos_Plaus_8.5[,i]<-as.numeric(as.character(SL_Pos_Plaus_8.5[,i]))
  
  SL_Pos_Plaus_4.5_herring[,i]<-as.numeric(as.character(SL_Pos_Plaus_4.5_herring[,i]))
  SL_Pos_Plaus_8.5_herring[,i]<-as.numeric(as.character(SL_Pos_Plaus_8.5_herring[,i]))
  
  SL_Pos_Red_WSW_4.5[,i]<-as.numeric(as.character(SL_Pos_Red_WSW_4.5[,i]))
  SL_Pos_Red_WSW_8.5[,i]<-as.numeric(as.character(SL_Pos_Red_WSW_8.5[,i]))
  
  SL_Pos_Red_WSW_4.5_herring[,i]<-as.numeric(as.character(SL_Pos_Red_WSW_4.5_herring[,i]))
  SL_Pos_Red_WSW_8.5_herring[,i]<-as.numeric(as.character(SL_Pos_Red_WSW_8.5_herring[,i]))
  
  
  SL_Pos_Avg_Herr_4.5[,i]<-as.numeric(as.character(SL_Pos_Avg_Herr_4.5[,i]))
  SL_Pos_Avg_Herr_8.5[,i]<-as.numeric(as.character(SL_Pos_Avg_Herr_8.5[,i]))
  
  SL_Pos_Avg_Herr_4.5_herring[,i]<-as.numeric(as.character(SL_Pos_Avg_Herr_4.5_herring[,i]))
  SL_Pos_Avg_Herr_8.5_herring[,i]<-as.numeric(as.character(SL_Pos_Avg_Herr_8.5_herring[,i]))
  
  SL_Pos_Avg_Herr_WSW_4.5[,i]<-as.numeric(as.character(SL_Pos_Avg_Herr_WSW_4.5[,i]))
  SL_Pos_Avg_Herr_WSW_8.5[,i]<-as.numeric(as.character(SL_Pos_Avg_Herr_WSW_8.5[,i]))
  
  SL_Pos_Avg_Herr_WSW_4.5_herring[,i]<-as.numeric(as.character(SL_Pos_Avg_Herr_WSW_4.5_herring[,i]))
  SL_Pos_Avg_Herr_WSW_8.5_herring[,i]<-as.numeric(as.character(SL_Pos_Avg_Herr_WSW_8.5_herring[,i]))
  
}
#now calculate proportion of years with negative sand lance anomalies 
#values will change (minimally) when re-run given these are simulations
mean(colSums(SL_Pos_Plaus_4.5))/81
#0.963
mean(colSums(SL_Pos_Plaus_8.5))/81
#0.981

mean(colSums(SL_Pos_Red_WSW_4.5))/81
#0.342
mean(colSums(SL_Pos_Red_WSW_8.5))/81
#0.466

mean(colSums(SL_Pos_Avg_Herr_4.5))/81
#0.967
mean(colSums(SL_Pos_Avg_Herr_8.5))/81
#0.981

mean(colSums(SL_Pos_Avg_Herr_WSW_4.5))/81
#0.640
mean(colSums(SL_Pos_Avg_Herr_WSW_8.5))/81
#0.718

SL_Herr_Negative_4.5<-floor((SL_Pos_Plaus_4.5+SL_Pos_Plaus_4.5_herring)/2)
SL_Herr_Negative_8.5<-floor((SL_Pos_Plaus_8.5+SL_Pos_Plaus_8.5_herring)/2)

SL_negative_counts_Plaus_4.5<-colSums(SL_Herr_Negative_4.5)
SL_negative_counts_Plaus_8.5<-colSums(SL_Herr_Negative_8.5)

#proportion of neg anomalies for both
#SL_Neg_Mean_Plaus_4.5<-
mean(SL_negative_counts_Plaus_4.5)/81
#95.0
#SL_Neg_Mean_Plaus_8.5
mean(SL_negative_counts_Plaus_8.5)/81
#96.7

##### same for REDWSW
SL_Herr_Negative_4.5_RedWSW<-floor((SL_Pos_Red_WSW_4.5+SL_Pos_Red_WSW_4.5_herring)/2)
SL_Herr_Negative_8.5_RedWSW<-floor((SL_Pos_Red_WSW_8.5+SL_Pos_Red_WSW_8.5_herring)/2)

SL_negative_counts_Red_WSW_4.5<-colSums(SL_Herr_Negative_4.5_RedWSW)
SL_negative_counts_Red_WSW_8.5<-colSums(SL_Herr_Negative_8.5_RedWSW)
#SL_Neg_Mean_Red_WSW_4.5<-
#  SL_Neg_Mean_Red_WSW_8.5<-
mean(SL_negative_counts_Red_WSW_4.5)/81
#33.1
mean(SL_negative_counts_Red_WSW_8.5)/81
#45.5

###########same for avg herring
SL_Herr_Negative_4.5_herring<-floor((SL_Pos_Avg_Herr_4.5+SL_Pos_Avg_Herr_4.5_herring)/2)
SL_Herr_Negative_8.5_herring<-floor((SL_Pos_Avg_Herr_8.5+SL_Pos_Avg_Herr_8.5_herring)/2)

SL_negative_counts_Avg_Herr_4.5<-colSums(SL_Herr_Negative_4.5_herring)
SL_negative_counts_Avg_Herr_8.5<-colSums(SL_Herr_Negative_8.5_herring)

SL_Neg_Mean_Avg_Herr_4.5<-mean(SL_negative_counts_Avg_Herr_4.5)/81
mean(SL_negative_counts_Avg_Herr_4.5)/81
#46.6
SL_Neg_Mean_Avg_Herr_8.5<-mean(SL_negative_counts_Avg_Herr_8.5)/81
mean(SL_negative_counts_Avg_Herr_8.5)/81
#48.2

########same for optimistic 

SL_Herr_Negative_WSW_4.5_herring<-floor((SL_Pos_Avg_Herr_WSW_4.5+SL_Pos_Avg_Herr_WSW_4.5_herring)/2)
SL_Herr_Negative_WSW_8.5_herring<-floor((SL_Pos_Avg_Herr_WSW_8.5+SL_Pos_Avg_Herr_WSW_8.5_herring)/2)

SL_negative_counts_Avg_Herr_WSW_4.5<-colSums(SL_Herr_Negative_WSW_4.5_herring)
SL_negative_counts_Avg_Herr_WSW_8.5<-colSums(SL_Herr_Negative_WSW_8.5_herring)

SL_Neg_Mean_Avg_Herr_WSW_4.5<-mean(SL_negative_counts_Avg_Herr_WSW_4.5)/81
mean(SL_negative_counts_Avg_Herr_WSW_4.5)/81
#19.9
SL_Neg_Mean_Avg_Herr_WSW_8.5<-mean(SL_negative_counts_Avg_Herr_WSW_8.5)/81
mean(SL_negative_counts_Avg_Herr_WSW_8.5)/81
#25.9
 
SL_neg_counts_Sl_only_4.5<-mean(colSums(SL_Pos_Avg_Herr_WSW_4.5))/81
SL_neg_counts_Sl_only_8.5<-mean(colSums(SL_Pos_Avg_Herr_WSW_8.5))/81

SL_Neg_sd_Plaus_8.5<-sd(SL_negative_counts_Plaus_4.5)
SL_Neg_sd_Plaus_4.5<-sd(SL_negative_counts_Plaus_8.5)

#find mean sl value and interquartile range

SL_Predicted_Values_Estimates_Plaus_4.5<-rowMeans(SL_Predicted_Values_Plaus_4.5)
SL_Predicted_Values_Estimates_Plaus_8.5<-rowMeans(SL_Predicted_Values_Plaus_8.5)


SL_Quants_Plaus_4.5<-t(apply(SL_Predicted_Values_Plaus_4.5, 1, quantile))
SL_Quants_Plaus_8.5<-t(apply(SL_Predicted_Values_Plaus_8.5, 1, quantile))

SL_75_Plaus_4.5<-SL_Quants_Plaus_4.5[,4]
SL_25_Plaus_4.5<-SL_Quants_Plaus_4.5[,2]

SL_75_Plaus_8.5<-SL_Quants_Plaus_8.5[,4]
SL_25_Plaus_8.5<-SL_Quants_Plaus_8.5[,2]


#same for reduced WSW

SL_Predicted_Values_Estimates_Red_WSW_4.5<-rowMeans(SL_Predicted_Values_Red_WSW_4.5)
SL_Predicted_Values_Estimates_Red_WSW_8.5<-rowMeans(SL_Predicted_Values_Red_WSW_8.5)


SL_Quants_Red_WSW_4.5<-t(apply(SL_Predicted_Values_Red_WSW_4.5, 1, quantile))
SL_Quants_Red_WSW_8.5<-t(apply(SL_Predicted_Values_Red_WSW_8.5, 1, quantile))

SL_75_Red_WSW_4.5<-SL_Quants_Red_WSW_4.5[,4]
SL_25_Red_WSW_4.5<-SL_Quants_Red_WSW_4.5[,2]

SL_75_Red_WSW_8.5<-SL_Quants_Red_WSW_8.5[,4]
SL_25_Red_WSW_8.5<-SL_Quants_Red_WSW_8.5[,2]

###average herring
SL_Predicted_Values_Estimates_Avg_Herr_4.5<-rowMeans(SL_Predicted_Values_Avg_Herr_4.5)
SL_Predicted_Values_Estimates_Avg_herr_8.5<-rowMeans(SL_Predicted_Values_Avg_Herr_8.5)


SL_Quants_Avg_Herr_4.5<-t(apply(SL_Predicted_Values_Avg_Herr_4.5, 1, quantile))
SL_Quants_Avg_Herr_8.5<-t(apply(SL_Predicted_Values_Avg_Herr_8.5, 1, quantile))

SL_75_Avg_Herr_4.5<-SL_Quants_Avg_Herr_4.5[,4]
SL_25_Avg_Herr_4.5<-SL_Quants_Avg_Herr_4.5[,2]

SL_75_Avg_Herr_8.5<-SL_Quants_Avg_Herr_8.5[,4]
SL_25_Avg_Herr_8.5<-SL_Quants_Avg_Herr_8.5[,2]

#### optimistic
SL_Predicted_Values_Estimates_Avg_Herr_WSW_4.5<-rowMeans(SL_Predicted_Values_Avg_Herr_WSW_4.5)
SL_Predicted_Values_Estimates_Avg_Herr_WSW_8.5<-rowMeans(SL_Predicted_Values_Avg_Herr_WSW_8.5)


SL_Quants_Avg_Herr_WSW_4.5<-t(apply(SL_Predicted_Values_Avg_Herr_WSW_4.5, 1, quantile))
SL_Quants_Avg_Herr_WSW_8.5<-t(apply(SL_Predicted_Values_Avg_Herr_WSW_8.5, 1, quantile))

SL_75_Avg_Herr_WSW_4.5<-SL_Quants_Avg_Herr_WSW_4.5[,4]
SL_25_Avg_Herr_WSW_4.5<-SL_Quants_Avg_Herr_WSW_4.5[,2]

SL_75_Avg_Herr_WSW_8.5<-SL_Quants_Avg_Herr_WSW_8.5[,4]
SL_25_Avg_Herr_WSW_8.5<-SL_Quants_Avg_Herr_WSW_8.5[,2]

# now precent decreases
((exp(SL_Predicted_Values_Estimates_Plaus_4.5[81])-exp(SL_Predicted_Values_Estimates_Plaus_4.5[1]))/(exp(SL_Predicted_Values_Estimates_Plaus_4.5[1])))*100
#47.5%, corrected to 51.6
((exp(SL_Predicted_Values_Estimates_Plaus_8.5[81])-exp(SL_Predicted_Values_Estimates_Plaus_8.5[1]))/(exp(SL_Predicted_Values_Estimates_Plaus_8.5[1])))*100
#74.2%, corrected to 76.9

((exp(SL_Predicted_Values_Estimates_Red_WSW_4.5[81])-exp(SL_Predicted_Values_Estimates_Red_WSW_4.5[1]))/(exp(SL_Predicted_Values_Estimates_Red_WSW_4.5[1])))*100
#44.5% decrease, corrected to 38.5
((exp(SL_Predicted_Values_Estimates_Red_WSW_8.5[81])-exp(SL_Predicted_Values_Estimates_Red_WSW_8.5[1]))/(exp(SL_Predicted_Values_Estimates_Red_WSW_8.5[1])))*100
#70.8% decrease, corrected to 74.8

((exp(SL_Predicted_Values_Estimates_Avg_Herr_4.5[81])-exp(SL_Predicted_Values_Estimates_Avg_Herr_4.5[1]))/(exp(SL_Predicted_Values_Estimates_Avg_Herr_4.5[1])))*100
#48.2%, corrected to 52.0
((exp(SL_Predicted_Values_Estimates_Avg_herr_8.5[81])-exp(SL_Predicted_Values_Estimates_Avg_herr_8.5[1]))/(exp(SL_Predicted_Values_Estimates_Avg_herr_8.5[1])))*100
#75.4% decrease, corrected to 75.1

((exp(SL_Predicted_Values_Estimates_Avg_Herr_WSW_4.5[81])-exp(SL_Predicted_Values_Estimates_Avg_Herr_WSW_4.5[1]))/(exp(SL_Predicted_Values_Estimates_Avg_Herr_WSW_4.5[1])))*100
#40.9, corrected to 40.4
((exp(SL_Predicted_Values_Estimates_Avg_Herr_WSW_8.5[81])-exp(SL_Predicted_Values_Estimates_Avg_Herr_WSW_8.5[1]))/(exp(SL_Predicted_Values_Estimates_Avg_Herr_WSW_8.5[1])))*100
#70.5%, corrected to 71.4

