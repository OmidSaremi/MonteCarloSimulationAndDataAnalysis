# This is one of the main functions were an exploratory analysis of the autocorrelation data is performed. 
# Using autocor function, it plots the log of the normalized 
# autocorrelation function versus the Monte Carlo time. This is to inspect if at short-time 
# separations it can be modeled as a exponential. It also computes the integrated autocorralation.  

source('./autocor.R')
MCTimeSample<-vector()
normalizedAutoCorrelation<-vector()

for(i in 0:100){
  MCTimeSample[i+1]<-i
  normalizedAutoCorrelation[i+1]<- autocor(timeSeparation=i, fileName="EnergyMeasurementResults0.txt", numRows=12000, mySample=12000)
}

autoCor<-as.data.frame(cbind(MCTimeSample, normalizedAutoCorrelation)) #dataframe were autocor is stored 
logOfAutoCor<-as.data.frame(cbind(MCTimeSample, log(normalizedAutoCorrelation))) # dataframe where time samples and log of autocor are stored
colnames(logOfAutoCor)=c("time", "logOf" )
logOfAutoCor<-logOfAutoCor[complete.cases(logOfAutoCor),]


#Find outliers

library(outliers)
outlier_tf = outlier(logOfAutoCor$logOf,logical=TRUE)
find_outlier = which(outlier_tf==TRUE,arr.ind=TRUE)
data_new = logOfAutoCor[-find_outlier,]

 

par(mfrow=c(1,2))
plot(autoCor$MCTimeSample, autoCor$normalizedAutoCorrelation, pch=10, col="blue")
plot(MCTimeSample, normalizedAutoCorrelation, pch=10, col="blue")
linearModel<-lm(data_new$logOf ~ data_new$time)

plot(data_new$time, data_new$logOf, pch=10, col="blue")
abline(linearModel)

# Integrated autocorrlation

len<-length(normalizedAutoCorrelation)
subsetAuto<-normalizedAutoCorrelation[2:len]
print(sum(subsetAuto)+0.5)
