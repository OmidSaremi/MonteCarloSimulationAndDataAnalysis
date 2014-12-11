# It compute the autocorrelation function of the Markov chain time series generated in 
# my Monte Carlo simulations. Among other objects, it returns a data frame with
# columns giving the autocorrelations for each time separation 
# up to timeRange. Rows correspond to different observations.  

autocorrelationOfMCMCTimeSeries <-function(pathToFile, timeRange){
 
  normalizedAutoCors<-vector()
  indepVar<-vector()
  averages<-vector()
  standardDeviations<-vector()
  
  x<-readLines(pathToFile)
  
  numOfObservations<-length(x)
  x<-strsplit(x, "\ |\t")
  numOfMCMeasurements<-length(x[[1]])-1
  dfOfAutoCors<-matrix(0, numOfObservations, timeRange+1)
  dfOfIntegratedAutoCors<-matrix(0, numOfObservations, 1)

  for(i in 1:numOfObservations){
     for(timeSeparation in 0:timeRange){
        
        y<-as.numeric(x[[i]][1:numOfMCMeasurements])
        indepVar[i]<-as.numeric(x[[i]][numOfMCMeasurements+1])
        averages[i]<-mean(y)
        standardDeviations[i]<-sd(y)
        
        sp<-numOfMCMeasurements-timeSeparation
        tp<-timeSeparation+1
        x1<-y[1:sp]
        x2<-y[tp:numOfMCMeasurements]
        normalizedAutoCors[timeSeparation+1]<-(mean(x1*x2)-mean(x1)*mean(x2))/sd(y)^2
     }
        dfOfAutoCors[i, ]<-normalizedAutoCors
  }
  return (list(averages, standardDeviations, dfOfAutoCors, numOfMCMeasurements, indepVar))
}