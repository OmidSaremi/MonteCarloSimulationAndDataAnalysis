# This script implements the Jackknife resampling technique. It both estimates the confidence 
# interval (the error-bars) for the Monte Carlo Polyakov Loop Susceptibility data 
# and corrects for the expected bias   

jackknifeXL<-function(indepVar, integratedTau, data, latticeSize){
  
  numOfObservations<-length(data)
  numOfMCMeasurements <- length(data[[1]])-1
  
  theta<-vector()
  sJack<-vector()
  unbiasedMean<-vector()
  
   for(j in 1:numOfObservations){
       tau_int<-floor(integratedTau[j])
       index<-seq(1, numOfMCMeasurements, by=tau_int)
       L<-as.numeric(data[[j]][index])
       Xl_hat<-latticeSize*mean((L-mean(L))^2)
       len<-length(L)
       for(i in 1:len){
       theta[i]<-latticeSize*mean((L[-i]-mean(L[-i]))^2) 
       }
       theta_tilda<-mean(theta[1:len])
       unbiasedMean[j]<-len*Xl_hat-(len-1)*theta_tilda
       sJack[j]<-sqrt((len-1)*mean((theta[1:len]-Xl_hat)^2))
   }  
  return(list(Mean=unbiasedMean, JackknifeError=sJack))    
}