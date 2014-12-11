# This script implements the Jackknife resampling technique. It both estimates the confidence
# interval (the error-bars) for the Monte Carlo quark density data
# and corrects for the expected bias

jackknifeOrder<-function(indepVar, integratedTau, data){
  
  numOfObservations<-length(data)
  #data<-strsplit(data, "\ |\t")
  numOfMCMeasurements <- length(data[[1]])-1
  theta<-vector()
  sJack<-vector()
  unbiasedMean<-vector()
  
   for(j in 1:numOfObservations){
       
      tau_int<-floor(integratedTau[j])
       index<-seq(1, numOfMCMeasurements, by=tau_int)
       L2<-as.numeric(data[[j]][index])
       L2_hat<-mean(L2)
       len<-length(L2)
       for(i in 1:len){
       theta[i]<-mean(L2[-i]) 
       }
       theta_tilda<-mean(theta[1:len])
       #unbiasedMean[j]<-len*L2_hat-(len-1)*theta_tilda
       unbiasedMean[j]<-L2_hat
       sJack[j]<-sqrt((len-1)*mean((theta[1:len]-L2_hat)^2))
   }  
  return(list(Mean=unbiasedMean, JackknifeError=sJack))    
}