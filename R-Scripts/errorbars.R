# In this scripts error bars are computes and are returned as an array

prepareErrorBarsAndData<-function(pathToFile, timeRange){
  
  l<-autocorrelationOfMCMCTimeSeries(pathToFile, timeRange)
  
  averages<-l[[1]]
  standardDeviations<-l[[2]]
  dfOfAutoCors<-l[[3]]
  numOfMCMeasurements<-l[[4]]
  indepVar<-l[[5]]
  
  integratedAutoCors<-apply(dfOfAutoCors[, 2:timeRange+1], 1, sum)
  
  standardDeviationEstimates<-standardDeviations*sqrt((1.0+2.0*integratedAutoCors)/numOfMCMeasurements)
  yplus<-averages+standardDeviationEstimates
  yminus<-averages-standardDeviationEstimates
 
  return (list(indepVar, averages, yplus, yminus))
}