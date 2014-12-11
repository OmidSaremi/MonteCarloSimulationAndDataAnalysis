# A more refined function which computes autocorrelations. It is used with the scripts which handle 
# large volumes of Monte Carlo data

autocorrelation <-function(timeRange, numOfMeasurements, numPerFile, pathToFolder, fileNameRoot){
  
  x<-vector()
  normalizedAutoCor<-vector()
  for(i in 1: numOfMeasurements){
    x[i]<-getX(i, numPerFile, pathToFolder, fileNameRoot)  
  }
  for(timeSeparation in 0:timeRange){
  sp<-numOfMeasurements-timeSeparation
  tp<-timeSeparation+1
  x1<-x[1:sp]
  x2<-x[tp:numOfMeasurements]

  x1<-as.numeric(x1)
  x2<-as.numeric(x2)
  normalizedAutoCor[timeSeparation+1]<-(mean(x1*x2)-mean(x1)*mean(x2))/sd(x)^2
  } 
  integratedAutoCor<- sum(normalizedAutoCor[2:timeRange+1])
  result<-list()
  result[[1]]<-integratedAutoCor
  result[[2]]<-normalizedAutoCor
  
return (result)
}