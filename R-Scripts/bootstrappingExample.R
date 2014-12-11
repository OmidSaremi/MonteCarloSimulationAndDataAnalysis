# This script performs bootstrap resampling on my data set. This is one of the methods used 
# in this project to estimate variability of the computed averages (i.e., the error-bars). 

bootstrapMean<-function(data, numResampling=100){
   meanVec<-vector()
   bootstrapMean<-vector()
   originalMean<-vector()
   numRows<-nrow(data)
   dataSize<-ncol(data)
 for(i in 1:numRows){ 
   for(iter in 1:numResampling){
      resamplingIndices<-sample(dataSize, replace=TRUE);
      aRow<-data[i,]
      resampled<-aRow[resamplingIndices]
      meanVec[iter]<-mean(as.numeric(resampled))
   }
   originalMean[i]<-mean(as.numeric(data[i,]))
   bootstrapMean[i]<-mean(meanVec)
}

output<-list(bootstrapMean, originalMean, meanVec)
return (output)
}

