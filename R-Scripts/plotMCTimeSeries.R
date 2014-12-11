# Extracts the "Polyakov Loop data" from data files, collected during Monte Carlo measurements
# and plots the corresponding histograms. 

poly<-readLines('./MCDataToAnalyze/OldData2/MyAction500000_numP100_1202_1206.txt')
poly<-strsplit(poly, "\ |\t")

numOfObservations<-length(poly)
numOfMCMeasurements<-length(poly[[1]])-1
cents<-vector()

for(i in 1:numOfObservations){
  cents[i]<-as.numeric(poly[[i]][numOfMCMeasurements+1])
}
centIndex<-sort.int(cents, index.return=TRUE)$ix

par(mfrow=c(5, 6))
for(i in 1:30){
  x11()
   simulatedLambda<-as.numeric(poly[[i]][numOfMCMeasurements+1])
   print(simulatedLambda)
   s<-as.character(simulatedLambda)
   S<-as.numeric(poly[[i]][1:numOfMCMeasurements])
   hist(S, breaks=100, main=s, col='blue')
}