# A simple rescaling of the data, needed for computing the quark density from
# spatial average of |TrU|^2

f<-readLines('./MCDataToAnalyze/OldData/N14_L25_M125_Order_30_measur_45000_disc_50.txt')
f<-strsplit(f, "\ |\t")
numMeas<-length(f[[1]])-1
numObs<-length(f)
newMeas<-matrix(0, nrow=numObs, ncol=numMeas+1)


for(i in 1:numObs){
  var<-as.numeric(f[[i]][numMeas+1])
  print(var)
  n <-as.numeric(f[[i]][1:numMeas])
  n<-var*n
  newMeas[i, ]<-c(n, var)
}
write.table(newMeas, append=T, row.names=F, col.names=F, file="./MCDataToAnalyze/OldData/New_N14_L25_M125_Order_30_measur_45000_disc_50.txt")