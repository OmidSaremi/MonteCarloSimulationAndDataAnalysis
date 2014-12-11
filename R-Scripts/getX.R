# Very large Monte Carlo data sets are kept in multiple files (to 
# compute autocorrelations, reading the
# entire data set into the memory as a single file is not an option. 
# This function decides which of the smaller data files need to be loaded


getX<-function(index, numPerFile, pathToFolder, fileNameRoot){
  s<-floor((index-1)/numPerFile)
  if(index==segment*numPerFile+1 ){
    rm(x, envir=globalenv())
    x<<-list()
    s1<-paste(pathToFolder, "/",fileNameRoot,as.character(s),".txt", sep="")
    s2<-paste(pathToFolder, "/",fileNameRoot,as.character(s+1),".txt", sep="")
    x[[1]]<<-readLines(s1)
    x[[2]]<<-readLines(s2) 
    segment<<-segment+2
    return (x[[s %% 2+1 ]][(index-1) %% numPerFile+1])
  }
  
  else{
    return (x[[s %% 2 + 1 ]][(index-1) %% numPerFile+1])
  }
}