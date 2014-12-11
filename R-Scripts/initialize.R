# Given path to the folder of the project and a file name root, it 
# reads the first two of the files only.

initialize<-function(pathToResultsFolder, fileNameRoot){
  
  fileNames<-list.files(pathToResultsFolder)
  if (length(fileNames)==0) {
      print("Empty Folder!")
      return (0)
}
  else if (length(fileNames)==1){
    x<<-list()  
    completePath0=paste(pathToResultsFolder,"/",fileNameRoot,"0",".txt", sep="")
      x[[1]]<<-readLines(completePath0)
}
  else {
    x<<-list()
    completePath0=paste(pathToResultsFolder,"/",fileNameRoot,"0",".txt", sep="")
    x[[1]]<<-readLines(completePath0)
    completePath1=completePath=paste(pathToResultsFolder,"/",fileNameRoot,"1",".txt", sep="")
    x[[2]]<<-readLines(completePath1)
}
return (x)
}
  