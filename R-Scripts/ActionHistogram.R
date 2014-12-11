# Extracts the action data from data file and creates action histograms using ggplot2 package. 
# The output is saves as pdf.


action<-readLines('./MCDataToAnalyze/ActionN12L13Final121_128.txt')      
numOfObservations<-length(action)
numOfMCMeasurements<-length(action[[1]])-1


cents<-vector()
indi <- seq(1, 50, by=1)
for(i in indi ){cents[i]<-as.numeric(action[[i]][numOfMCMeasurements+1]) }
centIndex<-sort.int(cents, index.return=TRUE)$ix
print(centIndex)

MCTime<-seq(1, 100000, by=1)

for(i in centIndex[14:14]){       " The desirable range is 12:26, 19 and 22 good for presentation"
   simulatedLambda=as.numeric(action[[i]][numOfMCMeasurements+1])
   s<-as.character(simulatedLambda)
   print(simulatedLambda)
   S<-as.numeric(action[[i]][1:numOfMCMeasurements])
   df <-data.frame(MCTime, S[MCTime])
}
   pdf('./ActionHistogramLambda123.pdf') 
   ggplot(df, aes(x=S[MCTime]))+geom_histogram(aes(y=..density..), fill='skyblue', color='blue', binwidth=25)+ 
   labs(x="S", y="P(S)")+theme(axis.title.x = element_text(color='purple',size = 15, angle = 00, family="serif"), axis.title.y = element_text(color='purple', size = 15, angle = 00, family="serif")) 
   dev.off()

