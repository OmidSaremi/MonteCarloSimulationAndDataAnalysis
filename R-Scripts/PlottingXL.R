# This scripts extracts the Monte Carlo Polyakov Loop Susceptibility data 
# from its data file and computes it variability using the Jackknife procedure to set the error-bars.
# At the same time, it fits a second-order polynomial to the data in the wanted range. 
# Computes the prediction error of the regression model to set the 68% confidence region 
# around the regression curve.
# The mean values and corresponding error-bars are visualized using ggplot2 package


library(ggplot2)
library(Hmisc)
source('./autocorrelationOfMCMCTimeSeries.R')
source('./jackknifeXL.R')
source('./estimateError.R')
source('./ChiSquareFit.R')
source('./predictChiSquareFit.R')

len <- 9
numb <- as.character(len)
numb1 <- '9'
myLatticeSize <- 90

st <- 6        #3          #4            #4            #               #5 
en <- 9        #7 N=7 L6   #7  N=7 L=7   #8 N=7 L=7.5  #   N=7 L=8     #9 N=7 L=9

fileName <- paste("XL_N7_L", numb1,"_disc250_final", sep="")
out1 <- paste("FinalVersion_XL_N7_L", numb1, "_Zoom.pdf", sep="")
out2 <- paste("FinalVersion_XL_N7_L", numb1,".pdf", sep="")
myLabel <- paste("N=7 L=",numb, sep="")

f <- readLines(paste('./MCDataToAnalyze/XL_N7_disc250/', fileName, sep=""))
f <- strsplit(f, "\ |\t")
autocor<-autocorrelationOfMCMCTimeSeries(f, 2000, 50, 150)
tau_integ <- autocor[[2]]
lambda <- autocor[[4]]
jackknifed <- jackknifeXL(lambda, tau_integ, f, myLatticeSize)

meanVals<-jackknifed$Mean
errors<-jackknifed$JackknifeError

ind <- sort.int(lambda, index.return=TRUE)

lambda <- lambda[ind$ix]
meanVals<-meanVals[ind$ix]
errors<-errors[ind$ix]


fitInfo<- estimateError(lambda[st:en], meanVals[st:en], errors[st:en], 1000000)
print(fitInfo)
predict.object <- predictChiSquareFit(min(lambda[st:en]), max(lambda[st:en]), 50, fitInfo$meanPoints, fitInfo$stdErr, fitInfo$cov)

pdf(out1)

errbar(lambda[st:en], meanVals[st:en], meanVals[st:en]+errors[st:en], meanVals[st:en]-errors[st:en], axes=F,pch=17, col='red', xlab='', ylab='', ylim=c(0, 0.86))
axis(1, cex.axis=0.95, las=0)
axis(2, cex.axis=0.95, las=2)
mtext(side=1, expression(lambda), col='purple', cex=2, las=1, line=2.5)
mtext(side=2, expression(chi), col='purple', cex=2, las=2, line=2.5)
text(x=c(1.12), y=c(0.8), labels=myLabel, col='blue', family="serif")
box(col='blue')

lines(predict.object$lambdaRange, predict.object$prediction, col='darkgreen', pch=23, cex=0.6)
lines(predict.object$lambdaRange, predict.object$lowerConfidenceCurve, col='blue', pch=17, cex=0.4, lty=2)
lines(predict.object$lambdaRange, predict.object$upperConfidenceCurve, col='blue', pch=17, cex=0.4, lty=2)
polygon(c(predict.object$lambdaRange,rev(predict.object$lambdaRange)),c(predict.object$upperConfidenceCurve,rev( predict.object$lowerConfidenceCurve)), col=rgb(0.7,0.8,230/255, 0.3))
legend("topright", c("CI","Fit"), col = c("blue", "darkgreen"), lty=c(1,1))

dev.off()

maxUpper <- predict.object$lambdaRange[which.max(predict.object$upperConfidenceCurve)]
maxUpperValue <- max(predict.object$upperConfidenceCurve)

maxLower <- predict.object$lambdaRange[which.max(predict.object$lowerConfidenceCurve)]
maxLowerValue <- max(predict.object$lowerConfidenceCurve)

maxLambda <-predict.object$lambdaRange[which.max(predict.object$prediction)]
maxValue <- max(predict.object$prediction)

print(" LambdaMax and its error bar: ")

print(maxUpper)
print(maxLambda)
print(maxLower)

print("The maximum values are: ")

print(maxUpperValue)
print(maxLowerValue)
print(maxValue)

# calculate the Chi-squared error of the fit

coefs<-fitInfo$meanPoints
print(coefs[1]*lambda[st:en]^2+coefs[2]*lambda[st:en]+coefs[3]-meanVals[st:en])
print(errors[st:en]^2)
chiSq<-sum((coefs[1]*lambda[st:en]^2+coefs[2]*lambda[st:en]+coefs[3]-meanVals[st:en])^2/errors[st:en]^2)
print("The chi-squared value is: ")
print(chiSq)

pdf(out2)

ggplot(data.frame(lambda, meanVals),aes(y=meanVals, x=lambda)) +  
  geom_point(colour='blue', size=2, pch=23, bg='red') + geom_errorbar(aes(ymin=meanVals-errors, ymax=meanVals+errors), colour='red', width=0.01, size=0.5)+ 
  theme(panel.background = element_rect(color='blue', fill = 'white'))+ 
  labs(x=expression(lambda), y=expression(chi), size=10)+theme(axis.title.x = element_text(color='purple',size = 20, angle = 00), axis.title.y = element_text(color='purple', size = 20, angle = 00))+
  annotate("text", label=myLabel, x=0.95, y=0.73, family="serif", fontface="plain", color='blue')

dev.off()

