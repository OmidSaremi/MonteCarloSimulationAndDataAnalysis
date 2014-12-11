# This function finds the maximum likelihood estimates for the regression line. It also estimates the standard error for the regression coefficients.


estimateErrorLine<-function(x, ymax, ymin, numSamples){

 aCoef<-vector(length=numSamples)
 bCoef<-vector(length=numSamples)
 y <- (ymax+ymin)/2
 sigma <- (ymax-ymin)/2
 
 s0 <- sum(1/sigma^2)
 s1 <- sum(x/sigma^2)
 s2 <- sum(x^2/sigma^2)
 
 
 Delta <- s0*s2-s1^2
 
 for(i in  1 : numSamples){
     newY<-vector(length=length(x))
     for(j in  1: length(x)){
         newY[j]<-rnorm(1, mean=y[j], sd=sigma[j]) 
   }
 # ss<-ChiSquareFit(x, newY, sigma)
 
 b0 <- sum(newY/sigma^2)
 b1 <- sum(x*newY/sigma^2)

 
 aCoef[i] <- (b1*s0-b0*s1)/Delta
 bCoef[i] <- (b0*s2-b1*s1)/Delta
 }
 
 meanA <- mean(aCoef)
 meanB <- mean(bCoef)
 
 return (list(meanCoeffs=c(meanA, meanB), stdErr=c(sd(aCoef), sd(bCoef)))) 
}