# This function computes the boundaries of a 68% confidence region for the polynomial regression fit. 
# As input, it takes the regressors, their estimated standard deviations as well as their covariance matix, computed 
# using function estimateErrorLine. 

predictChiSquareFit<-function(from, to, numPoints, FitCoeff, stdCoeff, COVs){
   
   range <- seq(from, to, by=(to-from)/(numPoints-1))
   pred <- FitCoeff[1]*range^2+FitCoeff[2]*range+FitCoeff[3]
   stdErr <- sqrt(stdCoeff[1]^2*range^4+stdCoeff[2]^2*range^2+stdCoeff[3]^2+2*COVs[1]*range^3+2*COVs[2]*range^2+2*COVs[3]*range)
   ymax <- pred+stdErr
   ymin <- pred-stdErr
   
  return (list(lambdaRange=range, prediction=pred, lowerConfidenceCurve=ymin, upperConfidenceCurve=ymax))  
}