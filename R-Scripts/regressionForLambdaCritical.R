# It fits a linear model to the susceptibility maximum location data, estimates uncertainty of the regression coefficients and visualizes the data using ggplot2.



library(ggplot2)
  source('./estimateErrorLine.R')
  
  # Plotting the extrapolating line for the critical lambda (lambda ~ inverse size)
  ep <- 10^(-2)
  inverseSize <- c(1/6, 1/7, 1/7.5, 1/8, 1/9)
  
  lambdaCritUpper <- c(1.066179, 1.096791,  1.123031, 1.1449, 1.164144 )
  lambdaCritLower <- c(1.066179-ep, 1.096791-ep, 1.119533, 1.141401, 1.16152)
  
  fitLine <- estimateErrorLine(inverseSize, lambdaCritUpper, lambdaCritLower, 200000)
  
  regressors <- fitLine$meanCoeffs
  regressionError <- fitLine$stdErr
  
  print("Printing regression error: ")
  print(regressionError)
  
  print(regressors)
  
  yMean <- (lambdaCritUpper+lambdaCritLower)/2
  
  stdEr <- (lambdaCritUpper-lambdaCritLower)/2
  
  dtf_data <- data.frame(inverseSize, yMean)
  
  inverseSizeRange <- seq(0, 0.17, 0.17/299)
  pred <- regressors[1]*inverseSizeRange+regressors[2]
  dtf_reg <- data.frame(inverseSizeRange, pred)
  
  #plot the line and the error bars 
  
  pdf('./ContinuumExtrapolationLambdaCritical.pdf')
  
  ggplot()+ 
    geom_line(data=dtf_reg, aes(y=pred, x=inverseSizeRange, linetype="Fit"), color='blue', size=0.05)+
    geom_point(data=dtf_data, aes(y=yMean, x=inverseSize), shape=23, fill='red', color='blue', size=1.8)+
    geom_errorbar(data=dtf_data, aes(x=inverseSize, ymin=yMean-stdEr, ymax=yMean+stdEr), color='red', width=0.006, size=0.1)+ 
    theme(panel.background = element_rect(color='blue', fill = 'lightyellow'))+ 
    theme(legend.title=element_blank())+
    labs(x="1/L", y=expression(lambda), size=2)+ theme(axis.title.x = element_text(color='purple', size = 20, angle = 00), axis.title.y = element_text(color='purple', size = 20, angle = 00, family='serif', vjust=0.55))+
    theme(axis.text.x=element_text(angle=00, size=15, vjust=0.5), axis.text.y=element_text(angle=00, size=15, vjust=0.5))
  dev.off()
