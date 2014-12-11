# It fits a linear model to the susceptibility maximum value data, estimates uncertainty of the regression coefficients and visualizes the data using ggplot2.  

library(ggplot2)
source('./estimateErrorLine.R')

# Plotting the extrapolating line for the critical lambda (lambda ~ inverse size)

inverseSize <- c(1/6, 1/7, 1/7.5, 1/8, 1/9)

XLMaxUpper <- c(0.4277315, 0.5020015, 0.5562661, 0.6716513, 0.847757)
XLMaxLower <- c(0.4111787, 0.4393566, 0.5046158, 0.6365836, 0.7106627)

fitLine <- estimateErrorLine(inverseSize, XLMaxUpper, XLMaxLower, 200000)

regressors <- fitLine$meanCoeffs
regressionError <- fitLine$stdErr

print("Printing regression error: ")
print(regressionError)

print(regressors)

yMean <- (XLMaxUpper+XLMaxLower)/2

stdEr <- (XLMaxUpper-XLMaxLower)/2

dtf_data <- data.frame(inverseSize, yMean)

inverseSizeRange <- seq(0, 0.17, 0.17/299)
pred <- regressors[1]*inverseSizeRange+regressors[2]
dtf_reg <- data.frame(inverseSizeRange, pred)

#plot the line and the error bars 

pdf('./FitForXLMax.pdf')

ggplot()+ 
  geom_line(data=dtf_reg, aes(y=pred, x=inverseSizeRange), shape=19, fill='blue', color='blue', size=0.1)+
  geom_point(data=dtf_data, aes(y=yMean, x=inverseSize), shape=23, fill='red', color='blue', size=1.8)+
  geom_errorbar(data=dtf_data, aes(x=inverseSize, ymin=yMean-stdEr, ymax=yMean+stdEr), color='red', width=0.01, size=0.5)+ 
  theme(panel.background = element_rect(color='blue', fill = 'white'))+ 
  labs(x="1/L", y=expression(lambda), size=2)+ theme(axis.title.x = element_text(color='purple', size = 20, angle = 00), axis.title.y = element_text(color='purple', size = 20, angle = 00, family='serif', vjust=0.55))+
  theme(axis.text.x=element_text(angle=00, size=15, vjust=0.5), axis.text.y=element_text(angle=00, size=15, vjust=0.5))
dev.off()
