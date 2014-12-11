# This function computes a maximum likelihood estimate for the regression coefficients 
# for a polynomial (specialized to a quadratic) regression model.

ChiSquareFit<-function(x, y, sigma){
  
  s0 <- sum(1/sigma^2)
  s1 <- sum(x/sigma^2)
  s2 <- sum(x^2/sigma^2)
  s3 <- sum(x^3/sigma^2)
  s4 <- sum(x^4/sigma^2)
  
  b0 <- sum(y/sigma^2)
  b1 <- sum(x*y/sigma^2)
  b2 <- sum(x^2*y/sigma^2)

 Delta <- -s2^3+s2*s4*s0+2*s2*s3*s1-s3^2*s0-s4*s1^2
 
 numerA <- -s1^2*b2+s1*s3*b0+s1*s2*b1-s0*s3*b1+s0*b2*s2-b0*s2^2
 numerB <- -(-s1*s2*b2+s1*s4*b0+s3*s0*b2-s3*s2*b0+s2^2*b1-s0*s4*b1)
 numerC <- s2*s3*b1-b2*s2^2+s2*s4*b0+s3*s1*b2-s3^2*b0-s4*s1*b1
 
 a <- numerA/Delta
 b <- numerB/Delta
 c <- numerC/Delta 
 
 return(c(a, b, c))
}