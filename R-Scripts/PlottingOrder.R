# This scripts extracts the Monte Carlo quark density data 
# from its data file and computes it variability using the Jackknife procedure to set the error-bars.
# At the same time, it fits a second-order polynomial to the data in the wanted range. 
# Computes the prediction error of the regression model to set the 68% confidence region 
# around the regression curve. The results are visuaized using ggplot2 package

library(ggplot2)
source('./autocorrelationOfMCMCTimeSeries.R')
source('./jackknifeOrder.R')

f<-readLines('./MCDataToAnalyze/OldData/New_N14_L25_M63_Order_30_measur_15000_disc_30.txt')
f<-strsplit(f, "\ |\t")
g<-readLines('./MCDataToAnalyze/OldData/New_N14_L25_M125_Order_30_measur_45000_disc_50.txt')
g<-strsplit(g, "\ |\t")

rangeF<-9000         #20000
stepF<-floor(rangeF/40)
decayF<-8*stepF

rangeG<-21000     #50000
stepG<-floor(rangeG/50)
decayG<-8*stepG

autocorF <- autocorrelationOfMCMCTimeSeries(f, rangeF, stepF, decayF)
tauIntegF<- autocorF$integratedAutoCorr
lambda <- autocorF$lambda

autocorG <- autocorrelationOfMCMCTimeSeries(g, rangeG, stepG, decayG)
tauIntegG<- autocorG$integratedAutoCorr
lambda <- autocorG$lambda

jackknifedF <- jackknifeOrder(lambda, tauIntegF, f)
meanValsF <- jackknifedF$Mean
errorsF <- jackknifedF$JackknifeError

jackknifedG <- jackknifeOrder(lambda, tauIntegG, g)
meanValsG <- jackknifedG$Mean
errorsG <- jackknifedG$JackknifeError
dev.off()
dtf_F <- data.frame(x=lambda, y=meanValsF, df=factor("K=61",levels=c("K=61","K=123")))
dtf_G <- data.frame(x=lambda, y=meanValsG, df=factor("K=123",levels=c("K=61","K=123")))
dft <- rbind(dtf_F,dtf_G)
x11()
qplot(x, y,color=df,data=dft, color=df)+ geom_point(shape=18,  color=df, size=1.5)+
  theme(legend.position="right")+
  theme(panel.background = element_rect(color='blue', fill ='lightyellow'))+
  geom_errorbar(data=dtf_F, aes(x=lambda, ymin=meanValsF-errorsF, ymax=meanValsF+errorsF), color='red', width=0.015, size=0.2)+ 
  geom_errorbar(data=dtf_G, aes(x=lambda, ymin=meanValsG-errorsG, ymax=meanValsG+errorsG), color='darkgreen', width=0.015, size=0.2)+
  scale_colour_manual(values=c('red', 'darkgreen'))+
  theme(legend.title=element_blank())+
  labs(x=expression(lambda), y="n", size=1.0)+ theme(axis.title.x = element_text(color='purple', size = 20, angle = 00), axis.title.y = element_text(color='purple', size = 20, angle = 00, family='serif', vjust=0.55))+
  xlim(0.5, 1.6)+ylim(0, 0.65)+
  annotate("text", label="N=14, L=25", x=0.65, y=0.6, family="serif", fontface="plain", color='blue')
  ggsave("./KsOverlayedN14L25Color.pdf")
