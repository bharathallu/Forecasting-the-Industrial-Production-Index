# This script produces ip forecasts

# Clearing
rm(list=ls())

# Accessing necessary libraries
library(astsa)
library(forecast)
library(fGarch)
library(MTS)
library(vars)
library(urca)
library(tseries)
library(cwhmisc)

# Reading in data
data = read.csv(file="/Users/jonandr01/Dropbox/Stat 626 Project/seasonal_full.csv", header=TRUE, sep=",")

# Initializing necessary variables
retail=ts(data[,11], frequency = 12, start=c(1992,1), end=c(2016, 3))
ip=ts(data[,6], frequency = 12, start=c(1992,1), end=c(2016,3))
gdp=ts(data[,12], frequency = 4, start=c(1992,1), end=c(2016, 1))
ipq=aggregate(ip, nfrequency=4)/3

# Setting directory (ONLY NECESSARY IF PRODUCING PDF OUTPUTS)
setwd("/Users/jonandr01/Dropbox/Stat 626 Project/code/ipfcast_pdfs")

# Detrending retail data
pdf('retail.pdf')
plot(retail, main="Retail Sales Data (January 1992 - March 2016)", ylab='U.S. Dollars')
dev.off()

# Taking first difference
dretail=diff(retail)
pdf('dretail.pdf')
plot(dretail, main='Retail Sales Data, First Difference', ylab='U.S. Dollars') # Data appears stationary
dev.off()

# Checking the ACF and PACF of retail
pdf('dretail_acf_pacf.pdf')
par(mfrow=c(2,1))
acf(dretail, main='ACF of Retail Sales Data, First Difference' ) # Close shop!
pacf(dretail, main='PACF of Retail Sales Data, First Difference')
dev.off()

# Detrending ip data
ip=na.omit(ip)
pdf('ip.pdf')
plot(ip, main='Industrial Production, (May 1992 - March 2016)', ylab='Index (2012=100)')
dev.off()

# Taking first difference
dip=diff(ip)
pdf('dip.pdf')
plot(dip, main='Industrial Production, First Difference', ylab=NULL)
dev.off()

# Taking second difference
ddip=diff(dip)
pdf('ddip.pdf')
par(mfrow=c(2,1))
plot(ddip, main='Industrial Production, Second Difference', ylab='Change in Index')
plot(diff(diff(log(ip))), main='Industrial Production, Second Difference Log Transform', ylab = NULL) # We decided on ddip instead of the B2log(ip) transformation
dev.off()

# Checking the ACF and PACF
pdf('ddip_acf_pacf.pdf')
par(mfrow=c(2,1))
acf(ddip, main='ACF of IP, Second Difference')
pacf(ddip, main='PACF of IP, Second Difference')
dev.off()

# Checking scatterplot matrix
pdf('ddip_scatter.pdf')
lag1.plot(ddip,4)
dev.off()

# From the ACF, PACF, and the scatterplot matrix, we determine that our best models are AR(1), AR(2), AR(3), and MA(1)
sarima(ddip,1,0,0, no.constant=TRUE)
sarima(ddip,2,0,0, no.constant=TRUE)
sarima(ddip,3,0,0, no.constant=TRUE)
sarima(ddip,0,0,1, no.constant=TRUE)

# The psi weights do not imply identical models, so we merge to form ARMA models
sarima(ddip,1,0,1, no.constant=TRUE)
pdf('ddip_residuals.pdf')
sarima(ddip,2,0,1, no.constant=TRUE)
dev.off()
sarima(ddip,3,0,1, no.constant=TRUE)

# Holding to perform ljung-box tests
ar1=sarima(ddip,1,0,0, no.constant=TRUE)
ar2=sarima(ddip,2,0,0, no.constant=TRUE)
ar3=sarima(ddip,3,0,0, no.constant=TRUE)
ma1=sarima(ddip,0,0,1, no.constant=TRUE)
arma11=sarima(ddip,1,0,1, no.constant=TRUE) # Second Best
arma21=sarima(ddip,2,0,1, no.constant=TRUE) # Best Model
arma31=sarima(ddip,3,0,1, no.constant=TRUE)

# Checking to see this residuals of an ARMA(2,1) and ARMA(1,1)
Box.test(resid(arma21$fit), type = c("Box-Pierce", "Ljung-Box"), fitdf = 0)
Box.test(resid(arma11$fit), type = c("Box-Pierce", "Ljung-Box"), fitdf = 0)

# Resetting Working Directory
setwd("/Users/jonandr01/Dropbox/Stat 626 Project/code/ipfcast_final_presentation")

######################################### Creating Datasets ###########################################
# Initializing new vars
ddipq = diff(diff(ipq))
ddgdp = diff(diff(gdp))

# New datasets
part1 = ip[1 : as.integer(length(ip)*.9)]
part2 = ip[as.integer(length(ip)*.9) + 1 : length(ip)]
part2 = na.omit(part2)
part1a = ddip[1 : as.integer(length(ddip)*.9)]
part2b = ddip[as.integer(length(ddip)*.9) + 1 : length(ddip)]
part2b  = na.omit(part2b)
ipq90=ddipq[1 : as.integer(length(ddipq)*.9)]
ipq10=ddipq[as.integer(length(ddipq)*.9) + 1 : length(ddipq)]
ipq10 = na.omit(ipq10)
gdp90 = ddgdp[1 : as.integer(length(ddgdp)*.9)]
gdp10=ddipq[as.integer(length(gdp10)*.9) + 1 : length(gdp10)]
gdp10 = na.omit(gdp10)
data2 = cbind(ddipq[86:length(ddipq)],ddgdp[86:length(ddipq)])
IP = ipq90[3:85]
GDP = gdp90[3:85]
data1 = cbind(IP,GDP)
part1ipq = ipq[1 : as.integer(length(ipq)*.9)]
part2ipq = ipq[as.integer(length(ipq)*.9) + 1 : length(ipq)]
part2ipq = na.omit(part2ipq)
part1gdp = gdp[1 : as.integer(length(gdp)*.9)]

######################################### Arima forecast ###########################################
# Fitting model
adf.test(diff(diff(part1)))
fitst=Arima(part1,order=c(2,2,1))

# Forecasting
ip.for=forecast.Arima(fitst,29)
f1 = sarima.for(part1, n.ahead = 29, 2, 2, 1)
seu = f1$pred + f1$se
sel = f1$pred - f1$se
x = cbind(f1$pred, seu, sel,part2)

# Plotting
pdf("arimafit.pdf")
plot(fitst$x,col="red",type="l", ylab = "IP", main="ARIMA Fit")
lines(fitted(fitst),col="blue")
legend(5,100,c("Realization","Forecast"),lty=c(1,1),lwd=c(2,2,2),col=c("blue","red"))
dev.off()
pdf("arimafcast.pdf")
par(mfrow = c(2,1))
plot(x[,1], ylim=c(90,130),lty=2, col = "blue", main="Close Up of ARIMA Forecast", ylab = "IP")
lines(x[,2], col = "red",lty=2)
lines(x[,3], col = "red",lty=2)
lines(x[,4], col = "black")
legend(260,129,c("Realization","Forecast","95% CI"),lty=c(1,2,2),lwd=c(2,2,2),col=c("black","blue","red"), cex = 0.35)
plot(ip.for, main = "ARIMA Forecast") 
dev.off()

# Error
(cumsum((abs(ip.for$mean-part2)/part2))[1])*100 # 0.06676035
(cumsum((abs(ip.for$mean-part2)/part2)/4)[4])*100 # 0.2166906
(cumsum((abs(ip.for$mean-part2)/part2)/10)[10])*100 # 0.4146009

######################################### Garch forecast ##################################################
# Fitting
square11 = resid(arma11$fit)^2
pdf("garchacf.pdf")
par(mfrow=c(2,1))
acf(square11, main = "ACF of ARIMA(1, 2, 1) Squared Residuals") # Implies GARCH(1,0) or GARCH(1,1)
pacf(square11, main = "PACF of ARIMA(1, 2, 1) Squared Residuals")
dev.off()
summary(arma11.g10 <- garchFit(~arma(1,1) + garch(1,1), ddip, include.mean = F)) 
pdf("garchresiduals.pdf")
plot(arma11.g10@residuals, type ="l")
acf(arma11.g10@residuals) 
dev.off()

# Forecasting
pref2 = garchFit(~arma(1,1) + garch(1,1), part1a, include.mean = F)
f2 = predict(pref2, n.ahead = 29)
sdu = f2$meanForecast + f2$meanError
sdl = f2$meanForecast - f2$meanError

# Plotting
pdf("garchfcast.pdf")
par(mfrow=c(2,1))
plot(f2$meanForecast, col = "blue", type = "l", main = "Closeup of GARCH Forecast", ylab = "Double Difference IP")
plot(f2$meanForecast, col = "blue", lty = 2, ylim = c(-1.5,2), type = "l", main = "GARCH Forecast", ylab = "Double Difference IP")
lines(sdu, col = "red", lty = 2)
lines(sdl, col = "red", lty = 2)
lines(part2b, col = "black")
legend("topright",c("Realization","Forecast","95% CI"),lty=c(1,2,2),col=c("black","blue","red"), cex = 0.6, bty = "n", horiz = T)
dev.off()

# Error
((cumsum(abs((f2$meanForecast - part2b)/part2b)))[1])*100 # 13.05125
((cumsum(abs((f2$meanForecast - part2b)/part2b))/4)[4])*100 # 62.29327
((cumsum(abs((f2$meanForecast - part2b)/part2b))/10)[10])*100 # 86.33364

################################################## VAR Forecast ######################################### 
# Fitting
ddipq = na.omit(ddipq)
pdf("var_scatter.pdf")
lag2.plot(ddipq, ddgdp, 3)
dev.off()

VARselect(data1, lag.max=5, type="both")
summary(ipvar <- VAR(data1,p=1,type="both"))
pdf("acf_var.pdf")
acf(resid(ipvar))
dev.off()
serial.test(ipvar, lags.pt=20, type="PT.adjusted")

# Forecasting
(fit.pr=predict(ipvar, n.ahead=10, ci = 0.95))
fanchart(fit.pr)

# Plotting
pdf("varfcast.pdf")
par(mfrow=c(2,1))
plot(fit.pr$fcst$IP[,1],type="l", col = "blue", main = "Close Up of VAR Forecast", ylab = "Double Differenced IP", xlab="Steps Ahead")
plot(fit.pr$fcst$IP[,1], ylim = c(-2,5),type="l", lty=2, col = "blue", main = "VAR Forecast", ylab = "Double Differenced IP", xlab="Steps Ahead")
lines(fit.pr$fcst$IP[,2], type = "l", lty = 2, col = "red")
lines(fit.pr$fcst$IP[,3], type = "l", lty = 2, col = "red")
lines(ipq10, type="l", lty=1, col = "black")
legend(9,5,c("Realization","Forecast","95% CI"),lty=c(1,2,2),lwd=c(1,1,1),col=c("black","blue","red"), cex=0.4)
dev.off()

# Error
(cumsum(abs((fit.pr$fcst$IP[,1] - ipq10)/ipq10)))[1] * 100 # 56.81626
(cumsum(abs((fit.pr$fcst$IP[,1] - ipq10)/ipq10))/4)[4] * 100 # 74.39396
(cumsum(abs((fit.pr$fcst$IP[,1] - ipq10)/ipq10))/10)[10] * 100 # 89.28729

################################################## VECM ##################################################
# Testing for cointegration
IP = part1ipq
GDP = part1gdp
data3 = cbind(IP, GDP)
summary(test <- ca.jo(data3, type="trace", ecdet="none", K=2)) # Our test statistic for r <= 1 suggests first order cointegration

# Fitting
vecm1 = cajorls(test, r=1)
vecmconv=vec2var(test,r=1)
pdf("acf_vecm.pdf")
acf(resid(vecmconv))
dev.off()

# Forecasting 
vecm.pr=predict(vecmconv,n.ahead=10)

# Plotting
pdf("vecmfcast.pdf")
par(mfrow=c(2,1))
plot(vecm.pr$fcst$part1ipq[,1], type = "l", col="blue",main = "Closeup of VECM Forecast", ylab = "Double Differenced IP", xlab="Steps Ahead")
plot(vecm.pr$fcst$part1ipq[,1],type="l",lty=2,ylim=c(90,120), col="blue", main = "VECM Forecast", ylab = "Double Differenced IP", xlab="Steps Ahead")
lines(vecm.pr$fcst$part1ipq[,2],type="l",lty=2,col="red")
lines(vecm.pr$fcst$part1ipq[,3],type="l",lty=2,col="red")
lines(part2ipq, type="l", col="black")
legend(0.75,120,c("Realization","Forecast","95% CI"),lty=c(1,2,2),lwd=c(1,1,1),col=c("black","blue","red"), cex=0.35)
dev.off()

# Error
((cumsum((abs(vecm.pr$fcst$IP[,1]-part2ipq)/part2ipq)))[1])*100 # 0.4005696
((cumsum((abs(vecm.pr$fcst$IP[,1]-part2ipq)/part2ipq)/4))[4])*100 # 1.286753
((cumsum((abs(vecm.pr$fcst$IP[,1]-part2ipq)/part2ipq)/10))[10])*100 # 1.663953


######################################### Exponential Smoothing ######################################### 
# Fitting and forecasting
summary(ip.esm <- ses(part1, h=29, alpha=.9, initial="simple"))

# Plotting
pdf("ewmafcast.pdf")
plot(ip.esm, main = "EWMA Forecast", ylab = "IP")
dev.off()

# Error
(cumsum((abs(ip.esm$mean-part2)/part2))[1])*100 # 0.3352157
(cumsum((abs(ip.esm$mean-part2)/part2)/4)[4])*100 # 0.4277891
(cumsum((abs(ip.esm$mean-part2)/part2)/10)[10])*100 # 1.500849
