#LIBRARIES
#####
library(fGarch)
library(tseries)
library(readxl)
library(rugarch)
library(lmtest)
library(moments)
library(pracma)
library(fBasics)
library(FinTS)
library(zoo)
library(TTR)
library(GeneralizedHyperbolic)
library(fBasics)
library(data.table)
library(xts)
library(rumidas)
library(GAS)
library(MCS)
#####
#import the data:
#####
energy_data<- read_excel("~/R/energy_covid_data.xlsx", sheet = "daily data", skip = 4)
energy_data<-na.omit(energy_data)
energy_data<-energy_data[-36,]
covid_data<- read_excel("~/R/energy_covid_data.xlsx", sheet = "COVID deaths USA", skip = 4)
covid_data<-na.omit(covid_data)
#create the dates:
date_R<-strptime(energy_data$...1, "%Y-%m-%d",tz="GMT")
date_R2<-as.Date(date_R)
date_VaR3<-date_R2[217:316]
date_VaR3<-na.omit(date_VaR3)
length(date_VaR3)
date_covid<-strptime(covid_data$...1, "%Y-%m-%d",tz="GMT")
#Create variable for the Covid deaths:
dUSA<-covid_data$`Deaths USA`
dUSA_xts<-as.xts(dUSA,date_covid)
log_covid_xts<-diff(log(dUSA_xts))#we take the first log differences to take away the trend component
#####
#Create WTI Crude Oil time series:
#####
Price_WTI<-energy_data$WTI
Price_WTI_xts<-as.xts(price_WTI,date_R)
Price_WTI_xts<-na.omit(Price_WTI_xts)
WTI_ts<-ts(Price_WTI_xts)
basicStats(Price_WTI_xts)
acf(WTI_ts,lag.max=315, main="WTI PRICE TIME SERIES ACF")
pacf(WTI_ts,lag.max=315, main="WTI PRICE TIME SERIES PACF")
plot.xts(Price_WTI_xts,main="WTI PRICE TIME SERIES")
plot(density(WTI_ts),main="WTI PRICE EMPIRICAL DENSITY")
#####
#LOG RETURNS:
#####
logR_WTI<-diff(log(price_WTI))*100
logR_WTI_xts<-diff(log(Price_WTI_xts))*100
logR_WTI_xts <- na.omit(logR_WTI_xts)
logR_WTI_ts<-ts(logR_WTI_xts)
#DESCRIPTIVE STATISTICS, ACF and PACF, plots of returns: 
basicStats(logR_WTI_xts)
acf(logR_WTI_ts,lag.max=316,main="ACF RETURNS WTI TIME SERIES")
pacf(logR_WTI_ts,lag.max=316,main="PACF RETURNS WTI TIME SERIES")
plot.ts(logR_WTI_ts,main="RETURNS WTI TIME SERIES")
plot(density(logR_WTI_ts),main="WTI RETURNS EMPIRICAL DENSITY")
hist(logR_WTI_xts ,nclass=100, freq=F,col="green",main="HISTOGRAM OF DAILY RETURNS", ylab="Density",xlab="Daily returns") 
#absolute returns:
acf(abs(logR_WTI_ts),lag.max=316,main="ACF ABSOLUTE VALUE RETURNS WTI TIME SERIES")
pacf(abs(logR_WTI_ts),lag.max=316,main="PACF ABSOLUTE VALUE RETURNS WTI TIME SERIES")
#squared returns:
acf((logR_WTI_ts)^2,lag.max=316,main="ACF SQUARED RETURNS WTI TIME SERIES")
pacf((logR_WTI_ts)^2,lag.max=316,main="PACF SQUARED RETURNS WTI TIME SERIES")
#####
#DIAGNOSTIC TESTS:
#####
#DICKEY FULLER TEST (stationarity): 
adf.test(logR_WTI_ts, alternative="stationary")#the series is stationary(we reject H0=not stationary)
adf.test(logR_WTI_ts, alternative="explosive")
#Normality test:
qqnorm(logR_WTI_ts)
qqline(logR_WTI_ts,col="red")#from the graph we hp a student w/ big tails
jarque.bera.test(logR_WTI_ts)#reject normality
shapiro.test(logR_WTI_ts)#reject normality
agostino.test(logR_WTI_ts)#reject normality
anscombe.test(logR_WTI_ts)#reject normality
#other tests:
Box.test((logR_WTI_ts), lag =6, type = "Ljung-Box") #optimal lag=ln(n)=ln(315)=5.7 #there is correlation 
Box.test((logR_WTI_ts), lag = 6, type = "Box-Pierce")
Box.test(((logR_WTI_ts)^2), lag = 6, type = "Box-Pierce")
kpss.test(logR_WTI_ts)
kpss.test(WTI_ts)
bds.test(logR_WTI_ts)#not iid
ArchTest(logR_WTI_ts)# there is an autoregressive heteroskedastic component
#####
#Model GARCH(1,1)using various error distributons:
#####
#NORMAL
#####
#fit:
Gnorm_spec_WTI=ugarchspec(variance.model=list(garchOrder=c(1,1)), mean.model=list(armaOrder=c(0,0)),distribution.model = "norm")
Gnorm_fit_WTI<-ugarchfit(spec=Gnorm_spec_WTI, data=logR_WTI_xts,out.sample = 100)
Gnorm_fit_WTI
show(coef(Gnorm_fit_WTI))
show(confint(Gnorm_fit_WTI))
#VOLATILTIY:
plot(Gnorm_fit_WTI, which = "all")
sigma_WTI_xts<-as.xts(sigma(Gnorm_fit_WTI))
plot.xts(sigma_WTI_xts,ylab="sigma(t)", col="blue", main="ESTIMATED VOLATILITY TIME SRIES"))NORM
#residual analysis: 
resid_Gnorm_WTI<-Gnorm_fit_WTI@fit$residuals
acf(resid_Gnorm_WTI, lag.max=220, main="ACF OF RESIDUALS")#ts of the residuals(errors) beahves like a WN
pacf(resid_Gnorm_WTI, lag.max=220, main="PACF OF Norm RESIDUALS")
qqnorm(resid_Gnorm_WTI)
grid()
qqline(resid_Gnorm_WTI,lwd=2,col="red")
sq_resid_Gnorm_WTI <-resid_Gnorm_WTI^2
acf(sq_resid_Gnorm_WTI, lag.max = 250, main="ACF OF SQUARE RESIDUALS") #ts of the squared residuals is autocorrelated as expected
#FORECASTING:
forecast_Gnorm_fit_WTI <-ugarchforecast(Gnorm_fit_WTI,data=logR_WTI_xts,n.ahead = 100)
plot(forecast_Gnorm_fit_WTI, which = 3)#forecasted volatility graph (Sigma)
#####
#Value at risk for GARCH (norm)
#####
rollingforecast_Gnorm_WTI= ugarchroll(Gnorm_spec_WTI, data = logR_WTI_xts,n.ahead = 1, forecast.length = 100,refit.every =10, refit.window = c("recursive"), solver = "hybrid", calculate.VaR = TRUE, VaR.alpha = c(0.01, 0.05), cluster = NULL, keep.coef = TRUE)
show(rollingforecast_Gnorm_WTI)#we predict from 12 jan 2021 to 4 june 2022
rollingforecast_Gnorm_WTI@forecast$VaR
plot(rollingforecast_Gnorm_WTI,which="all")
#Backtesting
report(rollingforecast_Gnorm_WTI, type = "VaR", VaR.alpha = 0.01, conf.level = 0.99)
report(rollingforecast_Gnorm_WTI, type = "VaR", VaR.alpha = 0.05, conf.level = 0.95)#both test should have Reject Null: NO, if kubiec no e cristof yes is ok anyways
report(rollingforecast_Gnorm_WTI, type = "fpm") #con mean squared error, mean absolute error e directional accuracy
plot.ts(rollingforecast_Gnorm_WTI@forecast[["VaR"]])
show(basicStats(rollingforecast_Gnorm_WTI@forecast[["VaR"]]))

WTIVaR1.norm=rollingforecast_Gnorm_WTI@forecast$VaR[,1]#alpha=1%
WTIVaR5.norm=rollingforecast_Gnorm_WTI@forecast$VaR[,2]#alpha=5%
WTIReal.norm=rollingforecast_Gnorm_WTI@forecast$VaR[,3]#realized

WTIVaR5.norm.xts<-as.xts(WTIVaR5.norm, tail(date_VaR3,100))
WTIVaR1.norm.xts<-as.xts(WTIVaR1.norm, tail(date_VaR3,100))
WTIReal.norm.xts<-as.xts(WTIReal.norm, tail(date_VaR3,100))
plot(WTIReal.norm.xts, type = "l", pch = 16, cex = 0.8,  col = "blue",
     ylab = "Returns", main = "Realized vs VaR Forecast", xaxt = "n", ylim=c(-12,6))
lines(WTIVaR1.norm.xts, col = "purple")#1%
lines(WTIVaR5.norm.xts, col = "red")#5%
legend('topright', c("returns", "99% VaR", "95% VaR"),
       lty=1, col=c("blue", "purple","red"), bty='n', cex=.75)

WTItest1.norm=VaRTest(alpha = 0.01, WTIReal.norm, WTIVaR1.norm, conf.level =  0.99)
WTItest5.norm=VaRTest(alpha = 0.05, WTIReal.norm, WTIVaR5.norm, conf.level =  0.95)

logR_WTI_last100_xts<-as.xts(logR_WTI_xts["2021-01-12/2021-06-04"],date_VaR3)
BackTestVaR1_WTI_GARCH_norm<-BacktestVaR(data =logR_WTI_last100_ts,VaR=WTIVaralpha1.norm,alpha = 0.01,Lags =6)
BackTestVaR5_WTI_GARCH_norm<-BacktestVaR(data =logR_WTI_last100_ts,VaR=WTIVaralpha2.norm,alpha = 0.05,Lags =6)

plot(tail(logR_WTI_xts,100), type = "l", pch = 16, cex = 0.8,  col = "blue",
     ylab = "Returns", main = "Comparison between returns and Var 95% and 99%", xaxt = "n",ylim=c(-12,6))
lines(tail(WTIVaR1.norm.xts,100), col = "purple")
lines(tail(WTIVaR5.norm.xts,100), col = "red")
legend('topright', c("sGARCH(1,1)"),
       lty=1, col=c("purple"), bty='n', cex=.75)

preds_WTI<-as.data.frame(rollingforecast_Gnorm_WTI)
garchvol_WTI<-xts(preds_WTI$Sigma, order.by=as.Date(rownames(preds_WTI)))
plot.xts(garchvol_WTI)
#####
#SKEW NORMAL
#####
#fit:
Gsnorm_spec_WTI=ugarchspec(variance.model=list(garchOrder=c(1,1)), mean.model=list(armaOrder=c(0,0)),distribution.model = "snorm")
Gsnorm_fit_WTI<-ugarchfit(spec=Gsnorm_spec_WTI, data=logR_WTI_xts,out.sample = 100)
Gsnorm_fit_WTI
#VOLATILTIY:
plot(Gsnorm_fit_WTI, which = "all")
sigma_WTI_snorm_xts<-as.xts(sigma(Gsnorm_fit_WTI))
plot.xts(sigma_WTI_snorm_xts,ylab="sigma(t)", col="blue", main="ESTIMATED VOLATILITY TIME SRIES")#volatilit? stimata dei rendimenti con il GARCH(1,1)snorm
#error analysis: 
resid_Gsnorm_WTI<-Gsnorm_fit_WTI@fit$residuals
acf(resid_Gsnorm_WTI, lag.max=220, main="ACF OF RESIDUALS")#ts of the residuals(errors) is a WN
pacf(resid_Gsnorm_WTI, lag.max=220, main="PACF OF RESIDUALS")
qqnorm(resid_Gsnorm_WTI)
grid()
qqline(resid_Gsnorm_WTI,lwd=2,col="red")
sq_resid_Gsnorm_WTI <-resid_Gsnorm_WTI^2
acf(sq_resid_Gsnorm_WTI, lag.max = 250, main="ACF OF SQUARE RESIDUALS") #ts of the squared residuals is autocorrelated, then the ts has a conditional heterscedasticity
pacf(sq_resid_Gsnorm_WTI, lag.max = 250,main="PACF OF SQUARE RESIDUALS")
#FORECASTING:
forecast_Gsnorm_fit_WTI <-ugarchforecast(Gsnorm_fit_WTI,data=logR_WTI_xts,n.ahead = 100)
forecast_Gsnorm_fit_WTI
plot(forecast_Gsnorm_fit_WTI, which =4 ) #grafico dei forecast della nostra serie storica 
plot(forecast_Gsnorm_fit_WTI, which = 3)#grafico dei forecast della volatilità (Sigma)
#####
#Value at risk for GARCH (snorm)
#####
rollingforecast_Gsnorm_WTI= ugarchroll(Gsnorm_spec_WTI, data = logR_WTI_xts,n.ahead = 1, forecast.length = 100,refit.every =10, refit.window = c("recursive"), solver = "hybrid", calculate.VaR = TRUE, VaR.alpha = c(0.01, 0.05), cluster = NULL, keep.coef = TRUE)
show(rollingforecast_Gsnorm_WTI)#we predict from 12 jan 2021 to 4 june 2022
roll.WTI.snorm=resume(rollingforecast_Gsnorm_WTI)
rollingforecast_Gsnorm_WTI@forecast$VaR#la serie storica stimata (all'interno di forcrolling) dei VaR
plot(rollingforecast_Gsnorm_WTI,which="all")
plot(rollingforecast_Gsnorm_WTI,which=4)
coef(rollingforecast_Gsnorm_WTI)
#Backtesting
report(rollingforecast_Gsnorm_WTI, type = "VaR", VaR.alpha = 0.01, conf.level = 0.99)
report(rollingforecast_Gsnorm_WTI, type = "VaR", VaR.alpha = 0.05, conf.level = 0.95)#both test should have Reject Null: NO, if kubiec no e cristof yes is ok anyways
report(rollingforecast_Gsnorm_WTI, type = "fpm") #con mean squared error, mean absolute error e directional accuracy
plot.ts(rollingforecast_Gsnorm_WTI@forecast[["VaR"]])
analisys_VaR_WTI.snorm=basicStats(rollingforecast_Gsnorm_WTI@forecast[["VaR"]])
show(analisys_VaR_WTI.snorm)

WTIVaralpha1.snorm=rollingforecast_Gsnorm_WTI@forecast$VaR[,1]#1%
WTIVaralpha2.snorm=rollingforecast_Gsnorm_WTI@forecast$VaR[,2]#5%
WTIVaralpha3.snorm=rollingforecast_Gsnorm_WTI@forecast$VaR[,3]#realized
WTItest1.snorm=VaRTest(alpha = 0.01, WTIVaralpha3.snorm, WTIVaralpha1.snorm, conf.level =  0.99)
WTItest1.snorm$actual.exceed
WTItest5.snorm=VaRTest(alpha = 0.05, WTIVaralpha3.snorm, WTIVaralpha2.snorm, conf.level =  0.95)
WTItest5.snorm$expected.exceed

BackTestVaR1_WTI_GARCH_snorm<-BacktestVaR(data =logR_WTI_last100_ts,VaR=WTIVaralpha1.snorm,alpha = 0.01,Lags =6)
BackTestVaR1_WTI_GARCH_snorm
BackTestVaR5_WTI_GARCH_snorm<-BacktestVaR(data =logR_WTI_last100_ts,VaR=WTIVaralpha2.snorm,alpha = 0.05,Lags =6)
BackTestVaR5_WTI_GARCH_snorm

plot(logR_WTI_last100_ts, type = "l", pch = 16, cex = 0.8,  col = gray(0.2, 0.5),
     ylab = "Returns", main = "95% VaR Forecasting", xaxt = "n", ylim=c(-12,6))
lines(WTIVaralpha1.snorm, col = "purple")
legend('topright', c("sGARCH(1,1)"),
       lty=1, col=c("purple"), bty='n', cex=.75)

logR_WTI_last1002_xts<-as.xts(logR_WTI_xts["2021-01-12/2021-06-04"],date_VaR3)
WTIVaralpha1.snorm.xts2<-as.xts(WTIVaralpha1.snorm,date_VaR3)
WTIVaralpha2.snorm.xts2<-as.xts(WTIVaralpha2.snorm,date_VaR3)
plot(logR_WTI_last1002_xts, type = "l", pch = 16, cex = 0.8,  col = "blue",
     ylab = "Returns", main = "Comparison between returns and Var 95% and 99%", xaxt = "n",ylim=c(-12,6))
lines(WTIVaralpha1.snorm.xts2, col = "purple")
lines(WTIVaralpha2.snorm.xts2, col = "red")

plot(logR_WTI_last100_ts, type = "l", pch = 16, cex = 0.8,  col = gray(0.2, 0.5),
     ylab = "Returns", main = "95% VaR Forecasting", xaxt = "n")
lines(WTIVaralpha2.snorm, col = "purple")
legend('topright', c("sGARCH(1,1)"),
       lty=1, col=c("purple"), bty='n', cex=.75)

preds_WTI.snorm<-as.data.frame(rollingforecast_Gsnorm_WTI)
preds_WTI.snorm$Mu
preds_WTI.snorm$Sigma
preds_WTI.snorm$Realized
garchvol_WTI.snorm<-xts(preds_WTI$Sigma, order.by=as.Date(rownames(preds_WTI.snorm)))
plot.xts(garchvol_WTI.snorm)
#prediction error for the mean: 
e_WTI.snorm<-preds_WTI.snorm$Realized-preds_WTI$Mu
mean(e_WTI.snorm^2)
#prediction error for the variance: 
d_WTI.snorm<-(e_WTI.snorm^2-(preds_WTI.snorm$Sigma)^2)
mean(d_WTI.snorm^2)
#####
#STUDENT
#####
#fit:
Gstd_spec_WTI=ugarchspec(variance.model=list(garchOrder=c(1,1)), mean.model=list(armaOrder=c(0,0)),distribution.model = "std")
Gstd_fit_WTI<-ugarchfit(spec=Gstd_spec_WTI, data=logR_WTI_xts,out.sample = 100)
Gstd_fit_WTI
WTI.coeff.std = coef(Gstd_fit_WTI)
show(WTI.coeff.std)
WTI.intconf.std = confint(Gstd_fit_WTI)
show(WTI.intconf.std)
#VOLATILTIY:
sigma_WTI_std_xts<-as.xts(sigma(Gstd_fit_WTI))
plot.xts(sigma_WTI_std_xts,ylab="sigma(t)", col="blue", main="ESTIMATED VOLATILITY TIME SRIES")#volatilit? stimata dei rendimenti con il GARCH(1,1)std
#residual analysis: 
resid_Gstd_WTI<-Gstd_fit_WTI@fit$residuals
acf(resid_Gstd_WTI, lag.max=220, main="ACF OF RESIDUALS")#ts of the residuals(errors) is a WN
pacf(resid_Gstd_WTI, lag.max=220, main="PACF OF RESIDUALS")
qqnorm(resid_Gstd_WTI)
grid()
qqline(resid_Gstd_WTI,lwd=2,col="red")
sq_resid_Gstd_WTI <-resid_Gstd_WTI^2
acf(sq_resid_Gstd_WTI, lag.max = 250, main="ACF OF SQUARE RESIDUALS") #ts of the squared residuals is autocorrelated, then the ts has a conditional heterscedasticity
pacf(sq_resid_Gstd_WTI, lag.max = 250,main="PACF OF SQUARE RESIDUALS")
#FORECASTING:
forecast_Gstd_fit_WTI <-ugarchforecast(Gstd_fit_WTI,data=logR_WTI_xts,n.ahead = 100)
forecast_Gstd_fit_WTI
plot(forecast_Gstd_fit_WTI, which = 1) #grafico dei forecast della nostra serie storica 
plot(forecast_Gstd_fit_WTI, which = 3)#grafico dei forecast della volatilità (Sigma)
#####
#Value at risk for GARCH (std)
#####
rollingforecast_Gstd_WTI= ugarchroll(Gstd_spec_WTI, data = logR_WTI_xts,n.ahead = 1, forecast.length = 100,refit.every =10, refit.window = c("recursive"), solver = "hybrid", calculate.VaR = TRUE, VaR.alpha = c(0.01, 0.05), cluster = NULL, keep.coef = TRUE)
show(rollingforecast_Gstd_WTI)#we predict from 12 jan 2021 to 4 june 2022
roll.WTI.std=resume(rollingforecast_Gstd_WTI)
rollingforecast_Gstd_WTI@forecast$VaR#la serie storica stimata (all'interno di forcrolling) dei VaR
plot(rollingforecast_Gstd_WTI,which="all")
plot(rollingforecast_Gstd_WTI,which=4)
coef(rollingforecast_Gstd_WTI)
#Backtesting
report(rollingforecast_Gstd_WTI, type = "VaR", VaR.alpha = 0.01, conf.level = 0.99)
report(rollingforecast_Gstd_WTI, type = "VaR", VaR.alpha = 0.05, conf.level = 0.95)#both test should have Reject Null: NO, if kubiec no e cristof yes is ok anyways
report(rollingforecast_Gstd_WTI, type = "fpm") #con mean squared error, mean absolute error e directional accuracy
plot.ts(rollingforecast_Gstd_WTI@forecast[["VaR"]])
analisys_VaR_WTI.std=basicStats(rollingforecast_Gstd_WTI@forecast[["VaR"]])
show(analisys_VaR_WTI.std)

WTIVaralpha1.std=rollingforecast_Gstd_WTI@forecast$VaR[,1]
WTIVaralpha2.std=rollingforecast_Gstd_WTI@forecast$VaR[,2]
WTIVaralpha3.std=rollingforecast_Gstd_WTI@forecast$VaR[,3]
WTItest1.std=VaRTest(alpha = 0.01, WTIVaralpha3.std, WTIVaralpha1.std, conf.level =  0.99)
WTItest1.std$expected.exceed
WTItest5.std=VaRTest(alpha = 0.05, WTIVaralpha3.std, WTIVaralpha2.std, conf.level =  0.95)
WTItest5.std$expected.exceed

BackTestVaR1_WTI_GARCH_std<-BacktestVaR(data =logR_WTI_last100_ts,VaR=WTIVaralpha1.std,alpha = 0.01,Lags =6)
BackTestVaR1_WTI_GARCH_std
BackTestVaR5_WTI_GARCH_std<-BacktestVaR(data =logR_WTI_last100_ts,VaR=WTIVaralpha2.std,alpha = 0.05,Lags =6)
BackTestVaR5_WTI_GARCH_std

plot(logR_WTI_last100_ts, type = "l", pch = 16, cex = 0.8,  col = gray(0.2, 0.5),
     ylab = "Returns", main = "95% VaR Forecasting", xaxt = "n", ylim=c(-12,6))
lines(WTIVaralpha1.std, col = "purple")
legend('topright', c("sGARCH(1,1)"),
       lty=1, col=c("purple"), bty='n', cex=.75)

logR_WTI_last1002_xts<-as.xts(logR_WTI_xts["2021-01-12/2021-06-04"],date_VaR3)
WTIVaralpha1.std.xts2<-as.xts(WTIVaralpha1.std,date_VaR3)
WTIVaralpha2.std.xts2<-as.xts(WTIVaralpha2.std,date_VaR3)
plot(logR_WTI_last1002_xts, type = "l", pch = 16, cex = 0.8,  col = "blue",
     ylab = "Returns", main = "Comparison between returns and Var 95% and 99%", xaxt = "n",ylim=c(-12,6))
lines(WTIVaralpha1.std.xts2, col = "purple")
lines(WTIVaralpha2.std.xts2, col = "red")


plot(logR_WTI_last100_ts, type = "l", pch = 16, cex = 0.8,  col = gray(0.2, 0.5),
     ylab = "Returns", main = "95% VaR Forecasting", xaxt = "n")
lines(WTIVaralpha2.std, col = "purple")
legend('topright', c("sGARCH(1,1)"),
       lty=1, col=c("purple"), bty='n', cex=.75)

preds_WTI.std<-as.data.frame(rollingforecast_Gstd_WTI)
preds_WTI.std$Mu
preds_WTI.std$Sigma
preds_WTI.std$Realized
garchvol_WTI.std<-xts(preds_WTI$Sigma, order.by=as.Date(rownames(preds_WTI.std)))
plot.xts(garchvol_WTI.std)
#prediction error for the mean: 
e_WTI.std<-preds_WTI.std$Realized-preds_WTI$Mu
mean(e_WTI.std^2)
#prediction error for the variance: 
d_WTI.std<-(e_WTI.std^2-(preds_WTI.std$Sigma)^2)
mean(d_WTI.std^2)
#####
#SKEW STUDENT
#####
#fit:
Gsstd_spec_WTI=ugarchspec(variance.model=list(garchOrder=c(1,1)), mean.model=list(armaOrder=c(0,0)),distribution.model = "sstd")
Gsstd_fit_WTI<-ugarchfit(spec=Gsstd_spec_WTI, data=logR_WTI_xts,out.sample = 100)
Gsstd_fit_WTI
WTI.coeff.sstd = coef(Gsstd_fit_WTI)
show(WTI.coeff.sstd)
WTI.intconf.sstd = confint(Gsstd_fit_WTI)
show(WTI.intconf.sstd)
#optimal parameters,t value:i parametri stimati sono tutti significativi oltre ogni ragionevole intervallo di confidenza
#VOLATILTIY:
plot(Gsstd_fit_WTI, which = "all")
sigma_WTI_sstd_xts<-as.xts(sigma(Gsstd_fit_WTI))
plot.xts(sigma_WTI_sstd_xts,ylab="sigma(t)", col="blue", main="ESTIMATED VOLATILITY TIME SRIES")#volatilit? stimata dei rendimenti con il GARCH(1,1)sstd
#residual analysis: 
resid_Gsstd_WTI<-Gsstd_fit_WTI@fit$residuals
acf(resid_Gsstd_WTI, lag.max=220, main="ACF OF RESIDUALS")#ts of the residuals(errors) is a WN
pacf(resid_Gsstd_WTI, lag.max=220, main="PACF OF std RESIDUALS")
qqnorm(resid_Gsstd_WTI)
grid()
qqline(resid_Gsstd_WTI,lwd=2,col="red")
sq_resid_Gsstd_WTI <-resid_Gsstd_WTI^2
acf(sq_resid_Gsstd_WTI, lag.max = 250, main="ACF OF SQUARE RESIDUALS") #ts of the squared residuals is autocorrelated, then the ts has a conditional heterscedasticity
pacf(sq_resid_Gsstd_WTI, lag.max = 250,main="PACF OF SQUARE RESIDUALS")
#FORECASTING:
forecast_Gsstd_fit_WTI <-ugarchforecast(Gsstd_fit_WTI,data=logR_WTI_xts,n.ahead = 100)
forecast_Gsstd_fit_WTI
plot(forecast_Gsstd_fit_WTI, which = 1) #grafico dei forecast della nostra serie storica 
plot(forecast_Gsstd_fit_WTI, which = 3)#grafico dei forecast della volatilità (Sigma)
#Value at risk for GARCH (sstd)
#####
rollingforecast_Gsstd_WTI= ugarchroll(Gsstd_spec_WTI, data = logR_WTI_xts,n.ahead = 1, forecast.length = 100,refit.every =10, refit.window = c("recursive"), solver = "hybrid", calculate.VaR = TRUE, VaR.alpha = c(0.01, 0.05), cluster = NULL, keep.coef = TRUE)
show(rollingforecast_Gsstd_WTI)#we predict from 12 jan 2021 to 4 june 2022
roll.WTI.sstd=resume(rollingforecast_Gsstd_WTI)
rollingforecast_Gsstd_WTI@forecast$VaR#la serie storica stimata (all'interno di forcrolling) dei VaR
plot(rollingforecast_Gsstd_WTI,which="all")
plot(rollingforecast_Gsstd_WTI,which=4)
coef(rollingforecast_Gsstd_WTI)
#Backtesting
report(rollingforecast_Gsstd_WTI, type = "VaR", VaR.alpha = 0.01, conf.level = 0.99)
report(rollingforecast_Gsstd_WTI, type = "VaR", VaR.alpha = 0.05, conf.level = 0.95)#both test should have Reject Null: NO, if kubiec no e cristof yes is ok anyways
report(rollingforecast_Gsstd_WTI, type = "fpm") #con mean squared error, mean absolute error e directional accuracy
plot.ts(rollingforecast_Gsstd_WTI@forecast[["VaR"]])
analisys_VaR_WTI.sstd=basicStats(rollingforecast_Gsstd_WTI@forecast[["VaR"]])
show(analisys_VaR_WTI.sstd)

WTIVaralpha1.sstd=rollingforecast_Gsstd_WTI@forecast$VaR[,1]
WTIVaralpha2.sstd=rollingforecast_Gsstd_WTI@forecast$VaR[,2]
WTIVaralpha3.sstd=rollingforecast_Gsstd_WTI@forecast$VaR[,3]
WTItest1.sstd=VaRTest(alpha = 0.01, WTIVaralpha3.sstd, WTIVaralpha1.sstd, conf.level =  0.99)
WTItest1.sstd$expected.exceed
WTItest5.sstd=VaRTest(alpha = 0.05, WTIVaralpha3.sstd, WTIVaralpha2.sstd, conf.level =  0.95)
WTItest5.sstd$expected.exceed

logR_WTI_last1002_xts<-as.xts(logR_WTI_xts["2021-01-12/2021-06-04"],date_VaR3)
WTIVaralpha1.sstd.xts2<-as.xts(WTIVaralpha1.sstd,date_VaR3)
WTIVaralpha2.sstd.xts2<-as.xts(WTIVaralpha2.sstd,date_VaR3)
WTIVaralpha3.sstd.xts2<-as.xts(WTIVaralpha3.sstd,date_VaR3)
plot(logR_WTI_last1002_xts, type = "l", pch = 16, cex = 0.8,  col = "blue",
     ylab = "Returns", main = "Comparison between returns and Var 95% and 99%", xaxt = "n",ylim=c(-12,6))
lines(WTIVaralpha1.sstd.xts2, col = "purple")
lines(WTIVaralpha2.sstd.xts2, col = "red")
lines(WTIVaralpha3.sstd.xts2, col = "red")

BackTestVaR1_WTI_GARCH_sstd<-BacktestVaR(data =logR_WTI_last100_ts,VaR=WTIVaralpha1.sstd,alpha = 0.01,Lags =6)
BackTestVaR1_WTI_GARCH_sstd
BackTestVaR5_WTI_GARCH_sstd<-BacktestVaR(data =logR_WTI_last100_ts,VaR=WTIVaralpha2.sstd,alpha = 0.05,Lags =6)
BackTestVaR5_WTI_GARCH_sstd

plot(logR_WTI_last100_ts, type = "l", pch = 16, cex = 0.8,  col = gray(0.2, 0.5),
     ylab = "Returns", main = "95% VaR Forecasting", xaxt = "n",ylim=c(-15,6))
lines(WTIVaralpha1.sstd, col = "purple")
legend('topright', c("sGARCH(1,1)"),
       lty=1, col=c("purple"), bty='n', cex=.75)



plot(logR_WTI_last100_ts, type = "l", pch = 16, cex = 0.8,  col = gray(0.2, 0.5),
     ylab = "Returns", main = "95% VaR Forecasting", xaxt = "n")
lines(WTIVaralpha2.sstd, col = "purple")
legend('topright', c("sGARCH(1,1)"),
       lty=1, col=c("purple"), bty='n', cex=.75)

preds_WTI.sstd<-as.data.frame(rollingforecast_Gsstd_WTI)
preds_WTI.sstd$Mu
preds_WTI.sstd$Sigma
preds_WTI.sstd$Realized
garchvol_WTI.sstd<-xts(preds_WTI$Sigma, order.by=as.Date(rownames(preds_WTI.sstd)))
plot.xts(garchvol_WTI.sstd)
#prediction error for the mean: 
e_WTI.sstd<-preds_WTI.sstd$Realized-preds_WTI$Mu
mean(e_WTI.sstd^2)
#prediction error for the variance: 
d_WTI.sstd<-(e_WTI.sstd^2-(preds_WTI.sstd$Sigma)^2)
mean(d_WTI.sstd^2)
#Generalized error distribution:
#####
#fit:
Gged_spec_WTI=ugarchspec(variance.model=list(garchOrder=c(1,1)), mean.model=list(armaOrder=c(0,0)),distribution.model = "ged")
Gged_fit_WTI<-ugarchfit(spec=Gged_spec_WTI, data=logR_WTI_xts,out.sample = 100)
Gged_fit_WTI
WTI.coeff.ged = coef(Gged_fit_WTI)
show(WTI.coeff.ged)
WTI.intconf.ged = confint(Gged_fit_WTI)
show(WTI.intconf.ged)
#VOLATILTIY:
plot(Gged_fit_WTI, which = "all")
sigma_WTI_ged_xts<-as.xts(sigma(Gged_fit_WTI))
plot.xts(sigma_WTI_ged_xts,ylab="sigma(t)", col="blue", main="ESTIMATED VOLATILITY TIME SRIES")#volatilit? stimata dei rendimenti con il GARCH(1,1)ged
#residual analysis: 
resid_Gged_WTI<-Gged_fit_WTI@fit$residuals
acf(resid_Gged_WTI, lag.max=220, main="ACF OF RESIDUALS")#ts of the residuals(errors) is a WN
pacf(resid_Gged_WTI, lag.max=220, main="PACF OF std RESIDUALS")
qqnorm(resid_Gged_WTI)
grid()
qqline(resid_Gged_WTI,lwd=2,col="red")
sq_resid_Gged_WTI <-resid_Gged_WTI^2
acf(sq_resid_Gged_WTI, lag.max = 250, main="ACF OF SQUARE RESIDUALS") #ts of the squared residuals is autocorrelated, then the ts has a conditional heterscedasticity
pacf(sq_resid_Gged_WTI, lag.max = 250,main="PACF OF SQUARE RESIDUALS")
#FORECASTING:
forecast_Gged_fit_WTI <-ugarchforecast(Gged_fit_WTI,data=logR_WTI_xts,n.ahead = 100)
forecast_Gged_fit_WTI
plot(forecast_Gged_fit_WTI, which = 1) #grafico dei forecast della nostra serie storica 
plot(forecast_Gged_fit_WTI, which = 3)#grafico dei forecast della volatilità (Sigma)
#Value at risk for GARCH (ged)
#####
rollingforecast_Gged_WTI= ugarchroll(Gged_spec_WTI, data = logR_WTI_xts,n.ahead = 1, forecast.length = 100,refit.every =10, refit.window = c("recursive"), solver = "hybrid", calculate.VaR = TRUE, VaR.alpha = c(0.01, 0.05), cluster = NULL, keep.coef = TRUE)
show(rollingforecast_Gged_WTI)#we predict from 12 jan 2021 to 4 june 2022
roll.WTI.ged=resume(rollingforecast_Gged_WTI)
rollingforecast_Gged_WTI@forecast$VaR#la serie storica stimata (all'interno di forcrolling) dei VaR
plot(rollingforecast_Gged_WTI,which="all")
plot(rollingforecast_Gged_WTI,which=4)
coef(rollingforecast_Gged_WTI)
#Backtesting


report(rollingforecast_Gged_WTI, type = "VaR", VaR.alpha = 0.01, conf.level = 0.99)
report(rollingforecast_Gged_WTI, type = "VaR", VaR.alpha = 0.05, conf.level = 0.95)#both test should have Reject Null: NO, if kubiec no e cristof yes is ok anyways
report(rollingforecast_Gged_WTI, type = "fpm") #con mean squared error, mean absolute error e directional accuracy
plot.ts(rollingforecast_Gged_WTI@forecast[["VaR"]])
analisys_VaR_WTI.ged=basicStats(rollingforecast_Gged_WTI@forecast[["VaR"]])
show(analisys_VaR_WTI.ged)

WTIVaralpha1.ged=rollingforecast_Gged_WTI@forecast$VaR[,1]
WTIVaralpha2.ged=rollingforecast_Gged_WTI@forecast$VaR[,2]
WTIVaralpha3.ged=rollingforecast_Gged_WTI@forecast$VaR[,3]
WTItest1.ged=VaRTest(alpha = 0.01, WTIVaralpha3.ged, WTIVaralpha1.ged, conf.level =  0.99)
WTItest1.ged$expected.exceed
WTItest5.ged=VaRTest(alpha = 0.05, WTIVaralpha3.ged, WTIVaralpha2.ged, conf.level =  0.95)
WTItest5.ged$expected.exceed

BackTestVaR1_WTI_GARCH_ged<-BacktestVaR(data =logR_WTI_last100_ts,VaR=WTIVaralpha1.ged,alpha = 0.01,Lags =6)
BackTestVaR1_WTI_GARCH_ged
BackTestVaR5_WTI_GARCH_ged<-BacktestVaR(data =logR_WTI_last100_ts,VaR=WTIVaralpha2.ged,alpha = 0.05,Lags =6)
BackTestVaR5_WTI_GARCH_ged

logR_WTI_last1002_xts<-as.xts(logR_WTI_xts["2021-01-12/2021-06-04"],date_VaR3)
WTIVaralpha1.ged.xts2<-as.xts(WTIVaralpha1.ged,date_VaR3)
WTIVaralpha2.ged.xts2<-as.xts(WTIVaralpha2.ged,date_VaR3)
plot(logR_WTI_last1002_xts, type = "l", pch = 16, cex = 0.8,  col = "blue",
     ylab = "Returns", main = "Comparison between returns and Var 95% and 99%", xaxt = "n",ylim=c(-12,6))
lines(WTIVaralpha1.ged.xts2, col = "purple")
lines(WTIVaralpha2.ged.xts2, col = "red")

plot(logR_WTI_last100_ts, type = "l", pch = 16, cex = 0.8,  col = gray(0.2, 0.5),
     ylab = "Returns", main = "95% VaR Forecasting", xaxt = "n",ylim=c(-12,6))
lines(WTIVaralpha1.ged, col = "purple")
legend('topright', c("sGARCH(1,1)"),
       lty=1, col=c("purple"), bty='n', cex=.75)

plot(logR_WTI_last100_ts, type = "l", pch = 16, cex = 0.8,  col = gray(0.2, 0.5),
     ylab = "Returns", main = "95% VaR Forecasting", xaxt = "n")
lines(WTIVaralpha2.ged, col = "purple")
legend('topright', c("sGARCH(1,1)"),
       lty=1, col=c("purple"), bty='n', cex=.75)

preds_WTI.ged<-as.data.frame(rollingforecast_Gged_WTI)
preds_WTI.ged$Mu
preds_WTI.ged$Sigma
preds_WTI.ged$Realized
garchvol_WTI.ged<-xts(preds_WTI$Sigma, order.by=as.Date(rownames(preds_WTI.ged)))
plot.xts(garchvol_WTI.ged)
#prediction error for the mean: 
e_WTI.ged<-preds_WTI.ged$Realized-preds_WTI$Mu
mean(e_WTI.ged^2)
#prediction error for the variance: 
d_WTI.ged<-(e_WTI.ged^2-(preds_WTI.ged$Sigma)^2)
mean(d_WTI.ged^2)
#####
#GARCH-MIDAS
#NORM
#####
#K=6:
#fit
covid_WTI_mv6<-mv_into_mat(logR_WTI_xts["2020-04-10/2021-06-04"],log_covid_xts["2020-03-02/2021-06-14"],K=6,"weekly")
GMIDAS_norm_WTI6<-ugmfit(model="GM",skew="NO",distribution = "norm",daily_ret=logR_WTI_xts["2020-04-11/2021-06-14"],mv_m=covid_WTI_mv6,K=6,out_of_sample = 100)
print(GMIDAS_norm_WTI6)
summary.rumidas(GMIDAS_norm_WTI6)#to see the coefficients
#ESTIMATED VOLATILTIY:
GMIDAS_norm_WTI_est_vol_in_s_WTI_xts<-as.xts(GMIDAS_norm_WTI6$est_vol_in_s)
plot.xts(GMIDAS_norm_WTI_est_vol_in_s_WTI_xts,ylab="ESTIMATED VOLATILITY IN SAMPLE", col="blue", main="ESTIMATED VOLATILITY TIME SRIES")#volatilit? stimata dei rendimenti IN SAMPLE con il GARCH-MIDAS norm
GMIDAS_norm_WTI_est_vol_oos_WTI_xts<-as.xts(GMIDAS_norm_WTI6$est_vol_oos)
plot.xts(GMIDAS_norm_WTI_est_vol_oos_WTI_xts,ylab="ESTIMATED VOLATILITY OUT SAMPLE", col="blue", main="ESTIMATED VOLATILITY TIME SRIES")#volatilit? stimata dei rendimenti OUT SAMPLE con il GARCH-MIDAS norm
#residual analysis: 
#FORECASTING GARCH-MIDAS (PREDICTED)
predicted_VaR_GMIDAS_norm_WTI6<-multi_step_ahead_pred(GMIDAS_norm_WTI6,h=100,logR_WTI_xts["2020-04-06/2021-06-14"])
plot.xts(xts(predicted_VaR_GMIDAS_norm_WTI6, date_VaR3))
#REALIZED the last 100 obs of the logR series:
nrow(logR_WTI_xts["2021-01-12/2021-06-14"])
Observations_WTI<-as.numeric(logR_WTI_xts["2021-01-12/2021-06-14"])#we consider the prediction for the last 100 obs (out of sample)
#####
#VaR GARCH-MIDAS NORM
#####
#ROLLING FORECAST: 
logR_WTI_xts1 = logR_WTI_xts["2020-04-11/2021-06-04"]
last.day=215
window.size=100
g=5
VaR5.GMIDAS.WTI<-vector(mode="numeric")
for (q in seq(1,last.day,g)){
  roll_GMIDAS_Norm_WTI7<-ugmfit(model="GM", skew="NO", distribution="norm",
                             daily_ret=logR_WTI_xts1[q:(q+window.size-1)],
                             mv_m=covid_WTI_mv6[,q:(q+window.size-1)],
                             K=6, out_of_sample = 100)
  summary.rumidas(roll_GMIDAS_Norm_WTI7)
  VaR5.GMIDAS.WTI = append(VaR5.GMIDAS.WTI, as.numeric(roll_GMIDAS_Norm_WTI7[q]$est_vol_oos))
}
#5%
Estimated_WTI_VaR5_norm<-as.numeric((GMIDAS_norm_WTI6$est_vol_oos)*qnorm(0.05))#we consider the prediction for the last 100 obs (out of sample)
Estimated_WTI_VaR5_norm_xts<-as.xts(Estimated_WTI_VaR5_norm,date_VaR3)
Estimated_WTI_VaR5_norm_ts<-ts(Estimated_WTI_VaR5_norm)
plot.xts(logR_WTI_xts["2021-01-19/2021-06-14"],col="blue",main="Comparison between returns and VaR 95% and 99%",ylim = c(-15,10))
lines(Estimated_WTI_VaR5_norm_xts,col="red")#we get a nice graph b/c the VaR is always under the returns
analisys_VaR5_GMIDAS_NORM_WTI=basicStats(Estimated_WTI_VaR5_norm_ts)
show(analisys_VaR5_GMIDAS_NORM_WTI)
#1%
Estimated_WTI_VaR1_norm<-as.numeric((GMIDAS_norm_WTI6$est_vol_oos)*qnorm(0.01))#we consider the prediction for the last 100 obs (out of sample)
Estimated_WTI_VaR1_norm_xts<-as.xts(Estimated_WTI_VaR1_norm,date_VaR3)
Estimated_WTI_VaR1_norm_ts<-ts(Estimated_WTI_VaR1_norm)
plot.xts(logR_WTI_xts["2021-01-12/2021-06-14"],col="blue",main="Comparison between returns and VaR 95% and 99% and 99%",ylim = c(-15,10))
lines(Estimated_WTI_VaR1_norm_xts,col="purple")#we get a nice graph b/c the VaR is always under the returns
analisys_VaR1_GMIDAS_NORM_WTI=basicStats(Estimated_WTI_VaR1_norm_ts)
show(analisys_VaR1_GMIDAS_NORM_WTI)
#BACKTESTING:
#5%
BVaR5_WTI_norm<-BacktestVaR(Observations_WTI,Estimated_WTI_VaR5_norm,alpha = 0.05, Lags = 6)
BVaR5_WTI_norm
Estimated_WTI_VaR1_norm.n<-as.numeric(Estimated_WTI_VaR1_norm)
logR_WTI_last100_ts.n<-as.numeric(logR_WTI_last100_ts)
WTIGMtest1.norm=VaRTest(alpha = 0.01,logR_WTI_last100_ts.n, Estimated_WTI_VaR1_norm.n, conf.level =  0.99)
WTIGMtest1.norm

Estimated_WTI_VaR5_snorm.n<-as.numeric(Estimated_WTI_VaR5_snorm)
logR_WTI_last100_ts.n<-as.numeric(logR_WTI_last100_ts)
WTIGMtest5.norm=VaRTest(alpha = 0.01,logR_WTI_last100_ts.n, Estimated_WTI_VaR1_norm.n, conf.level =  0.99)
WTIGMtest5.norm
#1%
BVaR1_WTI_norm<-BacktestVaR(Observations_WTI,Estimated_WTI_VaR1_norm,alpha = 0.01, Lags = 6)

#LOSS FUNCTION
Loss_WTI_VaR5_GMIDAS_norm<-LossVaR(realized=logR_WTI_xts["2021-01-12/2021-06-14"],evaluated=Estimated_WTI_VaR5_norm,which = 'asymmetricLoss', type = 'normal',tau = 0.05)
Loss_WTI_VaR1_GMIDAS_norm<-LossVaR(realized=logR_WTI_xts["2021-01-12/2021-06-14"],evaluated=Estimated_WTI_VaR1_norm,which = 'asymmetricLoss', type = 'normal',tau = 0.01)
#####
#STUD
#####
#K=6:
GMIDAS_std_WTI6<-ugmfit(model="GM",skew="NO",distribution = "std",daily_ret=logR_WTI_xts["2020-04-11/2021-06-14"],mv_m=covid_WTI_mv6,K=6,out_of_sample = 100)
print(GMIDAS_std_WTI6)
summary.rumidas(GMIDAS_std_WTI6)#to see the coefficients
#ESTIMATED VOLATILTIY:
GMIDAS_1_WTI_est_vol_in_s_WTI_xts<-as.xts(GMIDAS_std_WTI6$est_vol_in_s)
plot.xts(GMIDAS_std_WTI_est_vol_in_s_WTI_xts,ylab="ESTIMATED VOLATILITY IN SAMPLE", col="blue", main="ESTIMATED VOLATILITY TIME SRIES")#volatilit? stimata dei rendimenti IN SAMPLE con il GARCH-MIDAS std
GMIDAS_std_WTI_est_vol_oos_WTI_xts<-as.xts(GMIDAS_std_WTI6$est_vol_oos)
plot.xts(GMIDAS_std_WTI_est_vol_oos_WTI_xts,ylab="ESTIMATED VOLATILITY OUT SAMPLE", col="blue", main="ESTIMATED VOLATILITY TIME SRIES")#volatilit? stimata dei rendimenti OUT SAMPLE con il GARCH-MIDAS std
#residual analysis: 
#####
#VaR GARCH-MIDAS STUD
#####
#FORECASTING VaR GARCH-MIDAS (PREDICTED)
predicted_VaR_GMIDAS_std_WTI6<-multi_step_ahead_pred(GMIDAS_std_WTI6,h=100,logR_WTI_xts["2020-04-06/2021-06-14"])
predicted_VaR_GMIDAS_std_WTI6
#5%
Estimated_WTI_VaR5_std<-as.numeric((GMIDAS_std_WTI6$est_vol_oos)*qstd(0.05))#we consider the prediction for the last 100 obs (out of sample)
Estimated_WTI_VaR5_std_xts<-as.xts(Estimated_WTI_VaR5_std,date_VaR3)
Estimated_WTI_VaR5_std_ts<-ts(Estimated_WTI_VaR5_std)
plot.xts(logR_WTI_xts["2021-01-19/2021-06-04"],col="blue",main="Comparison between returns and VaR 95% and 99%",ylim = c(-15,10))
lines(Estimated_WTI_VaR5_std_xts,col="red")#we get a nice graph b/c the VaR is always under the returns
analisys_VaR5_GMIDAS_std_WTI=basicStats(Estimated_WTI_VaR5_std_ts)
show(analisys_VaR5_GMIDAS_std_WTI)
#1%
Estimated_WTI_VaR1_std<-as.numeric((GMIDAS_std_WTI6$est_vol_oos)*qstd(0.01))#we consider the prediction for the last 100 obs (out of sample)
Estimated_WTI_VaR1_std_xts<-as.xts(Estimated_WTI_VaR1_std,date_VaR3)
Estimated_WTI_VaR1_std_ts<-ts(Estimated_WTI_VaR1_std)
plot.xts(logR_WTI_xts["2021-01-19/2021-06-04"],col="blue",main="Comparison between returns and VaR 95% and 99%",ylim = c(-15,10))
lines(Estimated_WTI_VaR1_std_xts,col="purple")#we get a nice graph b/c the VaR is always under the returns
analisys_VaR1_GMIDAS_std_WTI=basicStats(Estimated_WTI_VaR1_std_ts)
show(analisys_VaR1_GMIDAS_std_WTI)
#BACKTESTING:
#5%
BVaR5_WTI_std<-BacktestVaR(Observations_WTI,Estimated_WTI_VaR5_std,alpha = 0.05, Lags = 6)
BVaR5_WTI_std
#1%
BVaR1_WTI_std<-BacktestVaR(Observations_WTI,Estimated_WTI_VaR1_std,alpha = 0.01, Lags = 6)
BVaR1_WTI_std
#LOSS FUNCTION
Loss_WTI_VaR5_GMIDAS_std<-LossVaR(realized=logR_WTI_xts["2021-01-12/2021-06-14"],evaluated=Estimated_WTI_VaR5_std,which = 'asymmetricLoss', type = 'normal',tau = 0.05)
Loss_WTI_VaR1_GMIDAS_std<-LossVaR(realized=logR_WTI_xts["2021-01-12/2021-06-14"],evaluated=Estimated_WTI_VaR1_std,which = 'asymmetricLoss', type = 'normal',tau = 0.01)
#####
#MODEL CONFIDENCE SET PROCEDURE:
#####
#create LOSS function for GARCH:
models <- c("sGARCH")
distributions <- c("norm","snorm" , "std", "sstd", "ged")
spec.comp <- list()
for( m in models ) {
  for( d in distributions ) {
    spec.comp[[paste( m, d, sep = "-" )]] <-
      ugarchspec(mean.model = list(armaOrder = c(0, 0)),
                 variance.model = list(model = m, garchOrder = c(1, 1)),
                 distribution.model=d)
  }
}

specifications <- names( spec.comp )
specifications


roll.comp_WTI <- list()
for( s in specifications ){
  roll.comp_WTI[[s]] <- ugarchroll(spec = spec.comp[[s]], data = logR_WTI,
                                   forecast.length = 100, refit.every = 10,refit.window = c("recursive"))
} 
#1%
VaR1.comp_WTI=list()
for( s in specifications ) {
  VaR1.comp_WTI[[s]] <- as.data.frame(roll.comp_WTI[[s]], which = "VaR")[, 1]
}
Loss1_WTI_GARCH<- do.call(cbind,lapply(specifications,
                                       function(s) LossVaR(tau=0.01, realized=tail(logR_WTI, 100),
                                                           evaluated=VaR1.comp_WTI[[s]])))

colnames(Loss1_WTI_GARCH) <- specifications
print(Loss1_WTI_GARCH)
#5%
VaR5.comp_WTI=list()
for( s in specifications ) {
  VaR5.comp_WTI[[s]] <- as.data.frame(roll.comp_WTI[[s]], which = "VaR")[, 2]
}
Loss5_WTI_GARCH<- do.call(cbind,lapply(specifications,
                                       function(s) LossVaR(tau=0.05, realized=tail(logR_WTI, 100),
                                                           evaluated=VaR5.comp_WTI[[s]])))

colnames(Loss5_WTI_GARCH) <- specifications
print(Loss5_WTI_GARCH)


#create LOSS function for GARCH-MIDAS:
modelss <- c("GARCH-MIDAS")
skew<-c("YES SKEW","NO SKEW")
distributionss <- c("norm","std")

spec.comps<-list()
for (m in modelss) {
  for (s in skew)  {
    for (d in distributionss) {
      spec.comps[[paste(m, s, d, sep="-")]]<-
        (list(modelss=m,skew=s,distributionss = d)) }}}
specificationss <- names(spec.comps)
#1%
LOSS1_WTI_GARCH_MIDAS<-cbind(Loss_WTI_VaR1_GMIDAS_snorm,Loss_WTI_VaR1_GMIDAS_sstd,Loss_WTI_VaR1_GMIDAS_norm,Loss_WTI_VaR1_GMIDAS_std) 
colnames(LOSS1_WTI_GARCH_MIDAS) <- specificationss#ORDER: SNORM,SSTUD,NORM,STUD
print(LOSS1_WTI_GARCH_MIDAS)
#5%
LOSS5_WTI_GARCH_MIDAS<-cbind(Loss_WTI_VaR5_GMIDAS_snorm,Loss_WTI_VaR5_GMIDAS_sstd,Loss_WTI_VaR5_GMIDAS_norm,Loss_WTI_VaR5_GMIDAS_std) 
colnames(LOSS5_WTI_GARCH_MIDAS) <- specificationss#ORDER: SNORM,SSTUD,NORM,STUD
print(LOSS5_WTI_GARCH_MIDAS)
#Create unique loss matrix for GARCH and GARCH-MIDAS togheter:
#1%
LOSS1_WTI_MCS<-cbind(Loss1_WTI_GARCH,LOSS1_WTI_GARCH_MIDAS)
print(LOSS1_WTI_MCS)
#5%
LOSS5_WTI_MCS<-cbind(Loss5_WTI_GARCH,LOSS5_WTI_GARCH_MIDAS)
print(LOSS5_WTI_MCS)
#MCS:
#1%
MCS1_WTI<-MCSprocedure(Loss=LOSS1_WTI_MCS, alpha=0.01,B=5000,cl=NULL,statistic = "Tmax")
#5%
MCS5_WTI<-MCSprocedure(Loss=LOSS5_WTI_MCS, alpha=0.05,B=5000,cl=NULL,statistic = "Tmax")

