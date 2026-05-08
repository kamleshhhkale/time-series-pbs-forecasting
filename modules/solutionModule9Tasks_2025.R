
library(TSA)
library(fUnitRoots)
library(forecast)
library(lmtest)
library(fGarch)
library(rugarch)
library(tseries)

sort.score <- function(x, score = c("bic", "aic")){
  if (score == "aic"){
    x[with(x, order(AIC)),]
  } else if (score == "bic") {
    x[with(x, order(BIC)),]
  } else {
    warning('score = "x" only accepts valid arguments ("aic","bic")')
  }
}

residual.analysis <- function(model, std = TRUE,start = 2, class = c("ARIMA","GARCH","ARMA-GARCH", "garch", "fGARCH")[1]){
  library(TSA)
  library(FitAR)
  if (class == "ARIMA"){
    if (std == TRUE){
      res.model = rstandard(model)
    }else{
      res.model = residuals(model)
    }
  }else if (class == "GARCH"){
    res.model = model$residuals[start:model$n.used]
  }else if (class == "garch"){
    res.model = model$residuals[start:model$n.used]
  }else if (class == "ARMA-GARCH"){
    res.model = model@fit$residuals
  }else if (class == "fGARCH"){
    res.model = model@residuals
  }else {
    stop("The argument 'class' must be either 'ARIMA' or 'GARCH' ")
  }
  par(mfrow=c(3,2))
  plot(res.model,type='o',ylab='Standardised residuals', main="Time series plot of standardised residuals")
  abline(h=0)
  hist(res.model,main="Histogram of standardised residuals")
  qqnorm(res.model,main="QQ plot of standardised residuals")
  qqline(res.model, col = 2)
  acf(res.model,main="ACF of standardised residuals")
  print(shapiro.test(res.model))
  k=0
  LBQPlot(res.model, lag.max = 30, StartLag = k + 1, k = 0, SquaredQ = FALSE)
  par(mfrow=c(1,1))
}



# --- TASK 1 ---.
#  Google stock data from August 1, 2014 to September 13, 2016 from the google dataset stored on the TSA package
data("google") # This is already returns series
head(google)
par(mfrow=c(1,1))
plot(google,type='o',main="Time series plot of daily returns of Google stock")
# There is sign of neither a trend nor seasonality. Observations are bouncing around the mean level. But changing variance is obvious.
acf(google)
summary(google) # Mean is very close to zero

adf.test(google)
pp.test(google)

McLeod.Li.test(y=google,main="McLeod-Li Test Statistics for Daily Google Returns")
# McLeod-Li test is significant at 5% level of significance for all lags. This gives a strong idea about existence of volatiliy clustering.

par(mfrow=c(1,2))
acf(google, main="The ACF plot for return series")
pacf(google, main="The PACF plot for return series")
# ACF, PACF show pattern of white noise for the correlation structure. However, there is an ARCH effect present in the series.

#So we'll use absolute value and square transformations to figure out this ARCH effect.
abs.google = abs(google)
sq.google = google^2

par(mfrow=c(1,2))
acf(abs.google, ci.type="ma",main="The ACF plot for absolute return series")
pacf(abs.google, main="The PACF plot for absolute return series")
# max(p, q) = 2, q = 3 ==> max(p, q = 3) = 2 ==> max(p, 3) = 2 <> No models

eacf(abs.google)
# From the EACF, the top left 'o', not distracted by an 'x' in the same row is at (1,1) for absolute 
# value series. 
# These models correspond to parameter settings:
# max(p,q) = 1 and q = 1 ==> max(p,1) = 1 and q = 1 ==> p can be either 0 or 1.
# max(p,q) = 1 and q = 2 ==> max(p,2) = 1 and q = 2 <> There is no model.
# max(p,q) = 2 and q = 2 ==> max(p,2) = 2 and q = 2 ==> p can be 0, 1 or 2.
# So the corresponding tentative GARCH models are 
# {GARCH(0,1), GARCH(1,1), GARCH(0,2), GARCH(1,2), and GARCH(2,2)}.

par(mfrow=c(1,2))
acf(sq.google, ci.type="ma",main="The ACF plot for squared return series")
pacf(sq.google, main="The PACF plot for squared return series")
# max(p, q) = 1, q = 2 ==> max(p, q = 2) = 1 ==> max(p, 2) = 1 <> No models

eacf(sq.google)
# From the EACF, the top left 'o', not distracted by an 'x' in the same row is at (1,1) for squared series. 
# These models correspond to parameter settings:
# max(p,q) = 1 and q = 1 ==> max(p,1) = 1 and q = 1 ==> p can be either 0 or 1.
# max(p,q) = 1 and q = 2 ==> max(p,2) = 1 and q = 2 <> There is no model.
# max(p,q) = 2 and q = 1 ==> max(p,1) = 2 and q = 1 ==> p can only be 2.
# max(p,q) = 2 and q = 2 ==> max(p,2) = 2 and q = 2 ==> p can be 0, 1 and 2.
# The corresponding tentative GARCH models are 
# {GARCH(0,1), GARCH(1,1), GARCH(2,1), GARCH(0,2), GARCH(1,2), and GARCH(2,2)}.

# Overall, we have {GARCH(0,1), GARCH(1,1), GARCH(2,1), GARCH(0,2), GARCH(1,2), GARCH(2,2)}.

m.01 = garch(google,order=c(0,1),trace = FALSE)
summary(m.01) 
residual.analysis(model = m.01, class= "garch") # use class = "garch" for garch() function from tseries package

m.11 = garch(google,order=c(1,1),trace = FALSE)
summary(m.11) 
residual.analysis(model = m.11, class= "garch")

m.21 = garch(google,order=c(2,1),trace = FALSE)
summary(m.21)
residual.analysis(model = m.21, class = "garch") # Does not work with the default "start = 2". So we increase "start".
residual.analysis(model = m.21, start = 3, class= "garch")

m.02 = garch(google,order=c(0,2),trace = FALSE)
summary(m.02)
residual.analysis(model = m.02, start = 3, class = "garch")

m.12 = garch(google,order=c(1,2),trace = FALSE)
summary(m.12)
residual.analysis(model = m.12, start = 3, class = "garch")

m.22 = garch(google,order=c(2,2),trace = FALSE)
summary(m.22)
residual.analysis(model = m.22, start = 3, class = "garch")

sc.AIC=AIC(m.01, m.11, m.21, m.02, m.12, m.22)
sc.BIC=AIC(m.01, m.11, m.21, m.02, m.12, m.22, k = log(length(google)))

sort.score(sc.AIC, score = "aic")
sort.score(sc.BIC, score = "aic")

# GARCH(1,1) model is the best one in terms of residual analysis and both AIC and BIC.

par(mfrow=c(1,1))
plot((fitted(m.11)[,1])^2,type='l',ylab='Conditional Variance',xlab='t',main="Estimated Conditional Variances of the Daily Returns")
# Changes in conditional variance at the beginning of the series and between observations 300 and 400, then the conditional variance settles down. 

# Forecasts for the confidence limits are based on the forecasts of conditional variance.
library(rugarch)
model<-ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)), 
                  mean.model = list(armaOrder = c(0, 0), include.mean = FALSE), 
                  distribution.model = "norm")
m.11R<-ugarchfit(spec=model,data=google)
par(mfrow=c(1,1))
plot(m.11R)
plot(m.11R, which = 2)

par(mfrow=c(1,1))
forc = ugarchforecast(m.11R, data = google, n.ahead = 10)
plot(forc, which = 1)
plot(forc, which = 3)


# --- TASK 2 ---
# Lg wave seismic data from an earthquake known as the Massachusetts Mountain Earthquake (August 5, 1971)
library(tswge)
data("mm.eq")
wave = mm.eq
wave = ts(wave)
head(wave)
par(mfrow=c(1,1))
plot(wave,main="Time series plot of seismic Lg wave series")
# There is sign of neither a trend nor seasonality. Observations are bouncing around the mean level. But changing variance is obvious especially 
# at the beginning of the observation period.
summary(wave)

qqnorm(wave, ylab="Seismic Lg wave ", xlab="Normal Scores")
qqline(wave)
shapiro.test(wave) 

par(mfrow=c(1,2))
acf(wave, main="The sample ACF plot for seismic 
    Lg wave series")
pacf(wave, main="The sample PACF plot for seismic 
     Lg wave series")
# From the ACF and PACF plots we observe significant correlations and there is no sign of 
# a white noise process. However, volatility clustering is obvious in the time series plot. So, 
# we will consider fitting an ARMA+GARCH model.

wave.positive = wave + min(abs(wave))+0.1
waveReturn = diff(log(wave.positive))

par(mfrow=c(1,1))
plot(waveReturn, main="Time series plot of returns series for seismic Lg wave series")

adf.test(waveReturn)

par(mfrow=c(1,1))
McLeod.Li.test(y=waveReturn, main="McLeod-Li Test Statistics for seismic 
               Lg wave series")
# McLeod-Li test is significant at 5% level of significance for all lags. This gives a strong idea about existence of volatiliy clustering.

# First ARMA part. 

par(mfrow=c(1,2))
acf(waveReturn, main="The ACF plot of returns series for seismic Lg wave series")
pacf(waveReturn, main="The PACF plot of returns series for seismic Lg wave series")
par(mfrow=c(1,1))
# Due to the changing variance, we do not have a clear picture in ACF and PACF.
# Referring very high autocorrelations {ARMA(0,3), ARMA(0,4), ARMA(0,5), ARMA(0,6)} can be identified

eacf(waveReturn)
#{ARMA(5,4), ARMA(5,5), ARMA(6,5)}

res = armasubsets(y=waveReturn,nar=14,nma=14,y.name='p',ar.method='ols')
plot(res)
# {ARMA(1,1), ARMA(2,1), ARMA(1,3), ARMA(2,3), ARMA(1,4), ARMA(2,4)}

#Overall, 
#{ARMA(0,3), ARMA(0,4), ARMA(0,5), ARMA(0,6), ARMA(5,4), ARMA(5,5), 
# ARMA(6,5), ARMA(1,1), ARMA(2,1), ARMA(1,3), ARMA(2,3), ARMA(1,4), 
# ARMA(2,4)}

# ARMA(0,3)
model_03_css = Arima(waveReturn,order=c(0,0,3),method='CSS')
coeftest(model_03_css)
residual.analysis(model = model_03_css)


model_03_ml = Arima(waveReturn,order=c(0,0,3),method='ML')
coeftest(model_03_ml)
residual.analysis(model = model_03_ml)

# ARMA(0,4)
model_04_css = Arima(waveReturn,order=c(0,0,4),method='CSS')
coeftest(model_04_css)
residual.analysis(model = model_04_css)

model_04_ml = Arima(waveReturn,order=c(0,0,4),method='ML')
coeftest(model_04_ml)
residual.analysis(model = model_04_ml)

# ARMA(0,5)
model_05_css = Arima(waveReturn,order=c(0,0,5),method='CSS')
coeftest(model_05_css)
residual.analysis(model = model_05_css)

model_05_ml = Arima(waveReturn,order=c(0,0,5),method='ML')
coeftest(model_05_ml)
residual.analysis(model = model_05_ml)

# ARMA(0,6)
model_06_css = Arima(waveReturn,order=c(0,0,6),method='CSS')
coeftest(model_06_css)
residual.analysis(model = model_06_css)

model_06_ml = Arima(waveReturn,order=c(0,0,6),method='ML')
coeftest(model_06_ml)
residual.analysis(model = model_06_ml)
# plot(forecast(model_06_ml,15))

# ARMA(5,4)
model_54_css = Arima(waveReturn,order=c(5,0,4),method='CSS')
coeftest(model_54_css)
residual.analysis(model = model_54_css)

model_54_ml = Arima(waveReturn,order=c(5,0,4),method='ML')
coeftest(model_54_ml)
residual.analysis(model = model_54_ml)

# ARMA(5,5)
model_55_css = Arima(waveReturn,order=c(5,0,5),method='CSS')
coeftest(model_55_css)
residual.analysis(model = model_55_css)

model_55_ml = Arima(waveReturn,order=c(5,0,5),method='ML')
coeftest(model_55_ml)
residual.analysis(model = model_55_ml)

# ARMA(6,5)
model_65_css = Arima(waveReturn,order=c(6,0,5),method='CSS')
coeftest(model_65_css)
residual.analysis(model = model_65_css)

model_65_ml = Arima(waveReturn,order=c(6,0,5),method='ML')
coeftest(model_65_ml)
residual.analysis(model = model_65_ml)

# ARMA(1,1)
model_11_css = Arima(waveReturn,order=c(1,0,1),method='CSS')
coeftest(model_11_css)
residual.analysis(model = model_11_css)

model_11_ml = Arima(waveReturn,order=c(1,0,1),method='ML')
coeftest(model_11_ml)
residual.analysis(model = model_11_ml)

# ARMA(2,1)
model_21_css = Arima(waveReturn,order=c(2,0,1),method='CSS')
coeftest(model_21_css)
residual.analysis(model = model_21_css)

model_21_ml = Arima(waveReturn,order=c(2,0,1),method='ML')
coeftest(model_21_ml)
residual.analysis(model = model_21_ml)

# ARMA(1,3)
model_13_css = Arima(waveReturn,order=c(1,0,3),method='CSS')
coeftest(model_13_css)
residual.analysis(model = model_13_css)

model_13_ml = Arima(waveReturn,order=c(1,0,3),method='ML')
coeftest(model_13_ml)
residual.analysis(model = model_13_ml)

# ARMA(2,3)
model_23_css = Arima(waveReturn,order=c(2,0,3),method='CSS')
coeftest(model_23_css)
residual.analysis(model = model_23_css)

model_23_ml = Arima(waveReturn,order=c(2,0,3),method='ML')
coeftest(model_23_ml)
residual.analysis(model = model_23_ml)

# ARMA(1,4)
model_14_css = Arima(waveReturn,order=c(1,0,4),method='CSS')
coeftest(model_14_css)
residual.analysis(model = model_14_css)

model_14_ml = Arima(waveReturn,order=c(1,0,4),method='ML')
coeftest(model_14_ml)
residual.analysis(model = model_14_ml)

# ARMA(2,4)
model_24_css = Arima(waveReturn,order=c(2,0,4),method='CSS')
coeftest(model_24_css)
residual.analysis(model = model_24_css)

model_24_ml = Arima(waveReturn,order=c(2,0,4),method='ML')
coeftest(model_24_ml)
residual.analysis(model = model_24_ml)


sc.AIC=AIC(model_03_ml, model_04_ml, model_05_ml, model_06_ml, model_54_ml, model_55_ml, 
           model_65_ml, model_11_ml, model_21_ml, model_13_ml, model_23_ml, model_14_ml, 
           model_24_ml)
sc.BIC=AIC(model_03_ml, model_04_ml, model_05_ml, model_06_ml, model_54_ml, model_55_ml, 
           model_65_ml, model_11_ml, model_21_ml, model_13_ml, model_23_ml, model_14_ml, 
           model_24_ml, k = log(length(waveReturn)))

sort.score(sc.AIC, score = "aic")
sort.score(sc.BIC, score = "aic")

sc.AIC=AIC(model_03_ml, model_04_ml, model_05_ml, model_06_ml, model_54_ml, model_55_ml, 
           model_65_ml, model_11_ml, model_21_ml, model_13_ml, model_23_ml, model_14_ml, 
           model_24_ml, model_16_ml)
sc.BIC=AIC(model_03_ml, model_04_ml, model_05_ml, model_06_ml, model_54_ml, model_55_ml, 
           model_65_ml, model_11_ml, model_21_ml, model_13_ml, model_23_ml, model_14_ml, 
           model_24_ml, model_16_ml, k = log(length(waveReturn)))

sort.score(sc.AIC, score = "aic")
sort.score(sc.BIC, score = "aic")

Smodel_03_ml <- accuracy(model_03_ml)[1:7]
Smodel_04_ml <- accuracy(model_04_ml)[1:7]
Smodel_05_ml <- accuracy(model_05_ml)[1:7]
Smodel_06_ml <- accuracy(model_06_ml)[1:7]
Smodel_54_ml <- accuracy(model_54_ml)[1:7]
Smodel_55_ml <- accuracy(model_55_ml)[1:7]
Smodel_65_ml <- accuracy(model_65_ml)[1:7]
Smodel_11_ml <- accuracy(model_11_ml)[1:7]
Smodel_21_ml <- accuracy(model_21_ml)[1:7]
Smodel_13_ml <- accuracy(model_13_ml)[1:7]
Smodel_23_ml <- accuracy(model_23_ml)[1:7]
Smodel_14_ml <- accuracy(model_14_ml)[1:7]
Smodel_24_ml <- accuracy(model_24_ml)[1:7]

df.Smodels <- data.frame(
  rbind(Smodel_03_ml, Smodel_04_ml, Smodel_05_ml, Smodel_06_ml, Smodel_54_ml, Smodel_55_ml, 
        Smodel_65_ml, Smodel_11_ml, Smodel_21_ml, Smodel_13_ml, Smodel_23_ml, Smodel_14_ml, 
        Smodel_24_ml)
)
colnames(df.Smodels) <- c("ME", "RMSE", "MAE", "MPE", "MAPE", 
                          "MASE", "ACF1")
rownames(df.Smodels) <- c("ARIMA(0,0,3)", "ARIMA(0,0,4)", "ARIMA(0,0,5)", "ARIMA(0,0,6)", "ARIMA(5,0,4)", "ARIMA(5,0,5)",
                          "ARIMA(6,0,5)", "ARIMA(1,0,1)","ARIMA(2,0,1)", "ARIMA(1,0,3)", "ARIMA(2,0,3)", "ARIMA(1,0,4)",
                          "ARIMA(2,0,4)")
round(df.Smodels,  digits = 3)


# ARMA(0,6) model is the best one in this set. Overfitting models are ARMA(1,6) and ARMA(0,7)

model_16_css = Arima(waveReturn,order=c(1,0,6),method='CSS')
coeftest(model_16_css)
residual.analysis(model = model_16_css)

model_16_ml = Arima(waveReturn,order=c(1,0,6),method='ML')
coeftest(model_16_ml)
residual.analysis(model = model_16_ml)


model_07_css = Arima(waveReturn,order=c(0,0,7),method='CSS')
coeftest(model_07_css)
residual.analysis(model = model_07_css)

model_07_ml = Arima(waveReturn,order=c(0,0,7),method='ML')
coeftest(model_07_ml)
residual.analysis(model = model_07_ml)


sc.AIC=AIC(model_03_ml, model_04_ml, model_05_ml, model_06_ml, model_54_ml, model_55_ml, 
           model_65_ml, model_11_ml, model_21_ml, model_13_ml, model_23_ml, model_14_ml, 
           model_24_ml, model_16_ml)
sc.BIC=AIC(model_03_ml, model_04_ml, model_05_ml, model_06_ml, model_54_ml, model_55_ml, 
           model_65_ml, model_11_ml, model_21_ml, model_13_ml, model_23_ml, model_14_ml, 
           model_24_ml, model_16_ml, k = log(length(waveReturn)))

sort.score(sc.AIC, score = "aic")
sort.score(sc.BIC, score = "aic")

Smodel_03_ml <- accuracy(model_03_ml)[1:7]
Smodel_04_ml <- accuracy(model_04_ml)[1:7]
Smodel_05_ml <- accuracy(model_05_ml)[1:7]
Smodel_06_ml <- accuracy(model_06_ml)[1:7]
Smodel_54_ml <- accuracy(model_54_ml)[1:7]
Smodel_55_ml <- accuracy(model_55_ml)[1:7]
Smodel_65_ml <- accuracy(model_65_ml)[1:7]
Smodel_11_ml <- accuracy(model_11_ml)[1:7]
Smodel_21_ml <- accuracy(model_21_ml)[1:7]
Smodel_13_ml <- accuracy(model_13_ml)[1:7]
Smodel_23_ml <- accuracy(model_23_ml)[1:7]
Smodel_14_ml <- accuracy(model_14_ml)[1:7]
Smodel_24_ml <- accuracy(model_24_ml)[1:7]
Smodel_16_ml <- accuracy(model_16_ml)[1:7]

df.Smodels <- data.frame(
  rbind(Smodel_03_ml, Smodel_04_ml, Smodel_05_ml, Smodel_06_ml, Smodel_54_ml, Smodel_55_ml, 
        Smodel_65_ml, Smodel_11_ml, Smodel_21_ml, Smodel_13_ml, Smodel_23_ml, Smodel_14_ml, 
        Smodel_24_ml, Smodel_16_ml)
)
colnames(df.Smodels) <- c("ME", "RMSE", "MAE", "MPE", "MAPE", 
                          "MASE", "ACF1")
rownames(df.Smodels) <- c("ARIMA(0,0,3)", "ARIMA(0,0,4)", "ARIMA(0,0,5)", "ARIMA(0,0,6)", "ARIMA(5,0,4)", "ARIMA(5,0,5)",
                          "ARIMA(6,0,5)", "ARIMA(1,0,1)","ARIMA(2,0,1)", "ARIMA(1,0,3)", "ARIMA(2,0,3)", "ARIMA(1,0,4)",
                          "ARIMA(2,0,4)","ARIMA(1,0,6)")
round(df.Smodels,  digits = 3)





# We will use ARMA(0,6) for the ARMA part of the model.

m06Residuals = model_06_ml$residuals
abs.res = abs(m06Residuals)
sq.res = m06Residuals^2

par(mfrow=c(1,2))
acf(abs.res, main="The ACF plot for absolute residual series")
pacf(abs.res, main="The PACF plot for absolute residual series")
par(mfrow=c(1,1))
# max(p,q) = 2, q = 0 ==>  max(p, q = 0) = 2 ==>  max(p, 0) = 2 ==> p can only be 2
# So we get GARCH(2,0)

eacf(abs.res)
# Top-left-o is at (3,1)
# max(p,q) = 3 and q = 1 ==> max(p,q = 1) = 3 ==> max(p,1) = 3 ==> p can only be 3.
# max(p,q) = 3 and q = 2 ==> max(p,q = 2) = 3 ==> max(p,2) = 3 ==> p can only be 3.

# The corresponding tentative GARCH models are {GARCH(2,0), GARCH(3,1), GARCH(3,2)}.


par(mfrow=c(1,2))
acf(sq.res, main="The ACF plot for square residual series")
pacf(sq.res, main="The PACF plot for square residual series")
par(mfrow=c(1,1))
# max(p,q) = 3, q = 2 ==>  max(p, q = 2) = 3 ==>  max(p, 2) = 3 ==> p can only be 3
# So we get GARCH(3,2)

eacf(sq.res)
# Top-left-o is at (0,2)
# max(p,q) = 0 and q = 2 ==> max(p,q = 2) = 0 ==> max(p,2) = 0 <> No models.
# max(p,q) = 0 and q = 3 ==> max(p,q = 3) = 0 ==> max(p,3) = 0 <> No models.
# max(p,q) = 1 and q = 3 ==> max(p,q = 3) = 1 ==> max(p,3) = 1 <> No models..

# Overall set of possible GARCH models is {GARCH(2,0), GARCH(3,1), GARCH(3,2)}.

# When combined with the best model for the ARMA part, we have
# {ARMA(0,6)+GARCH(2,0), ARMA(0,6)+GARCH(3,1), ARMA(0,6)+GARCH(3,2)}
#                                                                         q , p!
model1<-ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(0, 2)), 
                   mean.model = list(armaOrder = c(0, 6), include.mean = FALSE), 
                   distribution.model = "snorm") # The conditional density to use for the innovations. Try with "norm"
m.06_20<-ugarchfit(spec = model1, data = waveReturn, out.sample = 100)
m.06_20
residual.analysis(model=m.06_20, class = "ARMA-GARCH")

                                                                      #  q , p!
model2<-ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 3)), 
                   mean.model = list(armaOrder = c(0, 6), include.mean = FALSE), 
                   distribution.model = "norm")
m.06_31<-ugarchfit(spec = model2, data = waveReturn, out.sample = 100)
m.06_31
residual.analysis(model=m.06_31, class = "ARMA-GARCH")


model3<-ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(2, 3)), 
                   mean.model = list(armaOrder = c(0, 6), include.mean = FALSE), 
                   distribution.model = "norm")
m.06_32<-ugarchfit(spec = model3, data = waveReturn, out.sample = 100)
m.06_32
residual.analysis(model=m.06_32, class = "ARMA-GARCH")


# Final model that gives best diagnostics is ARMA(0,6) + GARCH(2,0)

AICs = c(infocriteria(m.06_20)[1], infocriteria(m.06_31)[1], infocriteria(m.06_32)[1])
BICs = c(infocriteria(m.06_20)[2], infocriteria(m.06_31)[2], infocriteria(m.06_32)[2])
df = data.frame(AIC = AICs, BIC = BICs)   
rownames(df) = c("m.06_20", "m.06_31", "m.06_32")

df[order(df$AIC),] # See sorted by AIC
df[order(df$BIC),] # See sorted by BIC

par(mfrow=c(1,1))
# plot(m.06_20)
plot(m.06_20, which = 1)

par(mfrow=c(1,1))
forc = ugarchforecast(m.06_20, data = waveReturn, n.ahead = 10)
plot(forc, which = 1)
plot(forc, which = 3)
forc

# All the transformations:
# wave.positive = wave + min(abs(wave))+0.1
# waveReturn = diff(log(wave.positive))
# Take them back
# for raw series
firstObs <- matrix(c(log(wave.positive)[1]),1)
log.wave.diff1.back = diffinv(waveReturn, xi = firstObs)
log.wave.diff1.back = exp(log.wave.diff1.back)
log.wave.diff1.back.original = log.wave.diff1.back - (min(abs(wave))+0.1)
log.wave.diff1.back.original - wave # Make sure you are doing it correctly!

# Then apply it to the forecasts:
frc <- forc@forecast$seriesFor
lastObs <- matrix(c(log(wave.positive)[454]),1)
log.wave.diff1.back = diffinv(frc, xi = lastObs)
log.wave.diff1.back = exp(log.wave.diff1.back)
log.wave.diff1.back.frc = log.wave.diff1.back - (min(abs(wave))+0.1)

plot(wave, xlim= c(1, 454+11), ylim = c(min(wave,log.wave.diff1.back.frc), 
                                        max(wave,log.wave.diff1.back.frc)), 
     ylab = "Seismic Lg wave series",
     main = "Forecasts from ARMA+GARCH model.")
lines(ts(as.vector(log.wave.diff1.back.frc), start = c(455)), col="blue", type="l")
legend("topleft", lty=1, pch=1, col=c("black","blue"), text.width = 18,
       c("Data", "Forecasts"))


