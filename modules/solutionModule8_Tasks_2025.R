
library(TSA)
library(fUnitRoots)
library(forecast)
# library(CombMSC)
library(lmtest)
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
residual.analysis <- function(model, std = TRUE,start = 2, class = c("ARIMA","GARCH","ARMA-GARCH", "fGARCH")[1]){
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
  seasonal_acf(res.model,main="ACF of standardised residuals")
  print(shapiro.test(res.model))
  k=0
  LBQPlot(res.model, lag.max = 30, StartLag = k + 1, k = 0, SquaredQ = FALSE)
  par(mfrow=c(1,1))
}


##################################################
# Following functions are developed by               #  
# MATH1318 students                              #
# Le Van Tra Tran and Tin Trung Pham             #
# in 2024. WE thank them for their contribution! #
##################################################


helper <- function(class = c("acf", "pacf"), ...) {
  
  # Capture additional arguments
  params <- match.call(expand.dots = TRUE)
  params <- as.list(params)[-1]
  
  # Calculate ACF/PACF values
  if (class == "acf") {
    acf_values <- do.call(acf, c(params, list(plot = FALSE)))
  } else if (class == "pacf") {
    acf_values <- do.call(pacf, c(params, list(plot = FALSE)))
  }
  
  # Extract values and lags
  acf_data <- data.frame(
    Lag = as.numeric(acf_values$lag),  
    ACF = as.numeric(acf_values$acf)   
  )
  
  # Identify seasonal lags to be highlighted
  seasonal_lags <- acf_data$Lag %% 1 == 0
  
  # Plot ACF/PACF values
  if (class == "acf") {
    do.call(acf, c(params, list(plot = TRUE)))
  } else if (class == "pacf") {
    do.call(pacf, c(params, list(plot = TRUE)))
  }
  
  # Add colored segments for seasonal lags
  for (i in which(seasonal_lags)) {
    segments(x0 = acf_data$Lag[i], y0 = 0, x1 = acf_data$Lag[i], y1 = acf_data$ACF[i], col = "red")
  }
}


# seasonal_acf ------------------------------------------------------------

seasonal_acf <- function(...) {
  helper(class = "acf", ...)
}


# seasonal_pacf -----------------------------------------------------------

seasonal_pacf <- function(...) {
  helper(class = "pacf", ...)
}
#########################################

#--- Task 1 ---

NMFS_Landings <- read.csv('/Users/kamleshhhkale/Documents/RMIT/SEM 3/TIME SERIES ANALYSIS/NMFS_Landings.csv')
class(NMFS_Landings)
head(NMFS_Landings)
# Convert data into a time series object
NMFS_Landings.ts = matrix(NMFS_Landings$Metric_Tons, nrow = 25, ncol = 12)
NMFS_Landings.ts = as.vector(t(NMFS_Landings.ts))
NMFS_Landings.ts = ts(NMFS_Landings.ts,start=c(1991,1), end=c(2015,12), frequency=12)
class(NMFS_Landings.ts)

par(mfrow=c(1,1))
plot(NMFS_Landings.ts,ylab='Landings in metric tons',xlab='Year',type='o', main = "Time series plot of monthly landings in metric tons.")
# There are two different periods of time that the series are bouncing around two different mean levels. Although there is not a trend for each 
# time period, as a whole the mean level is changing through the time. Seasonality and changing variance are obvious from the series. 
# To see the autocorrelation structure clearly, we need to filtrate out the effect of seasonality.

par(mfrow=c(1,2))
seasonal_acf(NMFS_Landings.ts,  lag.max = 64,main="The sample ACF of landings series.")
seasonal_pacf(NMFS_Landings.ts,  lag.max = 64,main="The sample PACF of landings series.")
par(mfrow=c(1,1))

# Seasonality and existence of trend are obvious from the ACF and PACF plots

adf.test(NMFS_Landings.ts)
pp.test(NMFS_Landings.ts)

BC <- BoxCox.ar(NMFS_Landings.ts) #,lambda = seq(-1, 0.5, 0.01) If you get an error.
BC$ci
lambda <- BC$lambda[which(max(BC$loglike) == BC$loglike)]
lambda # lambda = 0 corresponds to the log transformation!
# BC.value = ((value^lambda)-1)/lambda

# We will first work with the raw series and then work on the log-transformed series.

# First fit a plain model with only the first seasonal difference with order D = 1 
# and see if we can get rid of the seasonal trend effect
# by inspecting the autocorrelation structure of the residuals.
#                                           p,d,q                        P,D,Q
m1.landing = Arima(NMFS_Landings.ts,order=c(0,0,0),seasonal=list(order=c(0,1,0), period=12))
res.m1 = residuals(m1.landing);  
par(mfrow=c(1,1))
plot(res.m1,xlab='Time',ylab='Residuals',main="Time series plot of the residuals")
par(mfrow=c(1,2))
seasonal_acf(res.m1, lag.max = 48, main = "The sample ACF of the residuals")
seasonal_pacf(res.m1, lag.max = 48, main = "The sample PACF of the residuals")
# From the time series plot, we can conclude that we got rid of the trend. Seasonal autocorrelations are seen clearly in 
# ACF and PACF now at the lags corresponding to the periods. 
# We have one significant correlation at the first seasonal lag in both ACF and PACF. 

#So, we will add the SARMA(1,1) component and see if we get rid of seasonal component.
m2.landing = Arima(NMFS_Landings.ts,order=c(0,0,0),seasonal=list(order=c(1,1,1), period=12))
res.m2 = residuals(m2.landing)  
par(mfrow=c(1,1))
plot(res.m2,xlab='Time',ylab='Residuals',main="Time series plot of the residuals")
par(mfrow=c(1,2))
seasonal_acf(res.m2, lag.max = 36, main = "The sample ACF of the residuals")
seasonal_pacf(res.m2, lag.max = 36, main = "The sample PACF of the residuals")
# Although we have one significant autocorrelation at the first seasonal lag in PACF, we can conclude that the seasonality is filtrated out for now.
# The significant correlation in PACF would be due to the change point in the series.
# Now, we will specify the orders of ARIMA component. The first ordinary correlation in ACF is highly significant, and after a gap, 
# there are four significant autocorrelations before the first seasonal lag. Also, there is a jagged pattern in PACF. This is an indication of
# trend but it's not apparent. This would be due to the high changing variance. 
# So, it would help to apply a transformation to the raw series at this stage. 

# So, we will apply the log transformation and see if we can see the trend more clearly.
m3.landing = Arima(log(NMFS_Landings.ts),order=c(0,0,0),seasonal=list(order=c(1,1,1), period=12))
res.m3 = residuals(m3.landing);  
par(mfrow=c(1,1))
plot(res.m3,xlab='Time',ylab='Residuals',main="Time series plot of the residuals")
par(mfrow=c(1,2))
seasonal_acf(res.m3, lag.max = 36, main = "The sample ACF of the residuals")
seasonal_pacf(res.m3, lag.max = 36, main = "The sample PACF of the residuals")

# The variation in the residuals decreased after the transformation. We have a very high autocorrelation at the first lag of PACF and nearly all the 
# autocorrelations before the first seasonal lag are significant in ACF. So, we need to take the first ordinary difference to get rid of this trend 
# effect before going on with the specification of ARMA orders.

m4.landing = Arima(log(NMFS_Landings.ts),order=c(0,1,0),seasonal=list(order=c(1,1,1), period=12))
res.m4 = residuals(m4.landing);  
par(mfrow=c(1,1))
plot(res.m4,xlab='Time',ylab='Residuals',main="Time series plot of the residuals")
par(mfrow=c(1,2))
seasonal_acf(res.m4, lag.max = 36, main = "The sample ACF of the residuals")
seasonal_pacf(res.m4, lag.max = 36, main = "The sample PACF of the residuals")
# In ACF, there are 6 significant lags and the significant partial autocorrelations in PACF could be seen as a decreasing pattern 
# or we can set p = 6. So, we can add MA component up to order 5. We can use EACF over the residuals.
# Note that the significant partial autocorrelation at lag 1 of PACF has disappeared!

m5.landing = Arima(log(NMFS_Landings.ts),order=c(6,1,4),seasonal=list(order=c(1,1,1), period=12))
res.m5 = residuals(m5.landing);  
par(mfrow=c(1,1))
plot(res.m5,xlab='Time',ylab='Residuals',main="Time series plot of the residuals")
par(mfrow=c(1,2))
seasonal_acf(res.m5, lag.max = 36, main = "The sample ACF of the residuals")
seasonal_pacf(res.m5, lag.max = 36, main = "The sample PACF of the residuals")

eacf(res.m4) # Be careful - residuals of m4 is being used for EACF since there is leftover signal in them.

# The tentative models are specified as 
# SARIMA(0,1,4)x(1,1,1)_12
# SARIMA(0,1,5)x(1,1,1)_12
# SARIMA(1,1,4)x(1,1,1)_12  
# SARIMA(1,1,5)x(1,1,1)_12
# SARIMA(6,1,4)x(1,1,1)_12
par(mfrow=c(1,1))
bic_table = armasubsets(y=res.m4,nar=10,nma=10,y.name='p',ar.method='ols')
plot(bic_table)
# SARIMA(10,1,2)x(1,1,1)_12
# SARIMA(4,1,2)x(1,1,1)_12
# SARIMA(1,1,2)x(1,1,1)_12

# SARIMA(0,1,4)x(1,1,1)_12
m5_014.landing = Arima(log(NMFS_Landings.ts),order=c(0,1,4),seasonal=list(order=c(1,1,1), period=12),method = "ML")
coeftest(m5_014.landing)
residual.analysis(model = m5_014.landing)

m5_014.landingCSS = Arima(log(NMFS_Landings.ts),order=c(0,1,4),seasonal=list(order=c(1,1,1), period=12),method = "CSS")
coeftest(m5_014.landingCSS)
residual.analysis(model = m5_014.landingCSS)

# SARIMA(0,1,5)x(1,1,1)_12
m5_015.landing = Arima(log(NMFS_Landings.ts),order=c(0,1,5),seasonal=list(order=c(1,1,1), period=12),method = "ML")
coeftest(m5_015.landing)
residual.analysis(model = m5_015.landing)

m5_015.landingCSS = Arima(log(NMFS_Landings.ts),order=c(0,1,5),seasonal=list(order=c(1,1,1), period=12),method = "CSS")
coeftest(m5_015.landingCSS)
residual.analysis(model = m5_015.landingCSS)

# SARIMA(1,1,4)x(1,1,1)_12  
m5_114.landing = Arima(log(NMFS_Landings.ts),order=c(1,1,4),seasonal=list(order=c(1,1,1), period=12),method = "ML")
coeftest(m5_114.landing)
res.m5 = residuals(m5_114.landing);  
residual.analysis(model = m5_114.landing)

m5_114.landingCSS = Arima(log(NMFS_Landings.ts),order=c(1,1,4),seasonal=list(order=c(1,1,1), period=12),method = "CSS")
coeftest(m5_114.landingCSS)
res.m5 = residuals(m5_114.landingCSS);  
residual.analysis(model = m5_114.landingCSS)

m5_114.landingCSSML = Arima(log(NMFS_Landings.ts),order=c(1,1,4),seasonal=list(order=c(1,1,1), period=12),method = "CSS-ML")
coeftest(m5_114.landingCSSML)
res.m5 = residuals(m5_114.landingCSSML);  
residual.analysis(model = m5_114.landingCSSML)

# SARIMA(1,1,5)x(1,1,1)_12
m5_115.landing = Arima(log(NMFS_Landings.ts),order=c(1,1,5),seasonal=list(order=c(1,1,1), period=12),method = "ML")
coeftest(m5_115.landing)
res.m5 = residuals(m5_115.landing);  
residual.analysis(model = m5_115.landing)

m5_115.landingCSS = Arima(log(NMFS_Landings.ts),order=c(1,1,5),seasonal=list(order=c(1,1,1), period=12),method = "CSS")
coeftest(m5_115.landingCSS)
res.m5 = residuals(m5_115.landingCSS);  
residual.analysis(model = m5_115.landingCSS)

# SARIMA(6,1,4)x(1,1,1)_12
m5_614.landing = Arima(log(NMFS_Landings.ts),order=c(6,1,4),seasonal=list(order=c(1,1,1), period=12),method = "ML")
coeftest(m5_614.landing)
residual.analysis(model = m5_614.landing)

m5_614.landingCSS = Arima(log(NMFS_Landings.ts),order=c(6,1,4),seasonal=list(order=c(1,1,1), period=12),method = "CSS")
coeftest(m5_614.landingCSS)
residual.analysis(model = m5_614.landingCSS)

m5_614.landingCSSML = Arima(log(NMFS_Landings.ts),order=c(6,1,4),seasonal=list(order=c(1,1,1), period=12),method = "CSS-ML")
coeftest(m5_614.landingCSSML)
residual.analysis(model = m5_614.landingCSSML)

# SARIMA(10,1,2)x(1,1,1)_12
m5_1012.landing = Arima(log(NMFS_Landings.ts),order=c(10,1,2),seasonal=list(order=c(1,1,1), period=12),method = "ML")
coeftest(m5_1012.landing)
residual.analysis(model = m5_1012.landing)

m5_1012.landingCSS = Arima(log(NMFS_Landings.ts),order=c(10,1,2),seasonal=list(order=c(1,1,1), period=12),method = "CSS")
coeftest(m5_1012.landingCSS)
residual.analysis(model = m5_1012.landingCSS)

# SARIMA(4,1,2)x(1,1,1)_12
m5_412.landing = Arima(log(NMFS_Landings.ts),order=c(4,1,2),seasonal=list(order=c(1,1,1), period=12),method = "ML")
coeftest(m5_412.landing)
residual.analysis(model = m5_412.landing)

m5_412.landingCSS = Arima(log(NMFS_Landings.ts),order=c(4,1,2),seasonal=list(order=c(1,1,1), period=12),method = "CSS")
coeftest(m5_412.landingCSS)
residual.analysis(model = m5_412.landingCSS)

# SARIMA(1,1,2)x(1,1,1)_12
m5_112.landing = Arima(log(NMFS_Landings.ts),order=c(1,1,2),seasonal=list(order=c(1,1,1), period=12),method = "ML")
coeftest(m5_112.landing)
residual.analysis(model = m5_112.landing)

m5_112.landingCSS = Arima(log(NMFS_Landings.ts),order=c(1,1,2),seasonal=list(order=c(1,1,1), period=12),method = "CSS")
coeftest(m5_112.landingCSS)
residual.analysis(model = m5_112.landingCSS)


sc.AIC = AIC(m5_014.landing, m5_015.landing, m5_114.landing, 
             m5_115.landing, m5_1012.landing, m5_412.landing, m5_112.landing)

sc.BIC = BIC(m5_014.landing, m5_015.landing,  m5_114.landing, 
    m5_115.landing, m5_1012.landing, m5_412.landing, m5_112.landing)

sort.score(sc.AIC, score = "aic")
sort.score(sc.BIC, score = "bic")


Sm5_014.landing <- accuracy(m5_014.landing)[1:7]
Sm5_015.landing <- accuracy(m5_015.landing)[1:7]
Sm5_614.landing <- accuracy(m5_614.landing)[1:7]
Sm5_114.landing <- accuracy(m5_114.landing)[1:7]
Sm5_115.landing <- accuracy(m5_115.landing)[1:7]
Sm5_1012.landing <- accuracy(m5_1012.landing)[1:7]
Sm5_412.landing <- accuracy(m5_412.landing)[1:7]
Sm5_112.landing <- accuracy(m5_112.landing)[1:7]
df.Smodels <- data.frame(
  rbind(Sm5_014.landing, Sm5_015.landing, Sm5_114.landing, 
        Sm5_115.landing, Sm5_1012.landing, Sm5_412.landing, Sm5_112.landing)
)
colnames(df.Smodels) <- c("ME", "RMSE", "MAE", "MPE", "MAPE", 
                          "MASE", "ACF1")
rownames(df.Smodels) <- c("SARIMA(0,1,4)x(1,1,1)_12", "SARIMA(0,1,5)x(1,1,1)_12",  
                          "SARIMA(1,1,4)x(1,1,1)_12", "SARIMA(1,1,5)x(1,1,1)_12", "SARIMA(10,1,2)x(1,1,1)_12",
                          "SARIMA(4,1,2)x(1,1,1)_12", "SARIMA(1,1,2)x(1,1,1)_12")
round(df.Smodels,  digits = 3)

# Best model is SARIMA(10,1,2)x(1,1,1)_12 according to AIC and error measures
# SARIMA(1,1,2)x(1,1,1)_12 according to BIC. But SARIMA(1,1,2)x(1,1,1)_12 model has lots of significant autocorrelations in residuals.
# So SARIMA(10,1,2)x(1,1,1)_12 is the best model.
# Over-parameterisations for SARIMA(10,1,2)x(1,1,1)_12: SARIMA(11,1,2)x(1,1,1)_12 and SARIMA(10,1,3)x(1,1,1)_12

m5_1112.landing = Arima(log(NMFS_Landings.ts),order=c(11,1,2),seasonal=list(order=c(1,1,1), period=12))
coeftest(m5_1112.landing)
residual.analysis(model = m5_1112.landing)
# AR(11) turned significant.

m5_1013.landing = Arima(log(NMFS_Landings.ts),order=c(10,1,3),seasonal=list(order=c(1,1,1), period=12))
coeftest(m5_1013.landing)
residual.analysis(model = m5_1013.landing)
# MA(3) is insignificant

sc.AIC = AIC(m5_014.landing, m5_015.landing, m5_614.landing, m5_114.landing, 
             m5_115.landing, m5_1012.landing, m5_412.landing, m5_112.landing, m5_1112.landing)
sort.score(sc.AIC, score = "aic")

# SARIMA(10,1,2)x(1,1,1)_12 is better than SARIMA(11,1,2)x(1,1,1)_12 by AIC.

# Forecasting
m5_1012.landingA = Arima(NMFS_Landings.ts , order=c(10,1,2),seasonal=list(order=c(1,1,1), period=12), 
                         lambda = 0, method = "CSS")
# Notice that I use lambda = 0 and send NMFS_Landings.ts instead of log(NMFS_Landings.ts) to get Arima() function to do the transformation.
# This way, I will get the forecasts in the original scale.
preds1 = forecast(m5_1012.landingA, lambda = 0, h = 48)
preds1
plot(preds1)

# Notice the difference! In the above code I use lambda = 0 and send the raw series to Arima() function. 
# The resulting forecasts are in the original scale.
# In the below one, I send the log transformed series to Arima() function. The resulting forecasts are in the log scale! 
# Be careful to get the forecasts in the original scale!

m5_1012.landingA2 = Arima(log(NMFS_Landings.ts),order=c(10,1,2),seasonal=list(order=c(1,1,1), period=12), 
                        method = "CSS")
preds2 = forecast(m5_1012.landingA2, h = 48)
preds2
plot(preds2)

# ------------- Taking the first 8 years of data as the first part. -----------------

par(mfrow=c(1,1))
NMFS_Landings.p1 = ts(NMFS_Landings.ts[1:96],start=c(1991,1), frequency=12)

plot(NMFS_Landings.p1,ylab='Landings in metric tons',xlab='Year',type='o', main = "Time series plot of monthly landings in metric tons for the first 8 years.")

par(mfrow=c(1,2))
seasonal_acf(NMFS_Landings.p1, lag.max = 48,main="The sample ACF of landings series")
seasonal_pacf(NMFS_Landings.p1, lag.max = 48,main="The sample PACF of landings series")


# Start with the first seasonal difference
m1.landing = Arima(log(NMFS_Landings.p1),order=c(0,0,0),seasonal=list(order=c(0,1,0), period=12))
res.m1 = residuals(m1.landing);  
par(mfrow=c(1,1))
plot(res.m1,xlab='Time',ylab='Residuals',main="Time series plot of the residuals")
par(mfrow=c(1,2))
seasonal_acf(res.m1, lag.max = 48, main = "The sample ACF of the residuals")
seasonal_pacf(res.m1, lag.max = 48, main = "The sample PACF of the residuals")
# Q = 1 from ACF.
# And there is no sign of an ordinary trend as well. So we will go on with specification of ARMA component. 

m2.landing = Arima(log(NMFS_Landings.p1),order=c(0,0,0),
                   seasonal=list(order=c(0,1,1), period=12))
res.m2 = residuals(m2.landing);  
par(mfrow=c(1,1))
plot(res.m2,xlab='Time',ylab='Residuals',main="Time series plot of the residuals")
par(mfrow=c(1,2))
seasonal_acf(res.m2, lag.max = 36, main = "The sample ACF of the residuals")
seasonal_pacf(res.m2, lag.max = 36, main = "The sample PACF of the residuals")

adf.test(res.m2)
pp.test(res.m2)

# We have one significant correlation in ACF and one significant one in PACF in the section before the periods.
# So we can think of an ARMA(1,1) model.
m3.landing = Arima(log(NMFS_Landings.p1),order=c(1,0,1),
                   seasonal=list(order=c(0,1,1), period=12))
res.m3 = residuals(m3.landing);  
par(mfrow=c(1,1))
plot(res.m3,xlab='Time',ylab='Residuals',main="Time series plot of the residuals")
par(mfrow=c(1,2))
seasonal_acf(res.m3, lag.max = 36, main = "The sample ACF of the residuals")
seasonal_pacf(res.m3, lag.max = 36, main = "The sample PACF of the residuals")
# We got the white noise series with ARMA(1,1) component.


eacf(res.m2) # Be careful - residuals of m2 is being used for EACF since there is leftover signal in them.

# The tentative models are specified as 
# SARIMA(1,0,1)x(0,1,1)_12
# From the EACF, we will include 
# SARIMA(2,0,2)x(0,1,1)_12, SARIMA(2,0,3)x(0,1,1)_12 and SARIMA(3,0,3)x(0,1,1)_12 

par(mfrow=c(1,1))
bic_table = armasubsets(y=res.m2,nar=5,nma=5,y.name='p',ar.method='ols')
plot(bic_table)
# SARIMA(0,1,4)x(1,1,1)_12
# SARIMA(0,1,2)x(1,1,1)_12

#{SARIMA(1,0,1)x(0,1,1)_12, SARIMA(2,0,2)x(0,1,1)_12, 
# SARIMA(2,0,3)x(0,1,1)_12, SARIMA(3,0,3)x(0,1,1)_12, 
# SARIMA(0,1,4)x(1,1,1)_12, SARIMA(0,1,2)x(1,1,1)_12}

# SARIMA(1,0,1)x(0,1,1)_12
m3_101.landing = Arima(log(NMFS_Landings.p1),order=c(1,0,1),seasonal=list(order=c(0,1,1), period=12),method = "ML")
coeftest(m3_101.landing)
residual.analysis(model = m3_101.landing)

m3_101.landingCSS = Arima(log(NMFS_Landings.p1),order=c(1,0,1),seasonal=list(order=c(0,1,1), period=12),method = "CSS")
coeftest(m3_101.landingCSS)
residual.analysis(model = m3_101.landingCSS)

# SARIMA(2,0,2)x(0,1,1)_12
m3_202.landing = Arima(log(NMFS_Landings.p1),order=c(2,0,2),seasonal=list(order=c(0,1,1), period=12),method = "ML")
coeftest(m3_202.landing)
residual.analysis(model = m3_202.landing)

m3_202.landingCSS = Arima(log(NMFS_Landings.p1),order=c(2,0,2),seasonal=list(order=c(0,1,1), period=12),method = "CSS")
coeftest(m3_202.landingCSS)
residual.analysis(model = m3_202.landingCSS)

# SARIMA(2,0,3)x(0,1,1)_12
m3_203.landing = Arima(log(NMFS_Landings.p1),order=c(2,0,3),seasonal=list(order=c(0,1,1), period=12),method = "ML")
coeftest(m3_203.landing)
residual.analysis(model = m3_203.landing)

m3_203.landingCSS = Arima(log(NMFS_Landings.p1),order=c(2,0,3),seasonal=list(order=c(0,1,1), period=12),method = "CSS")
coeftest(m3_203.landingCSS)
residual.analysis(model = m3_203.landingCSS)

# SARIMA(3,0,3)x(0,1,1)_12 
m3_303.landing = Arima(log(NMFS_Landings.p1),order=c(3,0,3),seasonal=list(order=c(0,1,1), period=12),method = "ML")
coeftest(m3_303.landing)
residual.analysis(model = m3_303.landing)

m3_303.landingCSS = Arima(log(NMFS_Landings.p1),order=c(3,0,3),seasonal=list(order=c(0,1,1), period=12),method = "CSS")
coeftest(m3_303.landingCSS)
residual.analysis(model = m3_303.landingCSS)

# SARIMA(0,0,4)x(0,1,1)_12 
m3_004.landing = Arima(log(NMFS_Landings.p1),order=c(0,0,4),seasonal=list(order=c(0,1,1), period=12),method = "ML")
coeftest(m3_004.landing)
residual.analysis(model = m3_004.landing)

m3_004.landingCSS = Arima(log(NMFS_Landings.p1),order=c(0,0,4),seasonal=list(order=c(0,1,1), period=12),method = "CSS")
coeftest(m3_004.landingCSS)
residual.analysis(model = m3_004.landingCSS)

# SARIMA(0,0,2)x(0,1,1)_12 
m3_002.landing = Arima(log(NMFS_Landings.p1),order=c(0,0,2),seasonal=list(order=c(0,1,1), period=12),method = "ML")
coeftest(m3_002.landing)
residual.analysis(model = m3_002.landing)

m3_002.landingCSS = Arima(log(NMFS_Landings.p1),order=c(0,0,2),seasonal=list(order=c(0,1,1), period=12),method = "CSS")
coeftest(m3_002.landingCSS)
residual.analysis(model = m3_002.landingCSS)


sc.AIC=AIC(m3_101.landing, m3_202.landing, 
           m3_203.landing, m3_303.landing, m3_004.landing, m3_002.landing)
sc.BIC=BIC(m3_101.landing, m3_202.landing, 
           m3_203.landing, m3_303.landing, m3_004.landing, m3_002.landing)

sort.score(sc.AIC, score = "aic")
sort.score(sc.BIC, score = "bic")

Sm3_101.landing <- accuracy(m3_101.landing)[1:7]
Sm3_202.landing <- accuracy(m3_202.landing)[1:7]
Sm3_203.landing <- accuracy(m3_203.landing)[1:7]
Sm3_303.landing <- accuracy(m3_303.landing)[1:7]
Sm3_004.landing <- accuracy(m3_004.landing)[1:7]
Sm3_002.landing <- accuracy(m3_002.landing)[1:7]
df.Smodels <- data.frame(
  rbind(Sm3_101.landing, Sm3_202.landing, Sm3_203.landing,  
        Sm3_303.landing, Sm3_004.landing, Sm3_002.landing)
)
colnames(df.Smodels) <- c("ME", "RMSE", "MAE", "MPE", "MAPE", 
                          "MASE", "ACF1")
rownames(df.Smodels) <- c("SARIMA(1,0,1)x(0,1,1)_12", "SARIMA(2,0,2)x(0,1,1)_12", "SARIMA(2,0,3)x(0,1,1)_12", 
                          "SARIMA(3,0,3)x(0,1,1)_12", "SARIMA(0,0,4)x(0,1,1)_12", "SARIMA(0,0,2)x(0,1,1)_12")
round(df.Smodels,  digits = 3)

# Best model SARIMA(0,0,2)x(0,1,1)_12
# Over-parameterizations: SARIMA(1,0,2)x(0,1,1)_12 and 
# SARIMA(0,0,3)x(0,1,1)_12

m3_102.landing = Arima(log(NMFS_Landings.p1),order=c(1,0,2),seasonal=list(order=c(0,1,1), period=12),method = "ML")
coeftest(m3_102.landing)
residual.analysis(model = m3_102.landing)


m3_003.landing = Arima(log(NMFS_Landings.p1),order=c(0,0,3),seasonal=list(order=c(0,1,1), period=12),method = "ML")
coeftest(m3_003.landing)
residual.analysis(model = m3_003.landing)

m2.landingA = Arima(NMFS_Landings.p1,order=c(0,0,2),seasonal=list(order=c(0,1,1), period=12), 
                    lambda = 0, method = "ML")
future = forecast(m2.landingA, lambda = 0, h = 48)
future
plot(future) # Forecasts show what would had happened if there was no policy change.


# ------------ Taking the second period of data as the first part. -----------

NMFS_Landings.p2 =  ts(NMFS_Landings.ts[97:300],start=c(1999,1), frequency=12)

plot(NMFS_Landings.p2,ylab='Landings in metric tons',xlab='Year',type='o', main = "Time series plot of monthly landings in metric tons for the first 8 years.")
# Notice the changing variance in this plot. So, a log transformation would be useful.

par(mfrow=c(1,2))
seasonal_acf(NMFS_Landings.p2,lag.max = 48, main="The sample ACF of landings series")
seasonal_pacf(NMFS_Landings.p2,lag.max = 48, main="The sample PACF of landings series")



# Start with the first seasonal difference
m1.landing = Arima(log(NMFS_Landings.p2),order=c(0,0,0),
                   seasonal=list(order=c(0,1,0), period=12))
res.m1 = residuals(m1.landing);  
par(mfrow=c(1,1))
plot(res.m1,xlab='Time',ylab='Residuals',main="Time series plot of the residuals")
par(mfrow=c(1,2))
seasonal_acf(res.m1, lag.max = 60, main = "The sample ACF of the residuals")
seasonal_pacf(res.m1, lag.max = 60, main = "The sample PACF of the residuals")
# The highly significant first seasonal lag in ACF and slowly decreasing seasonal lags in PACF indicates a SMA(1) model.

m2.landing1 = Arima(log(NMFS_Landings.p2),order=c(0,0,0),
                   seasonal=list(order=c(0,1,1), period=12))
res.m21 = residuals(m2.landing1);  
par(mfrow=c(1,1))
plot(res.m2,xlab='Time',ylab='Residuals',main="Time series plot of the residuals")
par(mfrow=c(1,2))
seasonal_acf(res.m2, lag.max = 36, main = "The sample ACF of the residuals")
seasonal_pacf(res.m2, lag.max = 36, main = "The sample PACF of the residuals")
# There is no seasonal effect left in the residuals. From the time series plot and ACF-PACF pair, we do not observe an ordinary trend. 
# We consider 1 significant ordinary lags in ACF and 1 significant ans 1 slightly insignificant ordinary lags in PACF for model specification.
# Tentative models are 
# SARIMA(2,0,1)x(0,1,1)_12

adf.test(res.m2)
pp.test(res.m2)

m3.landing = Arima(log(NMFS_Landings.p2),order=c(2,0,1),
                   seasonal=list(order=c(0,1,1), period=12))
res.m3 = residuals(m3.landing);  
par(mfrow=c(1,1))
plot(res.m3,xlab='Time',ylab='Residuals',main="Time series plot of the residuals")
par(mfrow=c(1,2))
seasonal_acf(res.m3, lag.max = 36, main = "The sample ACF of the residuals")
seasonal_pacf(res.m3, lag.max = 36, main = "The sample PACF of the residuals")

# There are still significant lags in ACF and PACF. So, I'll increase q = 2 to see if there is any difference.

m3.landing = Arima(log(NMFS_Landings.p2),order=c(2,0,2),
                   seasonal=list(order=c(0,1,1), period=12))
res.m3 = residuals(m3.landing);  
par(mfrow=c(1,1))
plot(res.m3,xlab='Time',ylab='Residuals',main="Time series plot of the residuals")
par(mfrow=c(1,2))
seasonal_acf(res.m3, lag.max = 36, main = "The sample ACF of the residuals")
seasonal_pacf(res.m3, lag.max = 36, main = "The sample PACF of the residuals")

# no significant lags left. So, one of the promising models is SARIMA(2,0,2)x(0,1,1)_12.
# The tentative models are specified as 
# SARIMA(2,0,2)x(0,1,1)_12

eacf(res.m21)

# From the EACF, we will include 
# SARIMA(1,0,0)x(0,1,1)_12, SARIMA(1,0,1)x(0,1,1)_12 and 
# SARIMA(2,0,1)x(0,1,1)_12 

par(mfrow=c(1,1))
bic_table = armasubsets(y=res.m21,nar=5,nma=5,y.name='p',ar.method='ols')
plot(bic_table)
# SARIMA(1,0,0)x(1,1,1)_12
# SARIMA(4,0,0)x(1,1,1)_12


# SARIMA(2,0,2)x(0,1,1)_12
m3_202.landing = Arima(log(NMFS_Landings.p2),order=c(2,0,2),seasonal=list(order=c(0,1,1), period=12),method = "ML")
coeftest(m3_202.landing)
residual.analysis(model = m3_202.landing)

m3_202.landingCSS = Arima(log(NMFS_Landings.p2),order=c(2,0,2),seasonal=list(order=c(0,1,1), period=12),method = "CSS")
coeftest(m3_202.landingCSS)
residual.analysis(model = m3_202.landingCSS)

# SARIMA(1,0,0)x(0,1,1)_12
m3_100.landing = Arima(log(NMFS_Landings.p2),order=c(1,0,0),seasonal=list(order=c(0,1,1), period=12),method = "ML")
coeftest(m3_100.landing)
residual.analysis(model = m3_100.landing)

m3_100.landingCSS = Arima(log(NMFS_Landings.p2),order=c(1,0,0),seasonal=list(order=c(0,1,1), period=12),method = "CSS")
coeftest(m3_100.landingCSS)
residual.analysis(model = m3_100.landingCSS)

# SARIMA(1,0,1)x(0,1,1)_12
m3_101.landing = Arima(log(NMFS_Landings.p2),order=c(1,0,1),seasonal=list(order=c(0,1,1), period=12),method = "ML")
coeftest(m3_101.landing)
residual.analysis(model = m3_101.landing)

m3_101.landingCSS = Arima(log(NMFS_Landings.p2),order=c(1,0,1),seasonal=list(order=c(0,1,1), period=12),method = "CSS")
coeftest(m3_101.landingCSS)
residual.analysis(model = m3_101.landingCSS)

# SARIMA(2,0,1)x(0,1,1)_12
m3_201.landing = Arima(log(NMFS_Landings.p2),order=c(2,0,1),seasonal=list(order=c(0,1,1), period=12),method = "ML")
coeftest(m3_201.landing)
residual.analysis(model = m3_201.landing)

m3_201.landingCSS = Arima(log(NMFS_Landings.p2),order=c(2,0,1),seasonal=list(order=c(0,1,1), period=12),method = "CSS")
coeftest(m3_201.landingCSS)
residual.analysis(model = m3_201.landingCSS)

# SARIMA(4,0,0)x(0,1,1)_12
m3_400.landing = Arima(log(NMFS_Landings.p2),order=c(4,0,0),seasonal=list(order=c(0,1,1), period=12),method = "ML")
coeftest(m3_400.landing)
residual.analysis(model = m3_400.landing)

m3_400.landingCSS = Arima(log(NMFS_Landings.p2),order=c(4,0,0),seasonal=list(order=c(0,1,1), period=12),method = "CSS")
coeftest(m3_400.landingCSS)
residual.analysis(model = m3_400.landingCSS)

sc.AIC=AIC(m3_202.landing, m3_101.landing, m3_201.landing, 
           m3_100.landing, m3_400.landing)
sc.BIC=BIC(m3_202.landing, m3_101.landing, m3_201.landing, 
           m3_100.landing, m3_400.landing)

sort.score(sc.AIC, score = "aic")
sort.score(sc.BIC, score = "bic")


Sm3_101.landing <- accuracy(m3_101.landing)[1:7]
Sm3_202.landing <- accuracy(m3_202.landing)[1:7]
Sm3_201.landing <- accuracy(m3_201.landing)[1:7]
Sm3_100.landing <- accuracy(m3_100.landing)[1:7]
Sm3_400.landing <- accuracy(m3_400.landing)[1:7]
df.Smodels <- data.frame(
  rbind(Sm3_101.landing, Sm3_202.landing, Sm3_201.landing,  
        Sm3_100.landing, Sm3_400.landing)
)
colnames(df.Smodels) <- c("ME", "RMSE", "MAE", "MPE", "MAPE", 
                          "MASE", "ACF1")
rownames(df.Smodels) <- c("SARIMA(1,0,1)x(0,1,1)_12", "SARIMA(2,0,2)x(0,1,1)_12", "SARIMA(2,0,1)x(0,1,1)_12", 
                          "SARIMA(1,0,0)x(0,1,1)_12", "SARIMA(4,0,0)x(0,1,1)_12")
round(df.Smodels,  digits = 3)

# Best model SARIMA(1,0,0)x(0,1,1)_12
# Over-parameterizations: SARIMA(2,0,0)x(0,1,1)_12 and 
# SARIMA(1,0,1)x(0,1,1)_12

m3_200.landing = Arima(log(NMFS_Landings.p2),order=c(2,0,0),seasonal=list(order=c(0,1,1), period=12),method = "ML")
coeftest(m3_200.landing)
residual.analysis(model = m3_200.landing)

m3_101.landing = Arima(log(NMFS_Landings.p2),order=c(1,0,1),seasonal=list(order=c(0,1,1), period=12),method = "ML")
coeftest(m3_101.landing)
residual.analysis(model = m3_101.landing)


m3.landing = Arima(NMFS_Landings.p2,order=c(1,0,0),
                   seasonal=list(order=c(0,1,1), period=12), 
                   lambda = 0)
future = forecast(m3.landing, lambda = 0, h = 48)
future
plot(future)

# Compare it with the forecasts of the full analysis
preds1
plot(preds1)
