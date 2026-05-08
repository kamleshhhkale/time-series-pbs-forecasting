
library(TSA)
library(fUnitRoots)
library(lmtest)
library(tseries)
library(forecast)
library(FitAR)
library(Hmisc)
# Download the latest version 1.94 of FitAR from archive:https://cran.r-project.org/src/contrib/Archive/FitAR/
# Run install.packages("bestglm")
# And install it from the source with the following:
# install.packages("FOLDER WHERE YOU DOWNLOADED THE PACKAGE SOURCE FILE", repos = NULL, type = "source")

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
  acf(res.model,main="ACF of standardised residuals")
  print(shapiro.test(res.model))
  k=0
  LBQPlot(res.model, lag.max = 30, StartLag = k + 1, k = 0, SquaredQ = FALSE)
  par(mfrow=c(1,1))
}

setwd("~/Documents/MATH1318_TimeSeries/tasks/Task7")

#--- Task 1 ---

unemployment <- read.table("unemployment.csv", 
                           quote="\"", comment.char="")
head(unemployment)

unemployment = ts(unemployment,start = 1978)
class(unemployment)

plot(unemployment,type='o',ylab= "Unemployment",main='Time series plot of unemployment series.')

acf(unemployment, main = "ACF of unemployment series.")
pacf(unemployment, main = "PACF of unemployment series.")

adf.test(unemployment)

qqnorm(unemployment, ylab="Unemployment", xlab="Normal Scores")
qqline(unemployment)
shapiro.test(unemployment) 

BC = BoxCox.ar(unemployment)
BC$ci
lambda = BC$lambda[which(max(BC$loglike) == BC$loglike)]
lambda
unemploymentBC = (unemployment^lambda-1)/lambda

qqnorm(unemploymentBC, ylab="BC - Unemployment", xlab="Normal Scores")
qqline(unemploymentBC)
shapiro.test(unemploymentBC)

plot(unemploymentBC,type='o',main='Time series plot of BC transformation of unemployment series.')

# Although Box-Cox transformation does not improve normality, I'll go on with that to demonstrate how Arima() handles the transformation.

diff.unemploymentBC = diff(unemploymentBC, differences = 1) 
plot(diff.unemploymentBC,type='o',ylab= "Unemployment",main='The first difference of unemployment series.')

adf.test(diff.unemploymentBC)

pp.test(diff.unemploymentBC)

par(mfrow=c(1,2))
acf(diff.unemploymentBC, main = "ACF of the first difference 
    of unemployment series.")
pacf(diff.unemploymentBC, main = "PACF of the first difference 
     of unemployment series.")
par(mfrow=c(1,1))

# {ARIMA(1,1,1)} 

eacf(diff.unemploymentBC, ar.max = 5, ma.max = 5)
# {ARIMA(1,1,1), ARIMA(0,1,1), ARIMA(0,1,2), ARIMA(1,1,2)}

bic_ms = armasubsets(y=diff.unemploymentBC,nar=5,nma=5,y.name='p',ar.method='ols')
plot(bic_ms)
# {ARIMA(1,1,1), ARIMA(0,1,1), ARIMA(0,1,2), ARIMA(1,1,2), 
#  ARIMA(5,1,1), ARIMA(5,1,0)}


# ARIMA(1,1,1)
model_111_css = Arima(unemployment,order=c(1,1,1),method='CSS', lambda = 1.6 )
# Note that I use raw series and pass the BC transformation on to Arima() using , lambda = 1.6!

coeftest(model_111_css)

model_111_ml = Arima(unemployment,order=c(1,1,1),method='ML', lambda = 1.6 )
coeftest(model_111_ml)

# ARIMA(0,1,1)
model_011_css = Arima(unemployment,order=c(0,1,1),method='CSS', lambda = 1.6 )
coeftest(model_011_css)

model_011_ml = Arima(unemployment,order=c(0,1,1),method='ML', lambda = 1.6 )
coeftest(model_011_ml)

# ARIMA(0,1,2)
model_012_css = Arima(unemployment,order=c(0,1,2),method='CSS', lambda = 1.6 )
coeftest(model_012_css)

model_012_ml = Arima(unemployment,order=c(0,1,2),method='ML', lambda = 1.6 )
coeftest(model_012_ml)

# ARIMA(1,1,2) 
model_112_css = Arima(unemployment,order=c(1,1,2),method='CSS', lambda = 1.6 )
coeftest(model_112_css)

model_112_ml = Arima(unemployment,order=c(1,1,2),method='ML', lambda = 1.6 )
coeftest(model_112_ml)

# ARIMA(5,1,1)
model_511_css = Arima(unemployment,order=c(5,1,1),method='CSS', lambda = 1.6 )
coeftest(model_511_css)

model_511_ml = Arima(unemployment,order=c(5,1,1),method='ML', lambda = 1.6 )
coeftest(model_511_ml)

# ARIMA(5,1,0)
model_510_css = Arima(unemployment,order=c(5,1,0),method='CSS', lambda = 1.6 )
coeftest(model_510_css)

model_510_ml = Arima(unemployment,order=c(5,1,0),method='ML', lambda = 1.6 )
coeftest(model_510_ml)

# AIC and BIC values
sort.score(AIC(model_111_ml,model_011_ml,model_012_ml,model_112_ml,model_511_ml,model_510_ml), score = "aic")
sort.score(BIC(model_111_ml,model_011_ml,model_012_ml,model_112_ml,model_511_ml,model_510_ml), score = "bic" )

# The ARIMA(1,1,2) model is the best one according to both AIC and BIC

Smodel_111_ml <- accuracy(model_111_ml)[1:7]
Smodel_011_ml <- accuracy(model_011_ml)[1:7]
Smodel_012_ml <- accuracy(model_012_ml)[1:7]
Smodel_112_ml <- accuracy(model_112_ml)[1:7]
Smodel_511_ml <- accuracy(model_511_ml)[1:7]
Smodel_510_ml <- accuracy(model_510_ml)[1:7]
df.Smodels <- data.frame(
  rbind(Smodel_111_ml,Smodel_011_ml,Smodel_012_ml,
        Smodel_112_ml,Smodel_511_ml,Smodel_510_ml)
)
colnames(df.Smodels) <- c("ME", "RMSE", "MAE", "MPE", "MAPE", 
                          "MASE", "ACF1")
rownames(df.Smodels) <- c("ARIMA(1,1,1)", "ARIMA(0,1,1)", "ARIMA(0,1,2)", 
                          "ARIMA(1,1,2)", "ARIMA(5,1,1)", "ARIMA(5,1,0)")
round(df.Smodels,  digits = 3)

# AIC/BIC and most error measures support the ARIMA(1,1,2) model as well.

# Overfitting: To further assess the selected model ARIMA(1,1,2) by overfitting
# ARIMA(2,1,2) and ARIMA(1,1,3)

# ARIMA(2,1,2)
model_212_css = Arima(unemployment,order=c(2,1,2),method='CSS', lambda = 1.6 )
coeftest(model_212_css) # Gave us a significant AR(2) coefficient.

model_212_ml = Arima(unemployment,order=c(2,1,2),method='ML', lambda = 1.6 )
coeftest(model_212_ml)

model_212_CSSml = Arima(unemployment,order=c(2,1,2),method='CSS-ML', lambda = 1.6 )
coeftest(model_212_CSSml)

# ARIMA(1,1,3)
model_113_css = Arima(unemployment,order=c(1,1,3),method='CSS', lambda = 1.6 )
coeftest(model_113_css)

model_113_ml = Arima(unemployment,order=c(1,1,3),method='ML', lambda = 1.6 )
coeftest(model_113_ml)

# Now I will consider ARIMA(2,1,2) among the others.

sort.score(AIC(model_111_ml,model_011_ml,model_012_ml,model_112_ml,model_212_ml,model_511_ml,model_510_ml), score = "aic")
sort.score(BIC(model_111_ml,model_011_ml,model_012_ml,model_112_ml,model_212_ml,model_511_ml,model_510_ml), score = "bic" )


residual.analysis(model = model_112_ml) # Best model
residual.analysis(model = model_212_ml) # New model
residual.analysis(model = model_212_css) # New model


frcCSS = forecast::forecast(model_112_css, h=5 , lambda = 1.6 )
# Note that once you have "lambda = 1.6" in your model fitted with Arima(), you don't even need to include "lambda = 1.6" in forecast().
plot(frcCSS) 
lines(Lag(fitted(model_112_css),-1), col= "blue")
legend("topleft", lty=1, pch=1, col=c("blue","black"), text.width = 2, c("Data", "Fitted "))
frcCSS

frcML = forecast::forecast(model_112_ml,h=5)
plot(frcML) 
lines(Lag(fitted(model_112_ml),-1), col= "blue")
legend("topleft", lty=1, pch=1, col=c("blue","black"), text.width = 2, c("Data", "Fitted "))
frcML



#--- Task 2 ---
# arima.data1 = arima.sim(list(order = c(2,1,3),ar = c(-0.89,-0.26), ma=c(0.71,0.31,0.2)), n = 200) 
# This is how I have simulated this series. So you know the true values of parameters and the true orders. 
# You can compare your final results with the true values to see how good is your model.

arima.data1 <- read.csv("data.sim.csv", header=FALSE)
arima.data1 = ts(arima.data1,start = 1794)
class(arima.data1)

plot(arima.data1,type='o',main='Time series plots of simulated series.')

par(mfrow=c(1,2))
acf(arima.data1, main = "ACF of the simulated series.")
pacf(arima.data1, main = "PACF of the simulated series.")
par(mfrow=c(1,1))

adf.test(arima.data1)

qqnorm(arima.data1, ylab="Simulated series", xlab="Normal Scores")
qqline(arima.data1)
shapiro.test(arima.data1)

BC = BoxCox.ar(arima.data1+abs(min(arima.data1))+0.1) # There are negative values in the series.
BC$ci
lambda = BC$lambda[which(max(BC$loglike) == BC$loglike)]
lambda
arima.data1BC = ((arima.data1+abs(min(arima.data1))+0.1)^lambda-1)/lambda

qqnorm(arima.data1BC, ylab="BC - Simulated series", xlab="Normal Scores")
qqline(arima.data1BC)
shapiro.test(arima.data1BC)

plot(arima.data1BC,type='o',main='Time series plots of BC transformed simulated series.')
# Almost no improvement.

diff.arima.data1= diff(arima.data1,difference=1)
plot(diff.arima.data1,type='o',main='Time sereis plots of the first difference of the simulated series.')

adf.test(diff.arima.data1)

par(mfrow=c(1,2))
acf(diff.arima.data1, main = "ACF of the first difference 
    of the simulated series.")
pacf(diff.arima.data1, main = "PACF of the first difference 
     of the simulated series.")
par(mfrow=c(1,1))
#{ARIMA(2,1,0), ARIMA(2,1,4)}

eacf(diff.arima.data1)
#{ARIMA(2,1,0), ARIMA(2,1,4), 
# ARIMA(1,1,2), ARIMA(1,1,3), ARIMA(2,1,2), ARIMA(2,1,3)}

res = armasubsets(y=diff.arima.data1,nar=7,nma=7,y.name='p',ar.method='ols')
plot(res)
#{ARIMA(2,1,0), ARIMA(2,1,4), 
# ARIMA(1,1,2), ARIMA(1,1,3), ARIMA(2,1,2), ARIMA(2,1,3), 
# ARIMA(1,1,4), ARIMA(3,1,4)}


# ARIMA(2,1,0)
model_210_css = Arima(arima.data1,order=c(2,1,0),method='CSS')
coeftest(model_210_css)

model_210_ml = Arima(arima.data1,order=c(2,1,0),method='ML')
coeftest(model_210_ml)

# ARIMA(2,1,4)
model_214_css = Arima(arima.data1,order=c(2,1,4),method='CSS')
coeftest(model_214_css)

model_214_ml = Arima(arima.data1,order=c(2,1,4),method='ML')
coeftest(model_214_ml)

# ARIMA(1,1,2)
model_112_css = Arima(arima.data1,order=c(1,1,2),method='CSS')
coeftest(model_112_css)

model_112_ml = Arima(arima.data1,order=c(1,1,2),method='ML')
coeftest(model_112_ml)

# ARIMA(1,1,3)
model_113_css = Arima(arima.data1,order=c(1,1,3),method='CSS')
coeftest(model_113_css)

model_113_ml = Arima(arima.data1,order=c(1,1,3),method='ML')
coeftest(model_113_ml)

# ARIMA(2,1,2)
model_212_css = Arima(arima.data1,order=c(2,1,2),method='CSS')
coeftest(model_212_css)

model_212_ml = Arima(arima.data1,order=c(2,1,2),method='ML')
coeftest(model_212_ml)

model_212_CSSml = Arima(arima.data1,order=c(2,1,2),method='CSS-ML')
coeftest(model_212_CSSml)

# ARIMA(2,1,3)
model_213_css = Arima(arima.data1,order=c(2,1,3),method='CSS')
coeftest(model_213_css)

model_213_ml = Arima(arima.data1,order=c(2,1,3),method='ML')
coeftest(model_213_ml)

# ARIMA(1,1,4)
model_114_css = Arima(arima.data1,order=c(1,1,4),method='CSS')
coeftest(model_114_css)

model_114_ml = Arima(arima.data1,order=c(1,1,4),method='ML')
coeftest(model_114_ml)

# ARIMA(3,1,4)
model_314_css = Arima(arima.data1,order=c(3,1,4),method='CSS')
coeftest(model_314_css)

model_314_ml = Arima(arima.data1,order=c(3,1,4),method='ML')
coeftest(model_314_ml)


# AIC and BIC values
sort.score(AIC(model_210_ml,model_214_ml,model_213_ml,model_112_ml,model_212_ml,model_113_ml,model_114_ml,model_314_ml), score = "aic")
sort.score(BIC(model_210_ml,model_214_ml,model_213_ml,model_112_ml,model_212_ml,model_113_ml,model_114_ml,model_314_ml), score = "bic" )

# AIC favours ARIMA(2,1,3). But BIC picks ARIMA(1,1,2).

Smodel_210_ml <- accuracy(model_210_ml)[1:7]
Smodel_214_ml <- accuracy(model_214_ml)[1:7]
Smodel_213_ml <- accuracy(model_213_ml)[1:7]
Smodel_212_m <- accuracy(model_212_ml)[1:7]
Smodel_113_ml <- accuracy(model_113_ml)[1:7]
Smodel_114_ml <- accuracy(model_114_ml)[1:7]
Smodel_314_ml <- accuracy(model_314_ml)[1:7]
df.Smodels <- data.frame(
  rbind(Smodel_210_ml,Smodel_214_ml,Smodel_213_ml, Smodel_212_m,
        Smodel_113_ml,Smodel_114_ml,Smodel_314_ml)
)
colnames(df.Smodels) <- c("ME", "RMSE", "MAE", "MPE", "MAPE", 
                          "MASE", "ACF1")
rownames(df.Smodels) <- c("ARIMA(2,1,0)", "ARIMA(2,1,4)", "ARIMA(2,1,3)", 
                          "ARIMA(2,1,2)", "ARIMA(1,1,3)", "ARIMA(1,1,4)", "ARIMA(3,1,4)")
round(df.Smodels,  digits = 3)

# Error measures support the ARIMA(2,1,3) model. ARIMA(1,1,2) has slightly higher error measures.

residual.analysis(model = model_213_ml)
residual.analysis(model = model_212_ml)
par(mfrow=c(1,1))


# For overfitting, ARIMA(2,1,4) and ARIMA(3,1,3)

# ARIMA(2,1,3)
model_214_css = Arima(arima.data1,order=c(2,1,4),method='CSS')
coeftest(model_214_css) 

model_214_ml = Arima(arima.data1,order=c(2,1,4),method='ML')
coeftest(model_214_ml)

# ARIMA(3,1,3)
model_313_css = Arima(arima.data1,order=c(3,1,3),method='CSS')
coeftest(model_313_css) # Although AR(3) is significant, the rest turned insignificant.

model_313_ml = Arima(arima.data1,order=c(3,1,3),method='ML')
coeftest(model_313_ml)

library(forecast)

frc = forecast::forecast(model_213_ml,h=150) 
plot(frc)
lines(Lag(fitted(model_213_ml),-1), col= "blue")
legend("bottomleft", lty=1, pch=1, col=c("blue","black"), text.width = 11, c("Data", "Fitted "))
frc

