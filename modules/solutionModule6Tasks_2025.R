
library(TSA)
library(fUnitRoots)
library(lmtest)
library(tseries)
library(forecast)

sort.score <- function(x, score = c("bic", "aic")){
  if (score == "aic"){
    x[with(x, order(AIC)),]
  } else if (score == "bic") {
    x[with(x, order(BIC)),]
  } else {
    warning('score = "x" only accepts valid arguments ("aic","bic")')
  }
}


#--- Task 1 ---

#average weekly earnings in Novembers in Australia in 1994-2019 (27 years)
earnings <- read.table("~/Documents/MATH1318_TimeSeries/tasks/Task6/earnings.csv")
earnings = ts(earnings,start = 1994)
class(earnings)

plot(earnings,type='o',ylab = "Earnings", main='Time series plot of earnings series.')

par(mfrow=c(1,2))
acf(earnings, main="ACF of earnings series.")
pacf(earnings, main="PACF of earnings series.")
par(mfrow=c(1,1))


adf.test(earnings)

qqnorm(earnings, ylab="earnings", xlab="Normal Scores")
qqline(earnings)
shapiro.test(earnings) 

BC = BoxCox.ar(earnings) 
BC$ci
lambda <- BC$lambda[which(max(BC$loglike) == BC$loglike)]
lambda
BC.earnings = (earnings^lambda-1)/lambda

plot(BC.earnings,type='o',ylab = "Earnings", main='Time series plot of BC transformed
     earnings series.')
qqnorm(BC.earnings, ylab="BC.earnings", xlab="Normal Scores")
qqline(BC.earnings)
shapiro.test(BC.earnings) 

adf.test(BC.earnings)


diff.earnings = diff(earnings,differences = 1)
plot(diff.earnings,type='o',ylab='Earnings ', main='Time series plot of the first difference of earnings series.')

par(mfrow=c(1,2))
acf(diff.earnings, main="ACF of the first difference 
    of earnings series.")
pacf(diff.earnings, main="PACF of the first difference 
     of earnings series.")
par(mfrow=c(1,1))

adf.test(diff.earnings)

diff.earnings2 = diff(earnings, differences = 2)
plot(diff.earnings2,type='o',ylab='Earnings ', main='Time series plot of the second difference of earnings series.')

par(mfrow=c(1,2))
acf(diff.earnings2, main="ACF of the second difference 
    of earnings series.")
pacf(diff.earnings2, main="PACF of the second difference 
     of earnings series.")
par(mfrow=c(1,1))
adf.test(diff.earnings2)

pp.test(diff.earnings2)

# d = 2

# Model specification

par(mfrow=c(1,2))
acf(diff.earnings2, main="ACF of the second difference 
    of earnings series.")
pacf(diff.earnings2, main="PACF of the second difference 
     of earnings series.")
par(mfrow=c(1,1))
# {ARIMA(1,2,1)}

eacf(diff.earnings2, ar.max = 5, ma.max = 5)
# {ARIMA(1,2,1), ARIMA(0,2,1)}

res = armasubsets(y=diff.earnings2,nar=5,nma=5,y.name='p',ar.method='ols')
plot(res)

#The final set of possible models is 
# {ARIMA(1,2,1), ARIMA(0,2,1), ARIMA(2,2,1), 
#  ARIMA(3,2,1), ARIMA(4,2,1), ARIMA(5,2,1)}

# ARIMA(0,2,1)
model_021_css = Arima(earnings,order=c(0,2,1),method='CSS')
lmtest::coeftest(model_021_css)

model_021_ml = Arima(earnings,order=c(0,2,1),method='ML')
coeftest(model_021_ml)

# ARIMA(1,2,1)
model_121_css = Arima(earnings,order=c(1,2,1),method='CSS')
coeftest(model_121_css)

model_121_ml = Arima(earnings,order=c(1,2,1),method='ML')
coeftest(model_121_ml)

# ARIMA(2,2,1)
model_221_css = Arima(earnings,order=c(2,2,1),method='CSS')
coeftest(model_221_css)

model_221_ml = Arima(earnings,order=c(2,2,1),method='ML')
coeftest(model_221_ml)

# ARIMA(3,2,1) 
model_321_css = Arima(earnings,order=c(3,2,1),method='CSS')
coeftest(model_321_css)

model_321_ml = Arima(earnings,order=c(3,2,1),method='ML')
coeftest(model_321_ml)

model_321_CSSml = Arima(earnings,order=c(3,2,1),method='CSS-ML')
coeftest(model_321_CSSml)

# ARIMA(4,2,1)
model_421_css = Arima(earnings,order=c(4,2,1),method='CSS')
coeftest(model_421_css)

model_421_ml = Arima(earnings,order=c(4,2,1),method='ML')
coeftest(model_421_ml)

model_421_CSSml = Arima(earnings,order=c(4,2,1),method='CSS-ML')
coeftest(model_421_CSSml)

# ARIMA(5,2,1)
model_521_css = Arima(earnings,order=c(5,2,1),method='CSS')
coeftest(model_521_css)

model_521_ml = Arima(earnings,order=c(5,2,1),method='ML')
coeftest(model_521_ml)


# AIC and BIC values
# you need to source the sort.score() function, which is available in Canvas shell
sort.score(AIC(model_021_ml,model_121_ml,model_221_ml,model_321_ml,model_421_ml,model_521_ml), score = "aic")
sort.score(BIC(model_021_ml,model_121_ml,model_221_ml,model_321_ml,model_421_ml,model_521_ml), score = "bic" )

# Both AIC and BIC select ARIMA(0,2,1) model for this series.

Smodel_021_css <- accuracy(model_021_css)[1:7]
Smodel_121_css <- accuracy(model_121_css)[1:7]
Smodel_221_css <- accuracy(model_221_css)[1:7]
Smodel_321_css <- accuracy(model_321_css)[1:7]
Smodel_421_css <- accuracy(model_421_css)[1:7]
Smodel_521_css <- accuracy(model_521_css)[1:7]
df.Smodels <- data.frame(
  rbind(Smodel_021_css,Smodel_121_css,Smodel_221_css,
        Smodel_321_css,Smodel_421_css,Smodel_521_css)
)
colnames(df.Smodels) <- c("ME", "RMSE", "MAE", "MPE", "MAPE", 
                          "MASE", "ACF1")
rownames(df.Smodels) <- c("ARIMA(0,2,1)", "ARIMA(1,2,1)", "ARIMA(2,2,1)", 
                          "ARIMA(3,2,1)", "ARIMA(4,2,1)", "ARIMA(5,2,1)")
round(df.Smodels,  digits = 3)

# ARIMA(1,2,1) and ARIMA(0,2,2) are overparametrised models for ARIMA(0,2,1).

# ARIMA(0,2,2)
model_022_css = Arima(earnings,order=c(0,2,2),method='CSS')
coeftest(model_022_css)

model_022_ml = Arima(earnings,order=c(0,2,2),method='ML')
coeftest(model_022_ml)

# ARIMA(1,2,1)
model_121_css = Arima(earnings,order=c(1,2,1),method='CSS')
coeftest(model_121_css)

model_121_ml = Arima(earnings,order=c(1,2,1),method='ML')
coeftest(model_121_ml)




#--- Task 2 --- 
# Yearly averages of monthly total values of ATM cash withdrawals ($ million) in Australia between 1994 and 2016. 
data.cash <- read.csv("~/Documents/MATH1318_TimeSeries/tasks/Task6/data.cash.csv", header=FALSE)$V2
class(data.cash)
data.cash = ts(data.cash,start = 1994)
class(data.cash)

plot(data.cash,type='o',ylab="Cash", main='Time series plot of cash series.')

par(mfrow=c(1,2))
acf(data.cash, main="ACF of cash series.")
pacf(data.cash, main="PACF of cash series.")
par(mfrow=c(1,1))

adf.test(data.cash)

qqnorm(data.cash, ylab="data.cash", xlab="Normal Scores")
qqline(data.cash)
shapiro.test(data.cash) 

BC = BoxCox.ar(data.cash)
BC$ci
lambda <- BC$lambda[which(max(BC$loglike) == BC$loglike)]
lambda
BC.cash = (data.cash^lambda-1)/lambda

plot(BC.cash,type='o',ylab = "Cash", main='Time series plot of BC transformed
     cash series.')

qqnorm(BC.cash, ylab="BC data.cash", xlab="Normal Scores")
qqline(BC.cash)
shapiro.test(BC.cash) 

# BC has not helped with normality. 

diff.data.cash = diff(data.cash, differences=1)

plot(diff.data.cash,type='o',ylab="Cash", main='Time series plot of the first difference of cash series.')

adf.test(diff.data.cash)

diff.data.cash2 = diff(data.cash,differences=2)
plot(diff.data.cash2,type='o',ylab="Cash", main='Time series plot of the second difference of cash series.')

adf.test(diff.data.cash2)
pp.test(diff.data.cash2)

# Model specification

par(mfrow=c(1,2))
acf(diff.data.cash2, main="ACF of the second difference 
    of earnings series.")
pacf(diff.data.cash2, main="PACF of the second difference 
     of earnings series.")
par(mfrow=c(1,1))
# {ARIMA(1,2,1)}

eacf(diff.data.cash2,ar.max = 4, ma.max = 4)
# {ARIMA(1,2,1), ARIMA(0,2,1)}

res = armasubsets(y=diff.data.cash2,nar=5,nma=5,y.name='p',ar.method='ols')
plot(res)
# {ARIMA(1,2,1), ARIMA(3,2,1), ARIMA(1,2,4), 
#  ??ARIMA(3,2,4)??}

# Final set of possible models:
# {ARIMA(1,2,1), ARIMA(0,2,1), ARIMA(3,2,1), ARIMA(1,2,4), 
#  ??ARIMA(3,2,4)??}

# ARIMA(1,2,1)
model_121_css = Arima(data.cash,order=c(1,2,1),method='CSS')
coeftest(model_121_css)

model_121_ml = Arima(data.cash,order=c(1,2,1),method='ML')
coeftest(model_121_ml)

# ARIMA(0,2,1)
model_021_css = Arima(data.cash,order=c(0,2,1),method='CSS')
coeftest(model_021_css)

model_021_ml = Arima(data.cash,order=c(0,2,1),method='ML')
coeftest(model_021_ml)

# ARIMA(3,2,1)
model_321_css = Arima(data.cash,order=c(3,2,1),method='CSS')
coeftest(model_321_css)

model_321_ml = Arima(data.cash,order=c(3,2,1),method='ML')
coeftest(model_321_ml)

model_321_cssml = Arima(data.cash,order=c(3,2,1),method='CSS-ML')
coeftest(model_321_cssml)

# ARIMA(1,2,4) 
model_124_css = Arima(data.cash,order=c(1,2,4),method='CSS')
coeftest(model_124_css)

model_124_ml = Arima(data.cash,order=c(1,2,4),method='ML')
coeftest(model_124_ml)

model_124_CSSml = Arima(data.cash,order=c(1,2,4),method='CSS-ML')
coeftest(model_124_CSSml)

# ARIMA(3,2,4) 
model_324_css = Arima(data.cash,order=c(3,2,4),method='CSS')
coeftest(model_324_css)

model_324_ml = Arima(data.cash,order=c(3,2,4),method='ML')
coeftest(model_324_ml)

model_324_CSSml = Arima(data.cash,order=c(3,2,4),method='CSS-ML')
coeftest(model_324_CSSml)

# AIC and BIC values
# you need to source the sort.score() function, which is available in Canvas shell
sort.score(AIC(model_121_ml,model_021_ml,model_321_ml,model_124_ml,model_324_ml), score = "aic")
sort.score(BIC(model_121_ml,model_021_ml,model_321_ml,model_124_ml,model_324_ml), score = "bic" )

# ARIMA(0,2,1) model is the best one among the set of specified models according to AIC/BIC.


Smodel_021_css <- accuracy(model_021_css)[1:7]
Smodel_121_css <- accuracy(model_121_css)[1:7]
Smodel_321_css <- accuracy(model_321_css)[1:7]
Smodel_124_css <- accuracy(model_124_css)[1:7]
Smodel_324_css <- accuracy(model_324_css)[1:7]
df.Smodels <- data.frame(
  rbind(Smodel_021_css,Smodel_121_css,Smodel_321_css,
        Smodel_124_css,Smodel_324_css)
)
colnames(df.Smodels) <- c("ME", "RMSE", "MAE", "MPE", "MAPE", 
                          "MASE", "ACF1")
rownames(df.Smodels) <- c("ARIMA(0,2,1)","ARIMA(1,2,1)", "ARIMA(3,2,1)",
                          "ARIMA(1,2,4)", "ARIMA(3,2,4)")
round(df.Smodels,  digits = 3)

# Looking at the significance tests and and error measures ARIMA(3,2,4) is the best model.

# To check with the overparametrised models we will fit ARIMA(4,2,4) and ARIMA(3,2,5) models.

# ARIMA(4,2,4)
model_424_css = Arima(data.cash,order=c(4,2,4),method='CSS')
coeftest(model_424_css) 

model_424_ml = Arima(data.cash,order=c(4,2,4),method='ML')
coeftest(model_424_ml)

sort.score(AIC(model_121_ml,model_021_ml,model_321_ml,
               model_124_ml,model_324_ml, model_424_ml ), score = "aic")
sort.score(BIC(model_121_ml,model_021_ml,model_321_ml,
               model_124_ml,model_324_ml,model_424_ml ), score = "bic" )

Smodel_021_css <- accuracy(model_021_css)[1:7]
Smodel_121_css <- accuracy(model_121_css)[1:7]
Smodel_321_css <- accuracy(model_321_css)[1:7]
Smodel_124_css <- accuracy(model_124_css)[1:7]
Smodel_324_css <- accuracy(model_324_css)[1:7]
Smodel_424_css <- accuracy(model_424_css)[1:7]
df.Smodels <- data.frame(
  rbind(Smodel_021_css,Smodel_121_css,Smodel_321_css,
        Smodel_124_css,Smodel_324_css,Smodel_424_css)
)
colnames(df.Smodels) <- c("ME", "RMSE", "MAE", "MPE", "MAPE", 
                          "MASE", "ACF1")
rownames(df.Smodels) <- c("ARIMA(0,2,1)","ARIMA(1,2,1)", "ARIMA(3,2,1)",
                          "ARIMA(1,2,4)", "ARIMA(3,2,4)", "ARIMA(4,2,4)")
round(df.Smodels,  digits = 3)

# ARIMA(5,2,4)
model_524_css = Arima(data.cash,order=c(5,2,4),method='CSS')
coeftest(model_524_css) 

model_524_ml = Arima(data.cash,order=c(5,2,4),method='ML')
coeftest(model_524_ml)

# ARIMA(4,2,5)
model_425_css = Arima(data.cash,order=c(4,2,5),method='CSS')
coeftest(model_425_css) 

model_425_ml = Arima(data.cash,order=c(4,2,5),method='ML')
coeftest(model_425_ml)

# ARIMA(3,2,5)
model_325_css = Arima(data.cash,order=c(3,2,5),method='CSS')
coeftest(model_325_css) 

model_325_ml = Arima(data.cash,order=c(3,2,5),method='ML')
coeftest(model_325_ml)


#--- Task 3 ---
arima.data1 <- read.csv("~/Documents/MATH1318_TimeSeries/tasks/Task6/data.sim.csv", header=FALSE)
arima.data1 = ts(arima.data1,start = 1994)
class(arima.data1)
plot(arima.data1,type='o',ylab = "y",main='Time series plots of simulated series')

par(mfrow=c(1,2))
acf(arima.data1, main="ACF of simulated series.")
pacf(arima.data1, main="PACF of simulated series.")
par(mfrow=c(1,1))

adf.test(arima.data1)
pp.test(arima.data1)

qqnorm(arima.data1, ylab="arima.data1", xlab="Normal Scores")
qqline(arima.data1)
shapiro.test(arima.data1) 

summary(arima.data1)

BC = BoxCox.ar(arima.data1 + abs(min(arima.data1))+0.01, lambda = seq(1,2,0.01))
# BC$ci
# lambda <- BC$lambda[which(max(BC$loglike) == BC$loglike)]
# lambda
# BC.arima.data1 = (arima.data1^lambda-1)/lambda

BC.arima.data1 = log(arima.data1 + abs(min(arima.data1))+0.01)

plot(BC.arima.data1,type='o',ylab = "Cash", main='Time series plot of BC transformed
      arima.data1 series.')

qqnorm(BC.arima.data1, ylab="BC arima.data1", xlab="Normal Scores")
qqline(BC.arima.data1)
shapiro.test(BC.arima.data1)

diff.arima.data1= diff(arima.data1,differences=1)

plot(diff.arima.data1,type='o',ylab = "Value", main='Time series plot of the first differnce
      of arima.data1 series.')

par(mfrow=c(1,2))
acf(diff.arima.data1, main="ACF of the first difference 
    of the simulated series.")
pacf(diff.arima.data1, main="PACF of the first difference 
     of the simulated series.")
par(mfrow=c(1,1))

adf.test(diff.arima.data1)
pp.test(diff.arima.data1)

diff.arima.data2= diff(arima.data1, differences=2)

plot(diff.arima.data2,type='o',ylab = "Value", main='Time series plot of the second differnce
      of arima.data1 series.')

par(mfrow=c(1,2))
acf(diff.arima.data2, main="ACF of the second difference 
    of the simulated series.")
pacf(diff.arima.data2, main="PACF of the second difference 
     of the simulated series.")
par(mfrow=c(1,1))

adf.test(diff.arima.data2)
pp.test(diff.arima.data2)

# Model Specification

par(mfrow=c(1,2))
acf(diff.arima.data2, ci.type='ma', main="ACF of the second difference 
    of the simulated series.")
pacf(diff.arima.data2, main="PACF of the second difference 
     of the simulated series.")
par(mfrow=c(1,1))
#{ARIMA(2,2,2), ARIMA(3,2,2)} 

eacf(diff.arima.data2)
#{ARIMA(0,2,2), ARIMA(1,2,2), ARIMA(0,2,3), ARIMA(1,2,3)}

res = armasubsets(y=diff.arima.data1,nar=10,nma=10,y.name='p',ar.method='ols')
plot(res)
#{ARIMA(1,2,0), ARIMA(2,2,0), ARIMA(3,2,0)}

# Final set of possible models:
#{ARIMA(2,2,2), ARIMA(3,2,2), ARIMA(0,2,2), ARIMA(1,2,2), 
# ARIMA(0,2,3), ARIMA(1,2,3), 
# ARIMA(1,2,0), ARIMA(2,2,0), ARIMA(3,2,0)}

# ARIMA(0,2,2)
model_022_css = Arima(arima.data1,order=c(0,2,2),method='CSS')
coeftest(model_022_css)

model_022_ml = Arima(arima.data1,order=c(0,2,2),method='ML')
coeftest(model_022_ml)

# ARIMA(1,2,2)
model_122_css = Arima(arima.data1,order=c(1,2,2),method='CSS')
coeftest(model_122_css)

model_122_ml = Arima(arima.data1,order=c(1,2,2),method='ML')
coeftest(model_122_ml)

# ARIMA(2,2,2)
model_222_css = Arima(arima.data1,order=c(2,2,2),method='CSS')
coeftest(model_222_css)

model_222_ml = Arima(arima.data1,order=c(2,2,2),method='ML')
coeftest(model_222_ml)

# ARIMA(3,2,2)
model_322_css = Arima(arima.data1,order=c(3,2,2),method='CSS')
coeftest(model_322_css)

model_322_ml = Arima(arima.data1,order=c(3,2,2),method='ML')
coeftest(model_322_ml)

# ARIMA(0,2,3)
model_023_css = Arima(arima.data1,order=c(0,2,3),method='CSS')
coeftest(model_023_css)

model_023_ml = Arima(arima.data1,order=c(0,2,3),method='ML')
coeftest(model_023_ml)

# ARIMA(1,2,3)
model_123_css = Arima(arima.data1,order=c(1,2,3),method='CSS')
coeftest(model_123_css)

model_123_ml = Arima(arima.data1,order=c(1,2,3),method='ML')
coeftest(model_123_ml) 

# ARIMA(1,2,0)
model_120_css = Arima(arima.data1,order=c(1,2,0),method='CSS')
coeftest(model_120_css)

model_120_ml = Arima(arima.data1,order=c(1,2,0),method='ML')
coeftest(model_120_ml) 

# ARIMA(2,2,0)
model_220_css = Arima(arima.data1,order=c(2,2,0),method='CSS')
coeftest(model_220_css)

model_220_ml = Arima(arima.data1,order=c(2,2,0),method='ML')
coeftest(model_220_ml) 

# ARIMA(3,2,0)
model_320_css = Arima(arima.data1,order=c(3,2,0),method='CSS')
coeftest(model_320_css)

model_320_ml = Arima(arima.data1,order=c(3,2,0),method='ML')
coeftest(model_320_ml) 

# AIC and BIC values
#{ARIMA(2,2,2), ARIMA(3,2,2), ARIMA(0,2,2), ARIMA(1,2,2), 
# ARIMA(0,2,3), ARIMA(1,2,3), 
# ARIMA(1,2,0), ARIMA(2,2,0), ARIMA(3,2,0)}

sort.score(AIC(model_121_ml,model_022_ml,model_122_ml,model_222_ml,model_322_ml,
               model_023_ml,model_123_ml,model_120_ml,model_220_ml,model_320_ml), score = "aic")
sort.score(BIC(model_121_ml,model_022_ml,model_122_ml,model_222_ml,model_322_ml,
               model_023_ml,model_123_ml,model_120_ml,model_220_ml,model_320_ml), score = "bic" )

# ARIMA(1,2,2) and ARIMA(2,2,2) are the best models by AIC.
# ARIMA(2,2,0) and ARIMA(1,2,2) are the best model by BIC.


Smodel_022_css <- accuracy(model_022_css)[1:7]
Smodel_122_css <- accuracy(model_122_css)[1:7]
Smodel_222_css <- accuracy(model_222_css)[1:7]
Smodel_322_css <- accuracy(model_322_css)[1:7]
Smodel_023_css <- accuracy(model_023_css)[1:7]
Smodel_123_css <- accuracy(model_123_css)[1:7]
Smodel_120_css <- accuracy(model_120_css)[1:7]
Smodel_220_css <- accuracy(model_220_css)[1:7]
Smodel_320_css <- accuracy(model_320_css)[1:7]
df.Smodels <- data.frame(
  rbind(Smodel_022_css,Smodel_122_css,Smodel_222_css,
        Smodel_322_css,Smodel_023_css,Smodel_123_css,
        Smodel_120_css,Smodel_220_css,Smodel_320_css)
)
colnames(df.Smodels) <- c("ME", "RMSE", "MAE", "MPE", "MAPE", 
                          "MASE", "ACF1")
rownames(df.Smodels) <- c("ARIMA(0,2,2)", "ARIMA(1,2,2)", "ARIMA(2,2,2)", 
                          "ARIMA(3,2,2)", "ARIMA(0,2,3)", "ARIMA(1,2,3)",
                          "ARIMA(1,2,0)", "ARIMA(2,2,0)", "ARIMA(3,2,0)")
round(df.Smodels,  digits = 3)


# Overparameterised models for ARIMA(2,2,0) are ARIMA(2,2,1) and ARIMA(3,2,0)

# ARIMA(2,2,1)
model_221_css = Arima(arima.data1,order=c(2,2,1),method='CSS')
coeftest(model_221_css)

model_221_ml = Arima(arima.data1,order=c(2,2,1),method='ML')
coeftest(model_221_ml) 

# ARIMA(3,2,0)
model_320_css = Arima(arima.data1,order=c(3,2,0),method='CSS')
coeftest(model_320_css)

model_320_ml = Arima(arima.data1,order=c(3,2,0),method='ML')
coeftest(model_320_ml) 

