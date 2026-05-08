rm(list=ls())
library(TSA)
library(tseries)

#--- Task 1 ---
# The dataset “gold” contains the daily price of gold (in dollars per troy ounce) for 
# the 252 trading days of year 2005. 
data(gold) # Load data. Here no need to convert it to ts object because.
# it is already ts object. See below:
class(gold)

par(mfrow=c(1,1))
plot(gold,type='o',ylab="Gold price", main='Time series plot of gold series')
# In the time series plot there is trend, no seasonality, no obvious intervention, and no changing 
# variance and succeeding observations imply the existence of autoregressive 
# behavior.
summary(gold)

par(mfrow=c(1,2))
acf(gold, main = "ACF plot of gold series.")
pacf(gold, main = "PACF plot of gold series.")
par(mfrow=c(1,1))

adf.test(gold)

qqnorm(y=gold, main = "QQ plot of gold series.")
qqline(y=gold, col = 2, lwd = 1, lty = 2)
shapiro.test(gold)

BC <- BoxCox.ar(gold) #,lambda = seq(-1, 0.5, 0.01) If you get an error.
BC$ci
lambda <- BC$lambda[which(max(BC$loglike) == BC$loglike)]
lambda
BC.gold = ((gold^lambda)-1)/lambda

par(mfrow=c(1,1))
plot(BC.gold,type='o', ylab ="Gold price", main="Time series plot of BC transformed gold series.")

qqnorm(y=BC.gold, main = "QQ plot of BC transformed purchase values.")
qqline(y=BC.gold, col = 2, lwd = 1, lty = 2)
shapiro.test(BC.gold)

adf.test(BC.gold)

gold.diff <- diff(gold)
par(mfrow=c(1,1))
plot(gold.diff,type='o', ylab ="Gold price", main="Time series plot of the first difference of gold series.")


# Let's apply ADF unitroot test to the differenced series.
adf.test(gold.diff)
kpss.test(gold.diff)
pp.test(gold.diff)
# With a p-value of 0.01, we reject the null hypothesis stating that
# the series is non-stationary; hence, we conclude that the first differencing 
# make the series stationary.

# Specify the orders using ACF, PACF, EACF and BIC over the differenced series
par(mfrow=c(1,2))
acf(gold.diff, main = "ACF plot of the first 
    difference of gold series.")
pacf(gold.diff, main = "PACF plot of the first 
     difference of gold series.")
par(mfrow=c(1,1))
#{ARIMA(1,1,1)}

eacf(gold.diff)
#{ARIMA(0,1,1), ARIMA(1,1,1)}

res = armasubsets(y=gold.diff,nar=5,nma=5,y.name='p',ar.method='ols')
plot(res)
#{ARIMA(0,1,1), ARIMA(1,1,1), ARIMA(3,1,0), ARIMA(1,1,0)}

# Final set of possible models:
# {}

#--- Task 2 ---
# The dataset “JJ” contains quarterly earnings per share for the Johnson & Johnson Company.
data(JJ) # Load data. Here no need to convert it to ts object because
# it is already ts object. See below:
class(JJ)

par(mfrow=c(1,1))
plot(JJ,type='o',ylab='Earnings ', main='Time series plot of quarterly earnings.')
# In the time series plot there is trend, a repeating pattern (seasonality) and 
# changing (increasing variance), no intervention point, and bouncing observations around 
# the mean level imply the existence of moving average behavior.

par(mfrow=c(1,2))
acf(JJ, main='ACF plot of quarterly earnings.')
pacf(JJ,  main='PACF plot of quarterly earnings.')
par(mfrow=c(1,1))
# Slowly decaying pattern in ACF and very high first correlation in PACF
# implies the existence of trend and nonstationarity.

adf.test(JJ)


qqnorm(y=JJ, main = "QQ plot of earnings series.")
qqline(y=JJ, col = 2, lwd = 1, lty = 2)
shapiro.test(JJ)

# Let's first apply the box-Cox transformation.
BC = BoxCox.ar(JJ)
BC$ci
lambda <- BC$lambda[which(max(BC$loglike) == BC$loglike)]
lambda
BC.jj = (JJ^lambda-1)/lambda # apply the Box-Cox transformation

par(mfrow=c(1,1))
qqnorm(BC.jj)
qqline(BC.jj, col = 2)
shapiro.test(BC.jj)

# In the QQ plot the tails of the distribution is far from the normality.
# The p-value of Shapiro test is less than 0.05; hence, we have enough 
# evidence to reject the normality hypothesis. In conclusion, the Box-Cox
# transformation did not significantly help to improve normality of the observations.

plot(BC.jj,type='o',ylab='Earnings ', main='Box-Cox transformed Quarterly earnings.')
# In the time series plot, we observe that the variation in the series is 
# decreased after applying the Box-Cox transformation. But there is still trend.
# So let's apply the first difference and see if it helps.

adf.test(BC.jj)

diff.BC.jj = diff(BC.jj)

par(mfrow=c(1,1))
plot(diff.BC.jj,type='o',ylab='Quarterly earnings ', main='Time series plot of the first differenced, 
     Box-Cox transformed Quarterly earnings.')
# Now, there is only changing variance in the series after taking the first difference.
# Let's go on with the specification of the models.
adf.test(diff.BC.jj)
# The null hypothesis of non-stationarity is rejected with a p-value of 0.01; hence,
# we conclude that the first differencing make the series stationary.

par(mfrow=c(1,2))
acf(diff.BC.jj, main='ACF plot of first 
    differenced quarterly earnings.')
pacf(diff.BC.jj, main='PACF plot of first 
    differenced quarterly earnings.')
par(mfrow=c(1,1))
# {ARIMA(4,1,0), ARIMA(3,1,0)}

eacf(diff.BC.jj)
# {ARIMA(5,1,1), ARIMA(5,1,2), ARIMA(6,1,1)}
par(mfrow=c(1,1))
res2 = armasubsets(y=diff.BC.jj,nar=14,nma=14,y.name='p',ar.method='ols')
plot(res2)
# {ARIMA(4,1,3), ARIMA(4,1,11)}

#{ARIMA(4,1,0), ARIMA(3,1,0),ARIMA(4,1,3), ARIMA(4,1,11),ARIMA(5,1,1), ARIMA(5,1,2), ARIMA(6,1,1)}

#--- Task 3 ---
# The dataset given in “unemployment.csv” contains yearly averages of monthly total number of 
# unemployed people (in thousands) in Australia between 1978 and 2016.
# Read data into R
unemployment <- read.table("/Users/haydardemirhan/Documents/MATH1318_TimeSeries/tasks/Task5/unemployment.csv", 
                           quote="\"", comment.char="")
# Convert to ts object
unemployment = ts(unemployment, start = 1978)
class(unemployment)

par(mfrow=c(1,1))
plot(unemployment,type='o',ylab='Time series plot of yearly average unemployment numbers.')
# In the time series plot there is trend, no seasonality and no changing 
# variance, no intervention point, and succeeding observations imply the  
# existence of autoregressive behavior.

par(mfrow=c(1,2))
acf(unemployment, main='ACF plot of unemployment series.')
pacf(unemployment, main='PACF plot of unemployment series.')
par(mfrow=c(1,1))
# Slowly decaying pattern in ACF and very high first correlation in PACF
# implies the existence of trend and nonstationarity.

adf.test(unemployment)

par(mfrow=c(1,1))
qqnorm(unemployment)
qqline(unemployment, col = 2)
shapiro.test(unemployment)

# Let's first apply the box-Cox transformation.
par(mfrow=c(1,1))
BC = BoxCox.ar(unemployment)
BC$ci
lambda <- BC$lambda[which(max(BC$loglike) == BC$loglike)]
lambda
BC.unemployment = (unemployment^lambda-1)/lambda

plot(BC.unemployment,type='o',ylab='Time series plot of BC transformed 
     yearly average unemployment numbers.')

par(mfrow=c(1,1))
qqnorm(BC.unemployment)
qqline(BC.unemployment, col = 2)
shapiro.test(BC.unemployment)

diff.unemployment = diff(unemployment)

par(mfrow=c(1,1))
plot(diff.unemployment,type='o',ylab='Average unemployment numbers', main = "Time series plot of the first differenced
     yearly average unemployment numbers.")
# There is still some trend in there series
adf.test(diff.unemployment)

par(mfrow=c(1,2))
acf(diff.unemployment, main='ACF plot of first differenced
    unemployment series.')
pacf(diff.unemployment, main='PACF plot of first differenced
    unemployment series.')
par(mfrow=c(1,1))

# There is no slowly decaying pattern in ACF. But there is some trend in the time series plot of the first differenced series.
# So go for another unit-root test.
pp.test(diff.unemployment)
# PP test also leans towards stationarity. 

par(mfrow=c(1,2))
acf(diff.unemployment, main='ACF plot of first differenced
    unemployment series.')
pacf(diff.unemployment, main='PACF plot of first differenced
    unemployment series.')
par(mfrow=c(1,1))
# {ARIMA(1,1,1)}

eacf(diff.unemployment) # Returns an error message!
eacf(diff.unemployment, ar.max = 10, ma.max = 10) # R throws an error with ar.max = 10, ma.max = 10
eacf(diff.unemployment, ar.max = 5, ma.max = 5) # So reduce the numbers to avoid
# {ARIMA(0,1,1), ARIMA(0,1,2), ARIMA(1,1,1), ARIMA(1,1,2)}

par(mfrow=c(1,1))
res3 = armasubsets(y=diff.unemployment,nar=5,nma=5
                   ,y.name='p',ar.method='ols')
plot(res3)
# {ARIMA(4,1,0), ARIMA(5,1,0)}

# Final set of possible models:
# {ARIMA(0,1,1), ARIMA(0,1,2), ARIMA(1,1,1), ARIMA(1,1,2),ARIMA(4,1,0), ARIMA(5,1,0)}

#--- Task 4 ---
# The dataset given in “earnings.csv” file contains average weekly earnings (AUD) 
# in Australia on each November between 1994 and 2016.
# Read data into R
earnings <- read.table("/Users/haydardemirhan/Documents/MATH1318_TimeSeries/tasks/Task5/earnings.csv", 
                       quote="\"", comment.char="")
# Convert to ts object
earnings = ts(earnings,start = c(1994,11), frequency = 12)
class(earnings)

par(mfrow=c(1,1))
plot(earnings,type='o',ylab='Average weekly earnings',main = "Time series plot of average weekly earnings.")
# In the time series plot there is trend, no seasonality, no intervention, and no changing 
# variance and succeeding observations imply the existence of autoregressive 
# behavior.

par(mfrow=c(1,2))
acf(earnings, main='ACF plot of earnings series.')
pacf(earnings, main='PACF plot of earnings series.')
par(mfrow=c(1,1))
# Slowly decaying pattern in ACF and very high first correlation in PACF
# implies the existence of trend and nonstationarity.

adf.test(earnings)

par(mfrow=c(1,1))
qqnorm(earnings)
qqline(earnings, col = 2)
shapiro.test(earnings)

# Let's first apply the box-Cox transformation.
BC = BoxCox.ar(earnings)
BC$ci
lambda <- BC$lambda[which(max(BC$loglike) == BC$loglike)]
lambda
BC.earnings = (earnings^lambda-1)/lambda

par(mfrow=c(1,1))
plot(BC.earnings,type='o',ylab='Average weekly earnings',main = "Time series plot of BC transformed average weekly earnings.")
# Because there was no changing variance, transformation had no effect on the series.
# we will go on with the original series.

par(mfrow=c(1,1))
qqnorm(BC.earnings)
qqline(BC.earnings, col = 2)
shapiro.test(BC.earnings)

diff.earnings = diff(earnings)

par(mfrow=c(1,1))
plot(diff.earnings,type='o',ylab='Average weekly earnings', main = "Time series plot of the first differenced average weekly earnings.")
# There is still some trend in there series
adf.test(diff.earnings)
# ADF test concludes with the p-value of 0.7564 that the series is still non-stationary at 5% level of significance.
# So we will apply the second differencing.

diff.earnings2 = diff(earnings, differences = 2)

par(mfrow=c(1,1))
plot(diff.earnings2,type='o',ylab='Average weekly earnings', main = "Time series plot of the second differenced average weekly earnings.")
# There is still some trend in the series
adf.test(diff.earnings2)
# ADF test concludes with the p-value of 0.07604 that the series is still non-stationary at 5% level of significance.
# But this p-value is close to the threshold 0.05. So we will look into ACF and PACF and PP test.

par(mfrow=c(1,2))
acf(diff.earnings2, main='ACF plot of the second 
    differenced earnings series.')
pacf(diff.earnings2, main='PACF plot of the second 
     differenced earnings series.')
par(mfrow=c(1,1))

pp.test(diff.earnings2)
# From ACF and PACF and PP test, we conclude that the second difference of the series is stationary.


par(mfrow=c(1,2))
acf(diff.earnings2, main='ACF plot of the second 
    differenced earnings series.')
pacf(diff.earnings2, main='PACF plot of the second .
     differenced earnings series.')
par(mfrow=c(1,1))
# {ARIMA(0,2,1), ARIMA(0,2,2),ARIMA(1,2,1), ARIMA(1,2,2)}

eacf(diff.earnings2)
eacf(diff.earnings2, ar.max = 5, ma.max = 5) 
# We put these argument to limit the orders p and q at 5. Otherwise, the eacf()
# function returns an error and displays nothing. 

# {ARIMA(0,2,1), ARIMA(1,2,1)}

par(mfrow=c(1,1))
res4 = armasubsets(y=diff.earnings2,nar=5,nma=5,y.name='p',ar.method='ols')
plot(res4)

# {ARIMA(1,2,1), ARIMA(2,2,1), ARIMA(3,2,1), ARIMA(4,2,1), ARIMA(5,2,1),
# ARIMA(1,2,0), ARIMA(2,2,0), ARIMA(3,2,0), ARIMA(4,2,0), ARIMA(5,2,0)}

# Final set of possible models:
# {ARIMA(0,2,1), ARIMA(0,2,2),ARIMA(1,2,1), ARIMA(1,2,2),
# ARIMA(2,2,1), ARIMA(3,2,1), ARIMA(4,2,1), ARIMA(5,2,1),
# ARIMA(1,2,0), ARIMA(2,2,0), ARIMA(3,2,0), ARIMA(4,2,0), ARIMA(5,2,0)}

