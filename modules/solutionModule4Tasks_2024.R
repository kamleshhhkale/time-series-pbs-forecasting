library(TSA)
setwd("~/Documents/MATH1318_TimeSeries/tasks/Task4")
#--- Task 1 --- 
# The dataset given in “data.cash.csv” file contains yearly averages of monthly 
# total values of 
# ATM cash withdrawals ($ million) in Australia between 1994 and 2016. 
data.cash <- read.csv("data.cash.csv", header=FALSE)$V2
class(data.cash) 

data.cash = ts(data.cash,start = 1994) 
class(data.cash) 

plot(data.cash,type='o',ylab='ATM cash withdrawals series', main = " Time series plot of cash withdrawals series.")

par(mfrow=c(1,2))
acf(data.cash, main ="ACF plot of cash withdrawals series.")
pacf(data.cash, main ="PACF plot of cash withdrawals series.")
par(mfrow=c(1,1))

qqnorm(data.cash)
qqline(data.cash, col = 2)
shapiro.test(data.cash)

BC = BoxCox.ar(data.cash)
BC$ci
lambda <- BC$lambda[which(max(BC$loglike) == BC$loglike)]
lambda
BC.data.cash = (data.cash^lambda-1)/lambda

plot(BC.data.cash,type='o',ylab='ATM cash withdrawals series', main = " Time series plot of BC transformed cash 
     withdrawals series.")

qqnorm(BC.data.cash)
qqline(BC.data.cash, col = 2)
shapiro.test(BC.data.cash)

# Option: keep going on with the raw series instead of BC transformed series.

diff.BC.data.cash = diff(BC.data.cash)
plot(diff.BC.data.cash,type='o',ylab='First difference of the cash withdrawals series', 
     main ="Time series plot of the first difference of the cash withdrawals series.")

par(mfrow=c(1,2))
acf(diff.BC.data.cash, main ="ACF plot of the first difference of 
    the cash withdrawals series.")
pacf(diff.BC.data.cash, main ="PACF plot of the first difference of 
     the cash withdrawals series.")
par(mfrow=c(1,1))

diff.BC.data.cash2 = diff(BC.data.cash, differences = 2)
plot(diff.BC.data.cash2,type='o', ylab='Second difference of the cash withdrawals series.', 
     main ="Time series plot of the second difference of the cash withdrawals series.")

par(mfrow=c(1,2))
acf(diff.BC.data.cash2, main ="ACF plot of the second difference of 
    the cash withdrawals series.")
pacf(diff.BC.data.cash2, main ="PACF plot of the second difference of 
     the cash withdrawals series.")
par(mfrow=c(1,1))

# Set of possible models: {ARIMA(1,2,1), ARIMA(0,2,1), ARIMA(1,2,0)}


#--- Task 2 --- 
# The dataset “airpass” of TSA package contains international airline passenger monthly totals (in thousands) 
# flown from January 1960 through December 1971. 
data("airpass")
plot(airpass,type='o',ylab='Passenger monthly totals', main = " Time series plot of airline passenger series." )

par(mfrow=c(1,2))
acf(airpass, main ="ACF plot of airline passenger series.", lag.max = 60)
pacf(airpass, main ="PACF plot of airline passenger series.", lag.max = 60)
par(mfrow=c(1,1))

log.airpass = log(airpass)
plot(log.airpass,type='o',ylab='Log-transformed passenger monthly totals', main = " Time series plot of log-transformed 
     airline passenger series." )

par(mfrow=c(1,2))
acf(log.airpass, main ="ACF plot of airline passenger series.", lag.max = 60)
pacf(log.airpass, main ="PACF plot of airline passenger series.", lag.max = 60)
par(mfrow=c(1,1))

BC = BoxCox.ar(airpass,lambda = seq(0.5, 2, 0.02))
BC$ci
lambda <- BC$lambda[which(max(BC$loglike) == BC$loglike)]
lambda
BC.data.airpass = (airpass^lambda-1)/lambda

plot(BC.data.airpass,type='o',ylab='BC-transformed passenger monthly totals', main = " Time series plot of BC-transformed 
     airline passenger series." )

diff.log.airpass = diff(log(airpass))
plot(diff.log.airpass,type='o',ylab='The first difference of log-transformed passenger monthly totals', 
     main = " Time series plot of the first differenced, log-transformed airline passenger series.")

par(mfrow=c(1,2))
acf(diff.log.airpass, main ="ACF plot of the first differenced, 
    log-transformed airline passenger series.", lag.max = 160)
pacf(diff.log.airpass, main ="PACF plot of the first differenced, 
    log-transformed airline passenger series.", lag.max = 160)
par(mfrow=c(1,1))

# SEASONALITY!

#--- Task 3 --- 
# The dataset “larain”of TSA package contains the annual rainfall data for Los Angeles.
data(larain)

plot(larain,type='o', ylab='Los Angelos rainfall series', main = " Time series plot of LA rainfall series." )

par(mfrow=c(1,2))
acf(data.cash, main ="ACF plot of rainfall series.")
pacf(data.cash, main ="PACF plot of rainfall series.")
par(mfrow=c(1,1))

qqnorm(larain)
qqline(larain, col = 2)
shapiro.test(larain)

BC = BoxCox.ar(larain)
BC$ci
lambda <- BC$lambda[which(max(BC$loglike) == BC$loglike)]
lambda
BC.larain = (larain^lambda-1)/lambda


qqnorm(BC.larain)
qqline(BC.larain, col = 2)
shapiro.test(BC.larain)

plot(BC.larain,type='o', ylab='Los Angelos rainfall series', main = " Time series plot of BC transformed LA rainfall series." )

#As another transformation I apply the log transformation
log.larain = log(larain)
plot(log.larain,type='o',ylab='Log-transformed Los Angelos rainfall series', main = "Time series plot of log-transformed LA rainfall series.")

qqnorm(log.larain)
qqline(log.larain, col = 2)
shapiro.test(log.larain)

par(mfrow=c(1,2))
acf(BC.larain, main ="ACF plot of the BC transformed LA rainfall series.")#, lag.max = 60)
pacf(BC.larain, main ="PACF plot of the BC transformed LA rainfall series.")#, lag.max = 60)
par(mfrow=c(1,1))


#--- Task 4 --- 
# The dataset “gold”of TSA package contains the daily price of gold (in dollars per troy ounce) 
# for the 252 trading days of year 2005. 

data(gold, package = 'TSA')

plot(gold,type='o',ylab='Daily gold price', main = "Time series plot of daily price of gold series.")

summary(gold)

par(mfrow=c(1,2))
acf(gold, main ="ACF plot of the daily price of gold series.", lag.max = 160)
pacf(gold, main ="PACF plot of the daily price of gold series.", lag.max = 160)
par(mfrow=c(1,1))

qqnorm(gold)
qqline(gold, col = 2)
shapiro.test(gold)

BC = BoxCox.ar(gold)
BC$ci
lambda <- BC$lambda[which(max(BC$loglike) == BC$loglike)]
lambda
BC.gold = (gold^lambda-1)/lambda

qqnorm(BC.gold)
qqline(BC.gold, col = 2)
shapiro.test(BC.gold)

plot(BC.gold,type='o',ylab='Daily gold price', main = "Time series plot of BC transformed daily price of gold series.")


# Option: keep going on with the BC transformed series.

diff.gold = diff(gold)
plot(diff.gold,type='o',ylab='First difference of daily gold price', main ="Time series plot of the first difference of 
     the daily price of gold series.")

par(mfrow=c(1,2))
acf(diff.gold, main ="ACF plot of the first difference of the daily price of 
    gold series.", lag.max = 60)
pacf(diff.gold, main ="PACF plot of the first difference of the daily price of 
     gold series.", lag.max = 60)
par(mfrow=c(1,1))

# Set of possible models: {ARIMA(1,1,1)}

#--- Task 5 --- 
# arima.data1 = arima.sim(list(order = c(2,1,1), ar = c(0.2,-0.31), ma=-0.1), n = 100, sd=19) 
# plot(arima.data1)
# write.table(arima.data1,file="data.sim.csv",row.names = FALSE,col.names=FALSE)

data.sim <- read.table("data.sim.csv", quote="\"", comment.char="")
data.sim = ts(data.sim)

plot(data.sim,type='o',ylab='Simulated series', main = "Time series plot of the simulated series.")

par(mfrow=c(1,2))
acf(data.sim, main ="ACF plot of the simulated series.")
pacf(data.sim,  main ="PACF plot of the simulated series.")
par(mfrow=c(1,1))

qqnorm(data.sim)
qqline(data.sim, col = 2)
shapiro.test(data.sim)


BC = BoxCox.ar(data.sim) # We get and error message! What to do?

data.sim2 <- data.sim + abs(min(data.sim)) + 0.01

BC = BoxCox.ar(data.sim2) # We get another error message!

BC = BoxCox.ar(data.sim2, lambda = seq(-1, 1, 0.01))
BC$ci
lambda <- BC$lambda[which(max(BC$loglike) == BC$loglike)]
lambda
BC.data.sim = (data.sim2^lambda-1)/lambda

qqnorm(BC.data.sim)
qqline(BC.data.sim, col = 2)
shapiro.test(BC.data.sim)

diff.BC.data.sim = diff(BC.data.sim)
plot(diff.BC.data.sim,type='o',ylab='Simulated data')

par(mfrow=c(1,2))
acf(diff.BC.data.sim, main ="ACF plot of the first differenced, 
    BC transformed simulated series.")
pacf(diff.BC.data.sim, main ="PACF plot of the first differenced,
     BC transformed simulated series.")
par(mfrow=c(1,1))

# Set of possible models: {ARIMA(2,1,2), ARIMA(3,1,2)}
