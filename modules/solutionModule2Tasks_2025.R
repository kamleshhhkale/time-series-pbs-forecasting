rm(list=ls()) 
library(TSA)

#--- Task 1 --- 
# The dataset ‘wages’ of the TSA package contains monthly values of the average hourly wages (in dollars) for 
# workers in the U.S. apparel and textile products industry for July 1981 through June 1987.

data(wages)
class(wages)
head(wages)

plot(wages,type='o',ylab='Average hourly wages ')
points(y=wages,x=time(wages),pch=as.vector(season(wages)))


summary(wages)

y = wages 
x = zlag(wages) # lag of 1 year (t and t-1)
index = 2:length(x)
cor(y[index],x[index]) 
plot(y[index],x[index])

x = zlag(zlag(wages)) # lag of 2 years (t and t-2)
index = 3:length(x)
cor(y[index],x[index]) 
plot(y[index],x[index])

acf(wages, lag.max = 40) # lag till 40

# Fit the linear model
t <- time(wages) # taking the time out
dataFreq <- frequency(wages)
model1 = lm(wages ~ t) # label the linear trend model as model1
summary(model1)

plot(wages,type='o',ylab='y', ylim = c(7.5,10))
fitted.model1 <- -549.0061 + 0.2811 *t # high r squared so over fits
lines(fitted.model1)
coef(model1)[1]
coef(model1)[2]

# Residual analysis
res.model1 = rstudent(model1)
# win.graph(width=4.875, height=2.5,pointsize=8) # Use this for Windows computers
# x11() # Use this for Mac computers
par(mfrow=c(2,2))
plot(y = res.model1, x = as.vector(time(wages)),xlab = 'Time', ylab='Standardized Residuals',type='l',main = "Standardised residuals from linear model.")
hist(res.model1,xlab='Standardized Residuals', main = "Histogram of standardised residuals.")
qqnorm(y=res.model1, main = "QQ plot of standardised residuals.")
qqline(y=res.model1, col = 2, lwd = 1, lty = 2)
shapiro.test(res.model1)
acf(res.model1, main = "ACF of standardized residuals.")
# pacf(res.model1, main = "PACF of standardized residuals.")
par(mfrow=c(1,1))

h <- 5 # 5 steps ahead forecasts
lastTimePoint <- t[length(t)]
aheadTimes <- data.frame(t = seq(lastTimePoint+(1/dataFreq), lastTimePoint+5*(1/dataFreq), 1/dataFreq)) 

frc.model1 <- -549.0061 + 0.2811 *aheadTimes$t # This is what the next function does.

frcModel1 <- predict(model1, newdata = aheadTimes, interval = "prediction")

plot(wages, xlim= c(t[1],aheadTimes$t[nrow(aheadTimes)]), ylim = c(7,10), ylab = "Wages series",
     main = "Forecasts from the linear model fitted to the  wages series.")
lines(ts(fitted.model1,start = t[1],frequency = dataFreq)) #  Alternatively you can use abline(model1)
lines(ts(as.vector(frcModel1[,3]), start = aheadTimes$t[1],frequency = dataFreq), col="blue", type="l")
lines(ts(as.vector(frcModel1[,1]), start = aheadTimes$t[1],frequency = dataFreq), col="red", type="l")
lines(ts(as.vector(frcModel1[,2]), start = aheadTimes$t[1],frequency = dataFreq), col="blue", type="l")
legend("topleft", lty=1, pch=1, col=c("black","blue","red"), 
       c("Data","5% forecast limits", "Forecasts"))

# Fit the quadratic model
t = time(wages)
t2 = t^2
model2 = lm(wages~ t + t2) # label the quadratic trend model as model1
summary(model2)

fitted.model2 <- fitted(model2)

plot(ts(fitted.model2), ylim = c(min(c(fitted(model2), as.vector(wages))), max(c(fitted(model2),as.vector(wages)))),
     ylab='y' , main = "Fitted quadratic curve to wages series", type="l",lty=2,col="red")
lines(as.vector(wages),type="o")

# Residual analysis
res.model2 = rstudent(model2)
# win.graph(width=4.875, height=2.5,pointsize=8) # Use this for Windows computers
# x11() # Use this for Mac computers
par(mfrow=c(2,2))
plot(y = res.model2, x = as.vector(time(wages)),xlab = 'Time', ylab='Standardized Residuals',type='l',main = "Standardised residuals from quadratic model.")
hist(res.model2,xlab='Standardized Residuals', main = "Histogram of standardised residuals.")
qqnorm(y=res.model2, main = "QQ plot of standardised residuals.")
qqline(y=res.model2, col = 2, lwd = 1, lty = 2)
shapiro.test(res.model2)
acf(res.model2, main = "ACF of standardized residuals.")
# pacf(res.model2, main = "PACF of standardized residuals.")
par(mfrow=c(1,1))

h <- 5 # 5 steps ahead forecasts
lastTimePoint <- t[length(t)]
aheadTimes <- data.frame(t = seq(lastTimePoint+(1/dataFreq), lastTimePoint+h*(1/dataFreq), 1/dataFreq),
                         t2 =  seq(lastTimePoint+(1/dataFreq), lastTimePoint+h*(1/dataFreq), 1/dataFreq)^2) 

frcModel2 <- predict(model2, newdata = aheadTimes, interval = "prediction")

plot(wages, xlim= c(t[1],aheadTimes$t[nrow(aheadTimes)]), ylim = c(7,10), ylab = "Wages series",
     main = "Forecasts from the linear model fitted to the  wages series.")
lines(ts(fitted.model2,start = t[1],frequency = dataFreq)) #  Alternatively you can use abline(model1)
lines(ts(as.vector(frcModel2[,3]), start = aheadTimes$t[1],frequency = dataFreq), col="blue", type="l")
lines(ts(as.vector(frcModel2[,1]), start = aheadTimes$t[1],frequency = dataFreq), col="red", type="l")
lines(ts(as.vector(frcModel2[,2]), start = aheadTimes$t[1],frequency = dataFreq), col="blue", type="l")
legend("topleft", lty=1, pch=1, col=c("black","blue","red"), 
       c("Data","5% forecast limits", "Forecasts"))


#--- Task 2 --- 
# The dataset “beersales” of the TSA package contains monthly U.S. beer sales (in millions of barrels) 
# for the period January 1975 through December 1990.

data("beersales")
plot(beersales,type='o',ylab='Monthly U.S. beer sales')

 
plot(beersales,type='l',ylab='Sales')
points(y=beersales,x=time(beersales), pch=as.vector(season(beersales)))

y = beersales 
x = zlag(beersales) 
index = 2:length(x)
cor(y[index],x[index]) 
plot(y[index],x[index])

x = zlag(zlag(beersales))
index = 3:length(x)
cor(y[index],x[index]) 
plot(y[index],x[index])

acf(beersales, lag.max = 70)

summary(beersales)

# Fit seasonal model
month.=season(beersales) # period added to improve table display and this line sets up indicators
dataFreq <- frequency(beersales) # period added to improve table display and this line sets up indicators
model3=lm(beersales~month.-1) # -1 removes the intercept term
summary(model3)

plot(ts(fitted(model3)), ylim = c(min(c(fitted(model3), as.vector(beersales))), max(c(fitted(model3),as.vector(beersales)))),
     ylab='y' , main = "Fitted seasonal model to beersales series.", type="l",lty=2,col="red")
lines(as.vector(beersales),type="o")


res.model3 = rstudent(model3)
# win.graph(width=4.875, height=2.5,pointsize=8) # Use this for Windows computers
# x11() # Use this for Mac computers
par(mfrow=c(2,2))
plot(y = res.model3, x = as.vector(time(beersales)),xlab = 'Time', ylab='Standardized Residuals',type='l',main = "Standardised residuals from seasonal model.")
points(y=res.model3,x=time(beersales), pch=as.vector(season(beersales)))
hist(res.model3,xlab='Standardized Residuals', main = "Histogram of standardised residuals.")
qqnorm(y=res.model3, main = "QQ plot of standardised residuals.")
qqline(y=res.model3, col = 2, lwd = 1, lty = 2)
shapiro.test(res.model3)
acf(res.model3, main = "ACF of standardized residuals.")
# pacf(res.model3, main = "PACF of standardized residuals.")
par(mfrow=c(1,1))

# Fit quadratic model
t = time(beersales)
t2 = t^2
model4 = lm(beersales~ t + t2) # label the quadratic trend model as model1
summary(model4)

plot(ts(fitted(model4)), ylim = c(min(c(fitted(model4), as.vector(beersales))), max(c(fitted(model4),as.vector(beersales)))),
     ylab='y' , main = "Fitted quadratic curve to beersales series", type="l",lty=2,col="red")
lines(as.vector(beersales),type="o")

res.model4 = rstudent(model4)
# win.graph(width=4.875, height=2.5,pointsize=8) # Use this for Windows computers
# x11() # Use this for Mac computers
par(mfrow=c(2,2))
plot(y = res.model4, x = as.vector(time(beersales)),xlab = 'Time', ylab='Standardized Residuals',type='l',main = "Standardised residuals from quadratic model.")
points(y=res.model4,x=time(beersales), pch=as.vector(season(beersales)))
hist(res.model4,xlab='Standardized Residuals', main = "Histogram of standardised residuals.")
qqnorm(y=res.model4, main = "QQ plot of standardised residuals.")
qqline(y=res.model4, col = 2, lwd = 1, lty = 2)
shapiro.test(res.model4)
acf(res.model4, main = "ACF of standardized residuals.")
# pacf(res.model4, main = "PACF of standardized residuals.")
par(mfrow=c(1,1))

# Fit seasonal plus quadratic time trend model
model5 = lm(beersales~ month. + t + t2-1 ) # label the quadratic trend model as model1
summary(model5)

fitted.model5 <- fitted(model5)

plot(ts(fitted(model5)), ylim = c(min(c(fitted(model5), as.vector(beersales))), max(c(fitted(model5),as.vector(beersales)))),
     ylab='y' , main = "Fitted seasonal plus quadratic curve to beersales series", type="l",lty=2,col="red")
lines(as.vector(beersales),type="o")

res.model5 = rstudent(model5)
# win.graph(width=4.875, height=2.5,pointsize=8) # Use this for Windows computers
# x11() # Use this for Mac computers
par(mfrow=c(2,2))
plot(y = res.model5, x = as.vector(time(beersales)),xlab = 'Time', ylab='Standardized Residuals',type='l',main = "Standardised residuals from quadratic model.")
points(y=res.model5,x=time(beersales), pch=as.vector(season(beersales)))
hist(res.model5,xlab='Standardized Residuals', main = "Histogram of standardised residuals.")
qqnorm(y=res.model5, main = "QQ plot of standardised residuals.")
qqline(y=res.model5, col = 2, lwd = 1, lty = 2)
shapiro.test(res.model5)
acf(res.model5, main = "ACF of standardized residuals.",lag.max = 60)
# pacf(res.model5, main = "PACF of standardized residuals.")
par(mfrow=c(1,1))

h <- 15 # 15 steps ahead forecasts
t <- time(beersales)
lastTimePoint <- t[length(t)]
aheadTimes <- data.frame(month. = c("January", "February", "March", "April", "May", "June", "July", "August", "September", 
                                    "October", "November", "December","January", "February", "March"),
                         t = seq(lastTimePoint+(1/dataFreq), lastTimePoint+h*(1/dataFreq), 1/dataFreq),
                         t2 =  seq(lastTimePoint+(1/dataFreq), lastTimePoint+h*(1/dataFreq), 1/dataFreq)^2) 

frcModel5 <- predict(model5, newdata = aheadTimes, interval = "prediction")

plot(beersales, xlim= c(t[1],aheadTimes$t[nrow(aheadTimes)]), ylim = c(9,21), ylab = "Beer sales series",
     main = "Forecasts from the quadratic seasonal model fitted to the  beer sales series.")
lines(ts(fitted.model5,start = t[1],frequency = dataFreq), col = "green") #  Alternatively you can use abline(model1)
lines(ts(as.vector(frcModel5[,3]), start = aheadTimes$t[1],frequency = dataFreq), col="blue", type="l")
lines(ts(as.vector(frcModel5[,1]), start = aheadTimes$t[1],frequency = dataFreq), col="red", type="l")
lines(ts(as.vector(frcModel5[,2]), start = aheadTimes$t[1],frequency = dataFreq), col="blue", type="l")
legend("topleft", lty=1, pch=1, col=c("black","blue","red"), 
       c("Data","5% forecast limits", "Forecasts"))



#--- Task 3 ---
# Load the monthly commercial landings dataset of US-NMFS from the file “NMFS_Landings.csv  download” into the workspace of R. 

NMFS_Landings <- read.csv('/Users/kamleshhhkale/Documents/RMIT/SEM 3/TIME SERIES ANALYSIS/NMFS_Landings%282%29.csv')
class(NMFS_Landings)
head(NMFS_Landings)

# Convert data into a time series object
NMFS_Landings.ts = matrix(NMFS_Landings$Metric_Tons, nrow = 25, ncol = 12)
NMFS_Landings.ts = as.vector(t(NMFS_Landings.ts))
NMFS_Landings.ts = ts(NMFS_Landings.ts,start=c(1991,1), end=c(2015,12), frequency=12)
class(NMFS_Landings.ts)

NMFS_Landings.ts
plot(NMFS_Landings.ts,ylab='Landings in metric tons',xlab='Year',type='o', main = "Time series plot of monthly landings in metric tons.")

plot(NMFS_Landings.ts,ylab='Landings in metric tons',xlab='Year',main = "Time series plot of monthly landings in metric tons.")
points(y=NMFS_Landings.ts,x=time(NMFS_Landings.ts), pch=as.vector(season(NMFS_Landings.ts)))

plot(y=NMFS_Landings.ts,x=zlag(NMFS_Landings.ts),ylab='Landings in metric tons', xlab='Previous Year Landings in metric tons' , main = "Scatter plot of neighboring landings in metric tons")

summary(NMFS_Landings.ts) 

y = NMFS_Landings.ts              # Read the Landings in metric tons data into y
x = zlag(NMFS_Landings.ts)        # Generate first lag of the Spawners series
index = 2:length(x)          # Create an index to get rid of the first NA value and the last 5 missing values in x
cor(y[index],x[index]) 
plot(y[index],x[index])

x = zlag(zlag(NMFS_Landings.ts))
index = 3:length(x)
cor(y[index],x[index]) 
plot(y[index],x[index])

acf(NMFS_Landings.ts, lag.max = 70)

# Fit linear model
model.NMFS_Landings.ln = lm(NMFS_Landings.ts~time(NMFS_Landings.ts)) 
summary(model.NMFS_Landings.ln)

plot(NMFS_Landings.ts,type='o',ylab='y')
abline(model.NMFS_Landings.ln)

res.model.NMFS_Landings.ln = rstudent(model.NMFS_Landings.ln)

# win.graph(width=4.875, height=2.5,pointsize=8) # Use this for Windows computers
# x11() # Use this for Mac computers
par(mfrow=c(2,2))
plot(y = res.model.NMFS_Landings.ln, x = as.vector(time(NMFS_Landings.ts)),xlab = 'Time', ylab='Standardized Residuals',type='l',main = "Standardised residuals from linear model.")
points(y=res.model.NMFS_Landings.ln,x=time(NMFS_Landings.ts), pch=as.vector(season(NMFS_Landings.ts)))
hist(res.model.NMFS_Landings.ln,xlab='Standardized Residuals', main = "Histogram of standardised residuals.")
qqnorm(y=res.model.NMFS_Landings.ln, main = "QQ plot of standardised residuals.")
qqline(y=res.model.NMFS_Landings.ln, col = 2, lwd = 1, lty = 2)
shapiro.test(res.model.NMFS_Landings.ln)
acf(res.model.NMFS_Landings.ln, main = "ACF of standardized residuals.")
# pacf(res.model.NMFS_Landings.ln, main = "PACF of standardized residuals.")
par(mfrow=c(1,1))

# Fit quadratic model
t = time(NMFS_Landings.ts)
t2 = t^2
model.NMFS_Landings.qa = lm(NMFS_Landings.ts~ t + t2) 
summary(model.NMFS_Landings.qa)

plot(ts(fitted(model.NMFS_Landings.qa)), ylim = c(min(c(fitted(model.NMFS_Landings.qa),
    as.vector(NMFS_Landings.ts))), max(c(fitted(model.NMFS_Landings.qa),as.vector(NMFS_Landings.ts)))),
    ylab='y' , main = "Fitted quadratic curve to monthly landings in metric tons series", type="l",lty=2,col="red")
lines(as.vector(NMFS_Landings.ts),type="o")

res.model.NMFS_Landings.qa = rstudent(model.NMFS_Landings.qa)

# win.graph(width=4.875, height=2.5,pointsize=8) # Use this for Windows computers
# x11() # Use this for Mac computers
par(mfrow=c(2,2))
plot(y = res.model.NMFS_Landings.qa, x = as.vector(time(NMFS_Landings.ts)),xlab = 'Time', ylab='Standardized Residuals',type='l',main = "Standardised residuals from quadratic model.")
points(y=res.model.NMFS_Landings.qa,x=time(NMFS_Landings.ts), pch=as.vector(season(NMFS_Landings.ts)))
hist(res.model.NMFS_Landings.qa,xlab='Standardized Residuals', main = "Histogram of standardised residuals.")
qqnorm(y=res.model.NMFS_Landings.qa, main = "QQ plot of standardised residuals.")
qqline(y=res.model.NMFS_Landings.qa, col = 2, lwd = 1, lty = 2)
shapiro.test(res.model.NMFS_Landings.qa)
acf(res.model.NMFS_Landings.qa, main = "ACF of standardized residuals.")
# pacf(res.model.NMFS_Landings.qa, main = "PACF of standardized residuals.")
par(mfrow=c(1,1))

# Fit harmonic model
har.=harmonic(NMFS_Landings.ts,1)
data <- data.frame(NMFS_Landings.ts,har.)
model.NMFS_Landings.har = lm(NMFS_Landings.ts ~ cos.2.pi.t. + sin.2.pi.t. , data = data)
summary(model.NMFS_Landings.har)

fitted.model.NMFS_Landings.har <- fitted(model.NMFS_Landings.har)

plot(ts(fitted(model.NMFS_Landings.har)), ylim = c(min(c(fitted(model.NMFS_Landings.har),
     as.vector(NMFS_Landings.ts))), max(c(fitted(model.NMFS_Landings.har),as.vector(NMFS_Landings.ts)))),
     ylab='y' , main = "Fitted quadratic curve to random walk data", type="l",lty=2,col="red")
lines(as.vector(NMFS_Landings.ts),type="o")

res.model.NMFS_Landings.har = rstudent(model.NMFS_Landings.har)

# win.graph(width=4.875, height=2.5,pointsize=8) # Use this for Windows computers
# x11() # Use this for Mac computers
par(mfrow=c(2,2))
plot(y = res.model.NMFS_Landings.har, x = as.vector(time(NMFS_Landings.ts)),xlab = 'Time', ylab='Standardized Residuals',type='l',main = "Standardised residuals from harmonic model.")
points(y=res.model.NMFS_Landings.har,x=time(NMFS_Landings.ts), pch=as.vector(season(NMFS_Landings.ts)))
hist(res.model.NMFS_Landings.har,xlab='Standardized Residuals', main = "Histogram of standardised residuals.")
qqnorm(y=res.model.NMFS_Landings.har, main = "QQ plot of standardised residuals.")
qqline(y=res.model.NMFS_Landings.har, col = 2, lwd = 1, lty = 2)
shapiro.test(res.model.NMFS_Landings.har)
acf(res.model.NMFS_Landings.har, main = "ACF of standardized residuals.")
# pacf(res.model.NMFS_Landings.har, main = "PACF of standardized residuals.")
par(mfrow=c(1,1))

h <- 40 # 20 steps ahead forecasts
t <- time(NMFS_Landings.ts)
lastTimePoint <- t[length(t)]
t1 <- cos(2*pi*seq(lastTimePoint+(1/dataFreq), lastTimePoint+h*(1/dataFreq), 1/dataFreq)) # Create the variables in the model
t2 <- sin(2*pi*seq(lastTimePoint+(1/dataFreq), lastTimePoint+h*(1/dataFreq), 1/dataFreq))

#### Not to do!!! ####
aheadTimes <- data.frame(cos.t. = t1 , sin.t. = t2,
                         t = seq(lastTimePoint+(1/dataFreq), lastTimePoint+h*(1/dataFreq), 1/dataFreq),
                         t2 =  seq(lastTimePoint+(1/dataFreq), lastTimePoint+h*(1/dataFreq), 1/dataFreq)^2) 
summary(model.NMFS_Landings.har)
frc.model.NMFS_Landings.har <- predict(model.NMFS_Landings.har, newdata = aheadTimes, interval = "prediction")

# Check the names of columns you create in aheadTimes and  those in the model summary. 
# aheadTimes has cos.t., sin.t., t, t2
# However, the model has cos.2.pi.t., sin.2.pi.t., t, t2
# "cos.t. and cos.2.pi.t." and "sin.t. and sin.2.pi.t." do not match! 
# So you will get en error whe nyou run predict() function!
#### Not to do!!! ####

aheadTimes <- data.frame(cos.2.pi.t. = t1 , sin.2.pi.t. = t2,
                         t = seq(lastTimePoint+(1/dataFreq), lastTimePoint+h*(1/dataFreq), 1/dataFreq),
                         t2 =  seq(lastTimePoint+(1/dataFreq), lastTimePoint+h*(1/dataFreq), 1/dataFreq)^2) 

frc.model.NMFS_Landings.har <- predict(model.NMFS_Landings.har, newdata = aheadTimes, interval = "prediction")

plot(NMFS_Landings.ts, xlim= c(t[1],aheadTimes$t[nrow(aheadTimes)]), ylim = c(-2*10^5,8*10^5), ylab = "Landings series",
     main = "Forecasts from the quadratic seasonal model fitted to the  beer sales series.")
lines(ts(fitted.model.NMFS_Landings.har,start = t[1],frequency = dataFreq), col = "green") #  Alternatively you can use abline(model1)
lines(ts(as.vector(frc.model.NMFS_Landings.har[,3]), start = aheadTimes$t[1],frequency = dataFreq), col="blue", type="l")
lines(ts(as.vector(frc.model.NMFS_Landings.har[,1]), start = aheadTimes$t[1],frequency = dataFreq), col="red", type="l")
lines(ts(as.vector(frc.model.NMFS_Landings.har[,2]), start = aheadTimes$t[1],frequency = dataFreq), col="blue", type="l")
legend("bottomleft", lty=1, pch=1, col=c("black","blue","red"), 
       c("Data","5% forecast limits", "Forecasts"))


