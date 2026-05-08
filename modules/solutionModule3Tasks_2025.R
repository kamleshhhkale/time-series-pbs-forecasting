rm(list=ls())
library(TSA)

# NOTE THAT RESULTS WITH SIMULATED DATA ARE BASED ON A SPECIFIC SEED. 
# IF YOU USE A DIFFERENT SEED, YOU MAY GET SLIGHTLY DIFFERENT RESULTS.

n=100
#--- Task 1 --- 

# You can replace x11() with the code below for Windows computers:
# win.graph(width=10, height=10,pointsize=8)

# a) AR(1)  with  phi = 0.6
set.seed(645135)                    #p d q
# to get the result again and again


sim.dataa = arima.sim(list(order = c(1,0,0), ar = 0.6), n = n) 
# X11()
plot(sim.dataa, main ="Time series plot of simulated AR(1) process.")
# X11()
par(mfrow=c(2,1))
acf(sim.dataa, main ="ACF plot of simulated AR(1) process.")
pacf(sim.dataa, main ="PACF plot of simulated AR(1) process.")
par(mfrow=c(1,1))
# Set of possible models : {ARMA(1,2), ARMA(2,2), AR(1), AR(2) (since slowly decaying pattern)}
# p = 1 PACF q = 2 ACF

# trend downward in ACF so can put q = 0, AR(2) where p has 2 significant correlatio

# b) AR(1)  with  phi = -0.6
set.seed(645135)
sim.datab = arima.sim(list(order = c(1,0,0), ar = -0.6), n = n) 
# X11()
plot(sim.datab, main ="Time series plot of simulated AR(1) process.")
# X11()
par(mfrow=c(2,1))
acf(sim.datab, main ="ACF plot of simulated AR(1) process.")
pacf(sim.datab, main ="PACF plot of simulated AR(1) process.")
par(mfrow=c(1,1))
# Set of possible models : {AR(1)=ARMA(1,0) (since acf is getting up and down, up and down so it is a pattern), ARMA(1,3), ARMA(1,4)} 
# Pattern in ACF so put q = 0 AEMA(1,0) = AR(1)

# c) AR(2) with phi = [0.6,-0.6]
set.seed(645135)
sim.datac = arima.sim(list(order = c(2,0,0), ar = c(0.6, -0.6)), n = n) 
# X11()
plot(sim.datac, main ="Time series plot of simulated AR(2) process.")
# X11()
par(mfrow=c(2,1))
acf(sim.datac, main ="ACF plot of simulated AR(2) process.")
pacf(sim.datac, main ="PACF plot of simulated AR(2) process.")
par(mfrow=c(1,1))
# Set of possible models : {AR(2), ARMA(2,3), ARMA(2,5)} 


# d) MA(1) with theta = 0.5
set.seed(645135)
sim.datad = arima.sim(list(order = c(0,0,1), ma = 0.5), n = n) 
# X11()
plot(sim.datad, main ="Time series plot of simulated MA(1) process.")
# X11()
par(mfrow=c(2,1))
acf(sim.datad, main ="ACF plot of simulated MA(1) process.")
pacf(sim.datad, main ="PACF plot of simulated MA(1) process.")
par(mfrow=c(1,1))
# Set of possible models : {ARMA(0,1),AR(2,0), ARMA(2,1), ARMA(2,2)}
# Significant but not good I mean just 0.3
# Yes but can take 2 for PACF

# e) MA(1) with theta = -0.5
set.seed(645135)
sim.datae = arima.sim(list(order = c(0,0,1), ma = -0.5), n = n) 
# X11()
plot(sim.datae, main ="Time series plot of simulated MA(1) process.")
# X11()
par(mfrow=c(2,1))
acf(sim.datae, main ="ACF plot of simulated MA(1) process.")
pacf(sim.datae, main ="PACF plot of simulated MA(1) process.")
par(mfrow=c(1,1))
# Set of possible models : {MA(1) (slowly decaying pattern for pacf), ARMA(2,1), ARMA(3,1), ARMA(1,1)}
# MA(1)=(0,1)

# f) MA(2) with theta = [-0.5,0.5]
set.seed(645135)
sim.dataf = arima.sim(list(order = c(0,0,2), ma = c(-0.5, 0.5)), n = n) 
# X11()
plot(sim.dataf, main ="Time series plot of simulated MA(2) process.")
# X11()
par(mfrow=c(2,1))
a = acf(sim.dataf, main ="ACF plot of simulated MA(2) process.")
pacf(sim.dataf, main ="PACF plot of simulated MA(2) process.")
par(mfrow=c(1,1))
# Set of possible models : {ARMA(1,2), AR(1) }


# g) ARMA(1,1) with  phi = 0.6; theta = 0.5
set.seed(645135)
sim.datag = arima.sim(list(order = c(1,0,1), ar = 0.6, ma = 0.5), n = n) 
# X11()
plot(sim.datag, main ="Time series plot of simulated ARMA(1,1) process.")
# X11()
par(mfrow=c(2,1))
a = acf(sim.datag, main ="ACF plot of simulated ARMA(1,1) process.")
pacf(sim.datag, main ="PACF plot of simulated ARMA(1,1) process.")
par(mfrow=c(1,1))
# Set of possible models : {AR(1) (slowly decaying pattern in acf), AR(2), ARMA(2,2)}
a

graphics.off()

# h) ARMA(1,1) with  phi = 0.6; theta = -0.5
set.seed(11176)
sim.datah = arima.sim(list(order = c(1,0,1), ar = 0.6, ma = -0.5), n = n)
# X11()
plot(sim.datah, main ="Time series plot of simulated ARMA(1,1) process.")
# X11()
par(mfrow=c(2,1))
a = acf(sim.datah, main ="ACF plot of simulated ARMA(1,1) process.")
pacf(sim.datah, main ="PACF plot of simulated ARMA(1,1) process.")
par(mfrow=c(1,1))
# Set of possible models : {ARMA(1,1), ARMA(2,1), ARMA(2,2)}
# ARMA(2,2): NEAR MISS IN PACCF NEAR MISS IN ACF SO 2,2 INSTEAD 1,1

# i) ARMA(1,1) with  phi = -0.6; theta = 0.5
set.seed(23224343)
sim.datai = arima.sim(list(order = c(1,0,1), ar = -0.6, ma = 0.5), n = n)
# X11()
plot(sim.datai, main ="Time series plot of simulated ARMA(1,1) process.")
# X11()
par(mfrow=c(2,1))
acf(sim.datai, main ="ACF plot of simulated ARMA(1,1) process.")
pacf(sim.datai, main ="PACF plot of simulated ARMA(1,1) process.")
par(mfrow=c(1,1))
# Set of possible models : {AR(1), AR(2), AR(3), ARMA(1,4), ARMA(2,4), ARMA(3,4)}


# j) ARMA(1,1) with  phi = -0.6; theta = -0.5
set.seed(645135)
sim.dataj = arima.sim(list(order = c(1,0,1), ar = -0.6, ma = -0.5), n = n)
# X11()
plot(sim.dataj, main ="Time series plot of simulated ARMA(1,1) process.")
# X11()
par(mfrow=c(2,1))
acf(sim.dataj, main ="ACF plot of simulated ARMA(1,1) process.")
pacf(sim.dataj, main ="PACF plot of simulated ARMA(1,1) process.")
par(mfrow=c(1,1))
# Set of possible models : {AR(1), AR(2), ARMA(1,3), ARMA(2,3), ARMA(1,4), ARMA(2,4)}

# k) ARMA(2,2)  with  phi = [0.6,-0.6]; theta = [-0.5,0.5]
set.seed(910291)
sim.datak = arima.sim(list(order = c(2,0,2), ar = c(0.6,-0.6), ma = c(-0.5,0.5)), n = n)
# X11()
plot(sim.datak, main ="Time series plot of simulated ARMA(2,2) process.")
# X11()
par(mfrow=c(2,1))
acf(sim.datak, main ="ACF plot of simulated ARMA(2,2) process.")
pacf(sim.datak, main ="PACF plot of simulated ARMA(2,2) process.")
par(mfrow=c(1,1))

# Set of possible models : {ARMA(0,0), ARMA(1,1)}

# ARMA(0,0) = WHITE NOISE SERIES

graphics.off()

######################################################
#--- Task 2 --- 

# To reach out to the folder where the data is stored:
setwd("/Users/kamleshhhkale/Documents/RMIT/SEM 3/TIME SERIES ANALYSIS")

arma.data1 = read.table(file="data.sim1.csv") # Notice that I use read.table() here!
arma.data1 = ts(arma.data1)
# X11()
plot(arma.data1, type = "o", main ="Time series plot of simulated dataset - 1.")
# X11()
par(mfrow=c(2,1))
acf(arma.data1, main ="ACF plot of simulated dataset - 1.")
pacf(arma.data1, main ="PACF plot of simulated dataset - 1.")
par(mfrow=c(1,1))
# Set of possible models : {ARMA(3,3). But to extend the set of possible models: {ARMA(3,3), ARMA(2,3), ARMA(3,2)}
# graphics.off()



arma.data2 = read.table(file="data.sim2.csv")
arma.data2 = ts(arma.data2)
# X11()
plot(arma.data2, type = "o", main ="Time series plot of simulated dataset - 2.")
# X11()
par(mfrow=c(2,1))
acf(arma.data2, main ="ACF plot of simulated dataset - 2.")
pacf(arma.data2, main ="PACF plot of simulated dataset - 2.")
par(mfrow=c(1,1))
# Set of possible models : {AR(3)} - non-stationarity!
# graphics.off()



arma.data3 = read.table(file="data.sim3.csv")
arma.data3 = ts(arma.data3)
# X11()
plot(arma.data3, type = "o", main ="Time series plot of simulated dataset - 3.")
# X11()
par(mfrow=c(2,1))
a = acf(arma.data3, main ="ACF plot of simulated dataset - 3.")
pacf(arma.data3, main ="PACF plot of simulated dataset - 3.")
par(mfrow=c(1,1))
# Set of possible models : {???} - non-stationarity!
# graphics.off()

# This is slowly decaying pattern



arma.data4 = read.table(file="data.sim4.csv")
arma.data4 = ts(arma.data4)
# X11()
plot(arma.data4, type = "o", main ="Time series plot of simulated dataset - 4.")
# X11()
par(mfrow=c(2,1))
acf(arma.data4, main ="ACF plot of simulated dataset - 4.")
pacf(arma.data4, main ="PACF plot of simulated dataset - 4.")
par(mfrow=c(1,1))
# Set of possible models : {ARMA(1,0), ARMA(2,0), ARMA(1,2), ARMA(1,3), ARMA(2,2), ARMA(2,3)} - non-stationarity??
# graphics.off()



arma.data5 = read.table(file="data.sim5.csv")
arma.data5 = ts(arma.data5)
# X11()
plot(arma.data5, type = "o", main ="Time series plot of simulated dataset - 5.")
# X11()
par(mfrow=c(2,1))
acf(arma.data5, main ="ACF plot of simulated dataset - 5.")
pacf(arma.data5, main ="PACF plot of simulated dataset - 5.")
par(mfrow=c(1,1))
# Set of possible models : {ARMA(2,2), ARMA(3,2), ARMA(2,3) }
# graphics.off()




# True models: arma.data1: ARMA(2,3); arma.data2: Non-stationary; arma.data3: Non-stationary;
# arma.data4: ARMA(2,0); arma.data5: ARMA(0,2)
