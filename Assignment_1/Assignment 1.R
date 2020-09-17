# Assignment 1

# Niklas Landsberg 940421-T857
# Caroline Kuntz 940915-1462
# Amanda Paasila 970206-2689


##### 1. Load the data in the accompanying file indprod.txt. It contains percentage change in annual production in the US 1861-1988
data <- read.table('indprod.txt',
                                     header = FALSE, 
                                     col.names = '% change in annual production')

str(data)                          # 128 obversations, 1 variable and the class is data frame, hence we need to change the class to ts-object

ts_data <- ts(data, # Converting the data to a time series object
                             start = 1861,
                             frequency = 1)


summary(ts_data) # No missing values. The data set is complete.


###### 2. Estimate the five models AR(p) where p=1,2,3,4,5. 
library(dynlm)
library(lmtest)
library(portes)
library(tseries)

ar <- NULL
for (p in 1:5){                           
  ar[[p]] <- arima(ts_data, order = c(p,0,0)) 
  # Using [[]] so the whole model gets saved in arima-class
}

# 3. Use Akaike's information criterion (AIC) and compare the models. Which model should we choose according to AIC?
#install.packages('portes',lib = .libPaths())

aic <- NULL
for (p in 1:5){   
  aic[p] <- AIC(ar[[p]])
}
aic

###### 4. Perform Ljung-Box and Jarque-Bera tests on the residuals of the selected model. What are the results? Interpret the results

#### Jargue-Bera-Test
jarque.bera.test(residuals(ar[[5]]))

#### Ljung-Box-Test
LjungBox(               
         ar[[5]],       
                lags = 1:10)
# Because the residuals are not normally distributed, we perform ACF/PACF until the residuals are otmally distributed.
# After that we can forecast based on the model.