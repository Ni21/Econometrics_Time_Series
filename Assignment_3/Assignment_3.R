
#1. Construct recursive pseudo-out-of-sample forecasts for horizons h=1,12,24 months
# ahead with data 2000:M01 for evaluation with two different models: (i) AR(1) and
# (ii) ARIMA(p,d,q) using the function auto.arima() in package forecast.

# Loading the data

library(pxweb)
myDataSetName <- 
  get_pxweb_data(url = "http://api.scb.se/OV0104/v1/doris/en/ssd/PR/PR0101/PR0101A/KPI12MNy",
                 dims = list(ContentsCode = c('*'),
                             Tid = c('*')),
                 clean = TRUE)

pi <- ts(myDataSetName$values,
         start = c(1980 ,1),
         end = c(2019, 12),
         frequency = 12)
# (i) AR(2)

library(forecast)

AR.RE.1 <- ts(matrix(NA,
                      nrow = 240,
                      ncol = 1),
               start = c(2000,1),
               frequency = 12)

for(i in 1:240){                                               
  EndDate       <- 1999 + 11/12 + (i - 1) / 12                   
  
  data          <- window(pi,
                          end = EndDate)                    
  model        <- arima(data, 
                        order = c(1,0,0))                         
  AR.RE.1[i, ] <- forecast(model,
                           h = 1)$mean            
}

# For horizon, h = 1, we do not need the for loop to fill empty OUTCOME 
# because each 1 observation can directly be compared to the respective AR.RE.1 forecast
OUTCOME.1 <- window(pi,
                    start = c(2000, 1),
                    end = c(2019, 12))

AR.RE.12 <- ts(matrix(NA,
                     nrow = 229,
                     ncol = 12),
              start = c(2000, 1),
              frequency = 12)

for (i in 1:229){                                               
  EndDate       <- 1999 + 11/12 + (i - 1) / 12                   
  
  data          <- window(pi,
                          end = EndDate)                   
  model        <- arima(data, 
                        order = c(1,0,0))                        
  AR.RE.12[i, ] <- forecast(model,
                           h = 12)$mean             
}

OUTCOME.12 <- ts(matrix(NA,           
                       nrow = 229,  
                       ncol = 12),   
                start = 2000,
                frequency = 12) 

for(i in 1:229){
  StartDate <- 2000 + (i - 1)/12
  EndDate <- StartDate + 11/12
  OUTCOME.12[i, ] <- window(pi, 
                           start = StartDate, 
                           end = EndDate)
}

AR.RE.24 <- ts(matrix(NA,
                      nrow=217,
                      ncol=24),
               start = c(2000,1),
               frequency = 12)

for (i in 1:217){                                               
  EndDate       <- 1999 + 11/12 + (i - 1) / 12                   
  
  data          <- window(pi,
                          end = EndDate)                    
  model        <- arima(data, 
                        order = c(1,0,0))                         
  AR.RE.24[i, ] <- forecast(model,
                            h = 24)$mean            
}

OUTCOME.24 <- ts(matrix(NA,           
                        nrow = 217,  
                        ncol = 24),   
                 start = 2000, 
                 frequency = 12) 

for(i in 1:217){
  StartDate <- 2000 + (i - 1)/12
  EndDate <- StartDate + 23/12
  OUTCOME.24[i, ] <- window(pi, 
                            start = StartDate, 
                            end = EndDate)
}

# (ii) ARIMA(p,d,q)

AUTO.ARIMA.RE.1 <- ts(matrix(NA,
                          nrow = 240,
                          ncol = 1),
                     start = c(2000,1),
                     frequency = 12)

for (i in 1:240){                                               
  EndDate       <- 1999 + 11/12 + (i - 1) / 12                   
  
  data          <- window(pi,
                          end = EndDate)                    
  model        <- auto.arima(data)                         
  
  AUTO.ARIMA.RE.1[i, ] <- forecast(model,
                                   h = 1)$mean            
}


AUTO.ARIMA.RE.12 <- ts(matrix(NA,
                          nrow = 229,
                          ncol = 12),
                   start = c(2000, 1),
                   frequency = 12)

for (i in 1:229){                                               
  EndDate       <- 1999 + 11/12 + (i - 1) / 12                   
  
  data          <- window(pi,
                          end = EndDate)                    
  model        <- auto.arima(data)                         
  
  AUTO.ARIMA.RE.12[i, ] <- forecast(model,
                                    h = 12)$mean            
}

AUTO.ARIMA.RE.24 <- ts(matrix(NA,
                           nrow = 217, 
                           ncol = 24),
                    start = c(2000, 1),
                    frequency = 12)

for (i in 1:217){                                               
  EndDate       <- 1999 + 11/12 + (i - 1) / 12                   
  
  data          <- window(pi,
                          end = EndDate)                    
  model        <- auto.arima(data)                         
  
  AUTO.ARIMA.RE.24[i, ] <- forecast(model,
                                 h = 24)$mean            
}

# 2. Test the forecast for bias

# H0: a=0 (unbiasedness)
bias.test <- function(h,FE){
  require(dynlm); require(lmtest);require(sandwich)
  model <- dynlm(FE[,h] ~ 1)
  matrix <- NeweyWest(model,
                      lag = h-1)
  round(coeftest(model, 
                 vcov. = matrix)[[4]],
        2)
}

# h=1
AR.RE.FE.1 <- AR.RE.1 - OUTCOME.1
AR.RE.FE.p1 <- apply(t(1),
                    2, 
                    bias.test, 
                    AR.RE.FE.1) 
names(AR.RE.FE.p1) <- 1

AUTO.ARIMA.RE.FE.1 <- AUTO.ARIMA.RE.1 - OUTCOME.1

AUTO.ARIMA.RE.FE.p1 <- apply(t(1),
                            2, 
                            bias.test, 
                            AUTO.ARIMA.RE.FE.1) 
names(AUTO.ARIMA.RE.FE.p1) <- 1

rbind(AR.RE.FE.p1 ,AUTO.ARIMA.RE.FE.p1)

# Results: Both do not reject the null of unbiasedness
# Hence, falsely rejecting the null would be high.
# Both are unbiased.

# h=12
AR.RE.FE.12 <- AR.RE.12 - OUTCOME.12
AR.RE.FE.p12 <- apply(t(1:12),
                      2, 
                      bias.test, 
                      AR.RE.FE.12) 
names(AR.RE.FE.p12) <- 1:12

AUTO.ARIMA.RE.FE.12 <- AUTO.ARIMA.RE.12 - OUTCOME.12

AUTO.ARIMA.RE.FE.p12 <- apply(t(1:12),
                              2, 
                              bias.test, 
                              AUTO.ARIMA.RE.FE.12) 
names(AUTO.ARIMA.RE.FE.p12) <- 1:12

rbind(AR.RE.FE.p12 ,AUTO.ARIMA.RE.FE.p12)

# Results: Both do not reject the null of unbiasedness. 
# Hence, falsely rejecting the null would be high.
# Both are unbiased.

# h = 24
AR.RE.FE.24 <- AR.RE.24 - OUTCOME.24
AR.RE.FE.p24 <- apply(t(1:24),
                      2, 
                      bias.test, 
                      AR.RE.FE.24) 
names(AR.RE.FE.p24) <- 1:24

AUTO.ARIMA.RE.FE.24 <- AUTO.ARIMA.RE.24 - OUTCOME.24

AUTO.ARIMA.RE.FE.p24 <- apply(t(1:24),
                              2, 
                              bias.test, 
                              AUTO.ARIMA.RE.FE.24) 
names(AUTO.ARIMA.RE.FE.p24) <- 1:24

rbind(AR.RE.FE.p24, AUTO.ARIMA.RE.FE.p24)

# Results: Both do not reject the null of unbiasedness. 
# Hence, falsely rejecting the null would be high.
# Both are unbiased.

# 3. Compare forecast with the same horizon and test them for equal forecast accuracy

# H0: Equal accuracy
DM.TEST<- function(h,ld){
  require(dynlm); require(lmtest); require(sandwich)
  res <- dynlm(ld[, h] ~1)
  mat <- NeweyWest(res, 
                   lag = h -1)
  round(coeftest(res,
                 vcov. = mat)[4],
        2)
}

# h=1

ld1 <- AR.RE.FE.1^2 - AUTO.ARIMA.RE.FE.1^2


apply(t(1),
      2, 
      DM.TEST,
      ld1)

# Results: The null of equal accuracy is rejected with p-value egauls to 0.01

# h=12

ld12 <- AR.RE.FE.12^2 - AUTO.ARIMA.RE.FE.12^2


apply(t(1:12),
      2, 
      DM.TEST,
      ld12)

# Results: Only for the first/ second forecasts the null of equal precision is rejected and
# after that the forecasts seems to tends to increase their precision

# h=24

ld24 <- AR.RE.FE.24^2 - AUTO.ARIMA.RE.FE.24^2 

apply(t(1:24),
      2, 
      DM.TEST,
      ld24)

# For horizons 1 and 2 we would reject the null of equal forecast precision.
# After that the null is not rejected.