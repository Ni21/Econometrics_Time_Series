library(urca)
library(lmtest)
library(vars)
library(portes)

data(ChickEgg, 
     package = 'lmtest')

## Examine the data ##

str(ChickEgg)
summary(ChickEgg) #No missing values

chicken <- ChickEgg[,'chicken']
egg <- ChickEgg[,'egg']


boxplot(chicken, 
        xlab = 'chicken',
        col = 'green')

boxplot(egg, # eggs have a few outlier in the top percentile
        xlab = 'egg',
        col = 'yellow')



class(chicken) #Already class timeseries
class(egg) #Also a timeseries object


head(chicken, 
     n=10)
head(egg, 
     n=10)

par(mfrow=c(1,2))
plot.ts(chicken,
        ylab = "Chicken",
        col = 'green',
        lwd = 2)

abline(h=mean(chicken), 
       col ='red')

plot.ts(egg, # Process does not seem to be stationary from the chart. It looks like a structural change in 1940 - 1945.
        xlab = 'Time',
        ylab = "Egg",
        col = 'yellow',
        lwd = 2)

abline(h=mean(egg), 
       col ='red')


#1. Which order of integration do the variables have?

# In the examination, we come to the conclusion that the data might not be stationary.
# Hence, we perform an ADF-test.

### Getting the function in the memory ###
ADF.test <- function(data) {
        ADF <- function(type, data) {
                result1 <- ur.df(data, type = type, lags = 3 * frequency(data), 
                                 selectlags = "AIC")
                DETERM <- ifelse(type == "trend", 2, ifelse(type == "drift", 
                                                            1, 0))
                LAGS <- length(coefficients(result1@testreg)[, 1]) - 
                        DETERM - 1
                result2 <- cbind(t(result1@teststat), result1@cval, coefficients(result1@testreg)["z.lag.1", 
                                                                                                  1], LAGS)
                round(result2, 2)
        }
        types <- c("trend", "drift", "none")
        result3 <- apply(t(types), 2, ADF, data)
        cat(rep("#", 20), "\n")
        cat(rep("#", 20), "\n")
        cat("Augmented Dickey--Fuller test\n")
        cat(rep("#", 20), "\n")
        cat("type:", "  trend ", "drift ", "none\n")
        cat("AR1:   ", result3[[1]][1, 5] + 1, " ", result3[[2]][1, 
                                                                 5] + 1, " ", result3[[3]][5] + 1, "\n")
        cat("lags:  ", result3[[1]][1, 6], "   ", result3[[2]][1, 
                                                               6], "   ", result3[[3]][6], "\n")
        cat(rep("#", 20), "\n")
        result5 <- rbind(result3[[1]][c(1, 3), 1:4], result3[[2]][1:2, 
                                                                  1:4], result3[[3]][1:4])
        rownames(result5)[5] <- "tau1"
        result5
}

dchicken <- diff(chicken)
ADF.test(dchicken)
acf(chicken) # Fading away but high autocorrelation hence I(1)
pacf(chicken) # Spike at 1
# chicken is an arima(1,1,0)

acf(dchicken)
pacf(dchicken)

degg <- diff(egg)
ADF.test(degg)

acf(egg) # Fading away but high autocorrelation hence I(1)
pacf(egg) # Spike at 1
# egg is an arima(1,1,0)

acf(degg)
pacf(degg)


plot.ts(dchicken,
        ylab = "Chicken",
        col = 'green',
        lwd = 2)

abline(h=mean(dchicken), 
       col ='red')

plot.ts(degg, # Process does not seem to be stationary from the chart. It looks like a structural change in 1940 - 1945.
        xlab = 'Time',
        ylab = "Egg",
        col = 'yellow',
        lwd = 2)

abline(h=mean(degg), 
       col ='red')


#2. Show that the variables are not cointegrated?
# Perform an Engle-Granger test for cointegration

result <- lm(chicken ~ egg)
coeftest(result)
ADF.test(residuals(result)) 
# Read statistics against the Engle-Granger critical values

#3. Use information criteria to specify a VAR-model of an appropriate lag order.


var.data <- ts.intersect(dchicken,
                         degg)
colnames(var.data) <- c('Chicken','Egg')
VARselect(var.data, type='both')$selection
var.result <- VAR(var.data,
                  p=1, 
                  type='both')


#4. Perform Ljung-Box and Jarque-Bera tests. What are the conclusions?
# Perform multivariate Ljung Box test
serial.test(var.result) # We cannot reject the null of no autocorrelation -> good!

#Perform multivariate Jarque-Bera Test
normality.test(var.result)$jb.mul$JB # We cannot reject the null of normally distributed residuals -> good!

# No need to adjust for extreme values

roots(var.result) # The eigenvalues are in the unit root and the model is stable (stationary)
# The eigenvalues are positive and real -> exponetial decline of shocks

#5. Perform Granger's causality test. Do chicken Granger-cause eggs or eggs Granger-cause chicken?

causality(var.result,
          cause = 'Chicken')$Granger
# We cannot reject the null hypothesis, hence chicken does not Granger-cause Egg

causality(var.result, 
          cause ='Egg')$Granger
# We reject the null hypothesis: The egg Granger-cause the chicken. Hence, eggs have a predictive power to chicken.

#6. So which was first, the chicken or the egg?

# Granger-causality has nothing to do with causality but is about whether one variable can be used to predict another.
# does not infer anything about causality. Hence, we don't know who was first
