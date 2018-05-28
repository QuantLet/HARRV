[<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/banner.png" width="888" alt="Visit QuantNet">](http://quantlet.de/)

## [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/qloqo.png" alt="Visit QuantNet">](http://quantlet.de/) **futuresPriceSimulation** [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/QN2.png" width="60" alt="Visit QuantNet 2.0">](http://quantlet.de/)

```yaml

Name of QuantLet: futuresPriceSimulation


Published in: 'Realized volatility of CO2 futures'


Description: 'Simulates futures prices based on an estimated HAR model with an MC simulation and plots volatility structure'


Keywords: HAR, monte-carlo, volatility


See also: HARmodels

Author: Thijs Benschop

Submitted: Thijs Benschop

Datafile: simdata.txt

Example : 'The realized volatility structure based on the simulated option prices'
```

![Picture1](futuresPriceSimulation_1.png)

![Picture2](volSmile.png)

### R Code
```r

## Clear history
rm(list = ls(all = TRUE))
graphics.off()

# Set working directory
setwd("")

##### Simulation of option prices on futures #####
# Read data
realVol = read.table("simdata.txt")

# Step 1: Simulate volatility
startDate = as.numeric(as.Date("2011-04-01"))              # start date (year, month, day)
coef      = c(0.5333709, 0.1831072, 0.3975752, 0.1714059)  # coefficients estimated HAR model (intercept, RVd, RVw, RVm)
startVals = realVol[(which(realVol[,"Date"] == startDate) - 21):which(realVol[,"Date"] == startDate), "RVd"] 
                                                           # RVd over past month
S0        = realVol[which(realVol[,"Date"] == startDate), "Last"] 
                                                           # Price of future at time t = 0
meanInov  = -0.0000000000000002899414                      # mean of residuals of fitted model
sdInov    = 2.16635                                        # sd of residuals of fitted model

# Set seed for random number generator
set.seed(123524)

# Initialize object for simulated RVd
resSim = c()
# Simulate RVd for 252 days
for(i in 1:252){
  newVal    = coef[1] + coef[2] * startVals[21] + coef[3] * mean(tail(startVals, n = 5)) + coef[4] * mean(tail(startVals, n = 21)) + rnorm(1, meanInov, sdInov)
  resSim    = c(resSim, newVal)
  startVals = c(startVals[2:21], newVal)
}

# Plot simulated RVd
plot(resSim, type = "l")

# Step 2: Simulate futures prices with stochastic volatility model and compute option prices
# Function to simulate option prices
MCOptionPrice = function(r, gamma, resSim, k, numSim) { 
  # r,       risk free rate
  # gamma,   market price of risk
  # resSim,  simulated RVd
  # k,       strike price
  # numSim,  number of simulated MC paths
  
  v.sum = 0 
  for (i in 1:numSim) {
    st    = MCPath(r, gamma, resSim) 
    s.avg = sum(st) / length(st)           # arithmetic average 
    v     = exp(-r) * pmax(0, s.avg - k)   # discounted payoff (average price) 
    v.sum = v.sum + v                      # sum for MC estimate 
  } 
  v.mc = v.sum / numSim  # MC estimate 
  return(v.mc) 
} 

# Function for simulating futures prices
MCPath = function(r, gamma, IV) { 
  # r,     risk free rate
  # gamma, market price of risk
  # IV,    integrated volatility estimated by realized volatility
  st = S0 + sapply(resSim, function(x){rnorm(n = 1, mean = r + (gamma - 0.5) * abs(x), sd = abs(x))})
  return(st) 
} 

# Test MCOptionPrice
MCOptionPrice(r = 0.01, gamma =  2, resSim, k = 1, numSim =  1000)

# Step 3: Compute implied volatility for different K (strike prices)
# Function to compute implied volatility
# Price, (simulated) option price
# S0, futures price at time t=0
# K, trike price
# Rf, risk-free rate
# ttm, time to maturity
implVol = function(Price, S0, K, Rf, ttm){
  # Price at lower bound
  low    = 0.0001
  d1     = (log(S0 / K) + (low / 2) * ttm) / (sqrt(low) * sqrt(ttm))
  d2     = d1 - (sqrt(low) * sqrt(ttm))
  lowVal = exp(-Rf * ttm) * (S0 * pnorm(d1) - K * pnorm(d2))
  
  # Price at upperbound
  up    = 1000
  d1    = (log(S0 / K) + (up / 2) * ttm) / (sqrt(up) * sqrt(ttm))
  d2    = d1 - (sqrt(up) * sqrt(ttm))
  upVal = exp(-Rf * ttm) * (S0 * pnorm(d1) - K * pnorm(d2))
  
  # While loop until difference between up and low is smaller than 0.0001
  while(up - low > 0.0001){
    new  = (up + low) / 2
    d1   = (log(S0 / K) + (new / 2) * ttm) / (new * sqrt(ttm))
    d2   = d1 - (new * sqrt(ttm))
    tempVal = exp(-Rf * ttm) * (S0 * pnorm(d1) - K * pnorm(d2)) # futures price Black formula
    if(tempVal < Price){
      low    = new
      lowVal = tempVal
    } else {
      up     = new
      upVal  = tempVal
    }
  }
  return(up) 
}

# Simulated option prices for a set of strike prices
a = sapply(seq(from = 15, to = 20, by = 0.2), function(x){MCOptionPrice(r = 0.01, gamma = 1, resSim, k = x, numSim = 1000)})
# Add strike prices
a = cbind(a, seq(from = 15, to = 20, by = 0.2))

# Implied volatility for K = 15-20 with increments of 0.2
c = c()
for(i in 1:26){
  c = c(c, implVol(Price = a[i, 1], K = a[i, 2], Rf = 0.01, ttm = 1, S0 = S0))
}

# Step 4: Plot implied volatility smile
# Add moneyness
d = cbind (seq(from = 15, to = 20, by = 0.2) / S0, c)

# Plot moneyness and implied volatility (volatilty smile)
plot(d[,2] ~ d[,1], ylim = c(0,0.4), type = "l", xlab = "Moneyness", ylab = "Implied volatility", lwd = 2)

# Plot prices
png(file = "volSmile.png", width = 12000, height = 6000, res = 900)
plot(d[,2] ~ d[,1], ylim = c(0,0.4), type = "l", xlab = "Moneyness", ylab = "Implied volatility", lwd = 2)
dev.off()

```

automatically created on 2018-05-28