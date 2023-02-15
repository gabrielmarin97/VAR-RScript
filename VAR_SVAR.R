#------------------------------------------------------------------------------
# VAR, SVAR Estimation - Macroeconomic Policy, UA Spring 2023
# Gabriel Marin, IDB
# gmarinmunoz@iadb.org
# Note: This R-Script seeks to provide a guideline on the estimation of 
# VARs, SVARs, and some introduction to time-series econometrics. 
#------------------------------------------------------------------------------
# Clear all
rm(list = ls())
# Required libraries.
library(haven)
library(ggplot2) # cool graphs
library(tictoc)
library(tidyverse)
library(xts) # to work with time series data
library(lubridate)
library(stargazer) # Regression output in a pretty format
library(vars) # to estimate VAR, SVAR models
library(ggthemes) # cool graphs
library(tseries) # to work with time series data
# Set the working directory
setwd("/Users/gmm/Dropbox/BID/Macropolicy Course")

# We are using data from INEGI and SHCP to estimate fiscal multipliers
# based on Ilzetski, Mendoza, and Vegh (2011)
data <- read.csv(file = "Inputs.csv", header = TRUE, sep = ",")

# The dataset contains several variables from the regression posed in IMV(2011),
# the real GDP, the real exchange rate, the current account, and the government and tax revenues
# of several countries. 

# In our case, we will use as a case study our country: Mexico.

## 1.- Auto Regressive Processes

# Before we make use of the data, let us create visually an AR(1) process of 100 observations
set.seed(123)
n <- 100
simulated_ar <- simulated_ar1 <- simulated_ar2 <- rep(NA,n)

# Create a for loop to fill with an error term

rho <- 0.95 # High Persistence
rho1 <- 0.50 # Low Persistence
rho2 <- 1 # Random walk

simulated_ar[1] = simulated_ar1[1] = simulated_ar2[1] = 0

# Simulate the three stochastic processes with different Rho term.
for (i in 1:n){
  simulated_ar[i+1] = rho*simulated_ar[i] + rnorm(1)
  simulated_ar1[i+1] = rho1*simulated_ar1[i] + rnorm(1)
  simulated_ar2[i+1] = rho2*simulated_ar2[i] + rnorm(1)
}

# Transform to data frame for ggplot
simulated_ar <- as.data.frame(simulated_ar)
simulated_ar1 <- as.data.frame(simulated_ar1)
simulated_ar2 <- as.data.frame(simulated_ar2)


index = seq(1:(n+1))
simulated_ar95 <- ggplot(data = simulated_ar, aes(x =index, y=simulated_ar)) +
                  geom_path(color="darkblue") +
                  geom_path(data = simulated_ar1, aes(x =index, y=simulated_ar1), color="darkred") + 
                  geom_path(data = simulated_ar2, aes(x =index, y=simulated_ar2)) + 
                  theme_stata() +
                  ggtitle("Simulated AR(1) Process") +
                  xlab("Period") + ylab("") 
                  
                  
# Show the final plot
simulated_ar95

## 2.- Stationarity 
# We will evaluate whether our data is stationary, for that we use the Augmented Dickey Fuller
# Test, which is an extension of the Dickey-Fuller Test.

# IMV(2011) makes use of 4 main variables, we will evaluate each of them
names(data)
var_data <- ts(data = data[c("rgdp_mex", "gexp_mex", "ca_mex",
                      "d_reer_mex")],
                start = c(2000,1), frequency = 4)

# It appears visually some of our variables are not stationary!
# Which variables appear to not be stationary?
plot(var_data)

# ADF Test
adf.test(var_data[,"rgdp_mex"]) # Real GDP, not stationary 
adf.test(var_data[,"gexp_mex"]) # Gov Expenditure, not stationary
adf.test(var_data[,"ca_mex"]) # Current Account, not stationary 
adf.test(var_data[,"d_reer_mex"]) # First-difference of REER is stationary

# This suggest we might need to include a time dummy in our regression specification
# We do not demean the series and estimate the VAR because of interpretability
# It is hard to explain shocks to the change of a variable in the IRFs, better
# to use log transformations.


## 3.- VAR Models
# Select optimal lag by using information criterion.
VARselect(y = var_data, lag.max = 12, type = c("trend"))
# Lag selection chooses 1 lag as optimal.

# Estimate a VAR with 1 lag and include a trend
# We select the variables that way because ordering matters in SVAR models,
# R uses the VAR model as input
var_model <- VAR(var_data[,c("d_reer_mex","ca_mex",
                             "gexp_mex","rgdp_mex")], p = 4, type = "trend")

# type = "const", "both", "none"
summary(var_model)

## 4.- SVAR Models
# Recall ordering matters in SVAR models
# First, we construct the coefficient matrices
# Ilzetski Mendoza Vegh (Matrix is lower triangular by Choleski Decomp)
# Coefficient Matrix
A <- matrix(c(1,0,0,0,
               NA,1,0,0,
               NA,NA,1,0,
               NA,NA,NA,1),nrow=4,ncol=4,byrow=TRUE)
# Error Matrix
B <- matrix(c(NA,0,0,0,
               0,NA,0,0,
               0,0,NA,0,
               0,0,0,NA),nrow=4,ncol=4,byrow=TRUE)

# Estimate a SVAR model with the coefficient restriction we gave
svarimv_mex <- SVAR(var_model, Amat = A, Bmat=B)

svarimv_mex

## 5.- Impulse Response Functions
# Now, we use the IRF command to evaluate an orthogonal shock of the 
# government expenditure to the Real GDP up to 12 quarters ahead
t <- 12

# Bootstrap option is not working for me, its a bug on Macs
# Yet, try it in your computer, you should see the IRF
# clearly.
irf_svar <- irf(svarimv_mex,impulse = "gexp_mex", response = c("gexp_mex","rgdp_mex"),
    n.ahead = t, ortho = TRUE, boot=FALSE)
plot(irf_svar)


# Cholesky Decomposition
var_cov <-   svarimv_mex$Sigma.U

chol(var_cov)














