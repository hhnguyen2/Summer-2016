### IV_simulation.R
## Desc: This is the code to simulate the binary IV data to test if method works

## Dependencies: Call graphing libraries for histograms
library(ggplot2)

## Evaluation functions
source("Functions/expit.R")      # Evals exp(x)/(1+exp(x))
source("Functions/dexpit.R")     # Evals (expit(x))*(1-expit(x))

## Newton-Raphson functions
source("Functions/fr.R")            # Evaluates (X)*(zi - expit(xi %*% gamma))
source("Functions/grr.R")           # Evaluates gradient of fr
source("Functions/newtonRaphson.R") # Newton-Raphson root solver

## Optim functions
source("Functions/opt_fr.R")     # Evaluates squared of (X)*(zi - expit(xi %*% gamma))
source("Functions/opt_grr.R")    # Evaluates squared of (X)*(zi - expit(xi %*% ))
source("Functions/opt_hrrr.R")   # Evaluates squared of (1,Zi)*(Ri-B0-BaZi)(1/(f.hat.zx))
source("Functions/opt_grr_int.R")# Evaluates squared of (X)*(W - tanh(alpha_0))

## Data Generation Functions
source("Functions/gen_sim.data.R")  # Generates simulation data: 1,xi,zi,A,Y
source("Functions/gen_f.hat.zx.R")  # Generates f(z|x) from alpha_hat estimate
source("Functions/gen_W.R")         # Generates W = {A}^{1-z} / {f(z=1|x)}
source("Functions/gen_E.Wx.R")      # Generates E[W|X] = expitm1(xi %*% alpha_hat)
source("Functions/extract_xi.R")    # Extracts x_i matrix from data
source("Functions/gen_results.R")   # Generates mean, bias, %bias, var, and MSE of Ba estimates

## Read in Card Data
source("Functions/read_card_data.R")

## Run simulation
card.data <- read_card_data()
#source("Scripts/main_card_analysis.R",echo=TRUE)
source("Scripts/main_simulation_looped.R",echo=TRUE)
#source("Reference/main_card_analysis.R",echo=TRUE)