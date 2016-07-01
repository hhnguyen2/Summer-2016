### IV_simulation.R
## Desc: This is the code to simulate the binary IV data to test if method works

## Dependencies: Call graphing libraries for histograms
library(ggplot2)

## Evaluation functions
source("Functions/expit.R")      # Evals exp(x)/(1+exp(x))
source("Functions/expitm1.R")    # Evals {exp(x)-1} / {exp(x)+1}
source("Functions/dexpit.R")     # Evals (expit(x))/(1-expit(x))

## Solver
source("Functions/fr.R")            # Evaluates (X)*(zi - expit(xi %*% eta))
source("Functions/grr.R")           # Evaluates gradient of fr
source("Functions/newtonRaphson.R") # Newton-Raphson root solver

## Data Generation Functions
source("Functions/gen_sim.data.R")  # Generates simulation data: 1,xi,zi,A,Y
source("Functions/extract_xi.R")    # Extracts x_i matrix from data
source("Functions/gen_f.zx.R")      # Generates f(z|x) from alpha_hat estimate
source("Functions/gen_W.R")         # Generates W = {A}^{1-z} / {f(z=1|x)}
source("Functions/gen_E.Wx.R")      # Generates E[W|X] = expitm1(xi %*% omega_hat)


## CONDUCT SIMULATION
# Set initial conditions and initialize variables
iter <- 1000 
initial_eta <- rep(0,5) 
initial_omega <- rep(0,5)
B_hat <- matrix(NA,iter,2) # Bo, Ba 

####################
## BEGIN SIMULATION

### Set initial seed
set.seed(-23018349)

### Run simulation for (iter) times.
for(i in 1:iter){
  ## Generate simulation data
  sim.data <- gen_sim.data()
  
  ## Approximate alpha_hat by fitting logit Pr(z|X) = alpha'x
  alpha_hat <- newtonRaphson(initial_eta,sim.data,usem1 = FALSE)
  
  ## Generate f.zx & W into preallocated spot in sim.data
  sim.data$f.zx <- gen_f.zx(alpha_hat,sim.data)
  sim.data$W <- gen_W(sim.data)
  
  ## Fit E[W|X] = {exp(omega'x) - 1} / {exp(omega'x) + 1} for each i
  # Approximate omega_hat
  omega_hat <- newtonRaphson(initial_omega,sim.data,usem1 = TRUE)
  
  ## Generate E[W|X],R, M for each person i
  sim.data$E.Wx <- gen_E.Wx(omega_hat,sim.data)
  sim.data$R <- sim.data$Yi / sim.data$E.Wx
  sim.data$M <- 1 / (sim.data$f.zx)
  
  ## Define Z_bar, M_bar, R_bar
  Z_bar <- cbind(sim.data$ones,sim.data$zi)
  M_bar <- diag(sim.data$M)
  R_bar <- sim.data$R
  
  ## Compute B_bar
  B_hat[i,] <- solve(t(Z_bar) %*% M_bar %*% Z_bar) %*% (t(Z_bar) %*% M_bar %*% R_bar)
  
  ## Update progress (so we don't get bored)
  cat(100*i/iter, "% done", "\n")
}
## END SIMULATION
#################

## Results: avg. parameter, bias, variance, mean^2 error
true_mu <- mu[1]
Ba_hat <- B_hat[,2]
Ba_hat.avg <- mean(Ba_hat)
Ba_hat.bias <- Ba_hat.avg - true_mu
Ba_hat.pcbias <- (Ba_hat.bias / true_mu) * 100
Ba_hat.var <- var(Ba_hat)
Ba_hat.mean2 <- Ba_hat.bias^2 + Ba_hat.var

## Tabulated results
results.summary <- data.frame(mu=true_mu,
                              mean=Ba_hat.avg,
                              bias=Ba_hat.bias,
                              pc_bias=Ba_hat.pcbias,
                              var=Ba_hat.var,
                              mean2_err=Ba_hat.mean2)

## Make histogram
## Ba_hat_hist_full has all n values
## Ba_hat_hist is cleaned out outliers
## Ba_hat_outliers contains the outliers
Ba_hat_hist_full <- data.frame(Ba_hat=Ba_hat)
Ba_hat_hist <- data.frame(Ba_hat=Ba_hat[Ba_hat > Ba_hat.avg - 3*sd(Ba_hat) & Ba_hat < Ba_hat.avg + 3*sd(Ba_hat)])
Ba_hat_outliers <- Ba_hat[Ba_hat < Ba_hat.avg - 3*sd(Ba_hat) | Ba_hat > Ba_hat.avg + 3*sd(Ba_hat)]

## Plot histogram WITHOUT outliers
ggplot(Ba_hat_hist, aes(x=Ba_hat)) +
  geom_histogram(bins = 20) +
  geom_vline(aes(xintercept=mean(Ba_hat)),
             color="red", linetype="dashed", size=1, alpha=.5) +
  ggtitle("Histogram of Ba_hat (no outliers)")

## Plot histogram WITH outliers
ggplot(Ba_hat_hist_full, aes(x=Ba_hat)) +
  geom_histogram(bins = 20) +
  geom_vline(aes(xintercept=mean(Ba_hat)),
             color="red", linetype="dashed", size=1, alpha=.5) +
  ggtitle("Histogram of Ba_hat (all values)")