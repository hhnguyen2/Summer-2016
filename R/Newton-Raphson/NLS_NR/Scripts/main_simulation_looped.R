## CONDUCT SIMULATION
# Set initial conditions and initialize variables
init_gamma <- rep(0,5) 
init_alpha <- rep(0,5)
init_B <- rep(0,2)
# Set truth
my.gamma <- c(0.25,-0.25,0.10,-0.45,0.75)  # my.gamma: logistic parameters for Z
alph <- 0.40*c(1,.1,.1,.1,-3) 
mu <- c(1.62,0.73,0.58,-1.54,-0.65,0.55,1.26) # mu: logistic paramters for A, delta, x
psi <- c(0.1076,-0.00779,0.0044,-0.0054,0.608)
#psi <- 0.45
#psi <- c(0.57,-0.05,-0.01,-0.08,0.39)       # psi: used to coerce Pr.A.zx between (0,1)

#alph <- c(0.67,-0.34,-0.26,-0.73,0.28)     # alph: logistic parameters for A|Z,X

# Initialize loop conditions
iter <- 1000
Ba_hat <- rep(0,iter) # answer for each loop i
####################
## BEGIN SIMULATION

for(i in 1:iter){
  ### Set initial seed
  set.seed(-21000 - i) # Seed ranges from -23018350 to -23019349
  
  ## Generate simulation data
  sim.data <- gen_sim.data(my.gamma,alph,mu,psi)
  xi <- extract_xi(sim.data)
  
  ## Approximate alpha_hat by fitting logit Pr(z|X) = alpha'x
  gamma_hat <- newtonRaphson(init_gamma,sim.data,usem1 = FALSE)

  ## Generate f.zx & W into preallocated spot in sim.data
  sim.data$f.hat.zx <- gen_f.hat.zx(gamma_hat,sim.data)
  sim.data$W <- gen_W(sim.data)
  
  ## Fit E[W|X] = {exp(alpha'x) - 1} / {exp(alpha'x) + 1} for each i
  # Approximate alpha_hat
  alpha_hat <- newtonRaphson(init_alpha,sim.data,usem1 = TRUE)

  ## Generate E[W|X],R, M for each person i
  sim.data$E.Wx <- gen_E.Wx(alpha_hat,sim.data)
  sim.data$R <- sim.data$Yi / sim.data$E.Wx
  sim.data$M <- 1 / (sim.data$f.hat.zx)
  
  ## Define Z_bar, M_bar, R_bar
  Z_bar <- cbind(sim.data$ones,sim.data$zi)
  M_bar <- diag(sim.data$M)
  R_bar <- sim.data$R
  
  ## Compute B_bar
  B_hat <- solve(t(Z_bar) %*% M_bar %*% Z_bar) %*% (t(Z_bar) %*% M_bar %*% R_bar)
  Ba_hat[i] <- B_hat[2,1]
  
  ## Check done
  cat(i/iter*100,"% done\n")
}
## END SIMULATION
###################

## Results: avg. parameter, bias, variance, mean^2 error
true_mu <- mu[1]
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

Ba_hat_hist_full <- data.frame(Ba_hat=Ba_hat)
ggplot(Ba_hat_hist_full, aes(x=Ba_hat)) +
  geom_histogram(bins = 20) +
  geom_vline(aes(xintercept=mean(Ba_hat)),
             color="red", linetype="dashed", size=1, alpha=.5) +
  ggtitle("Histogram of Ba_hat (all values)")

Ba_hat_hist <- data.frame(Ba_hat=Ba_hat[Ba_hat > mean(Ba_hat) - 3*sd(Ba_hat) & Ba_hat < mean(Ba_hat) + 3*sd(Ba_hat)])
Ba_hat_outliers <- Ba_hat[Ba_hat < mean(Ba_hat) - 3*sd(Ba_hat) | Ba_hat > mean(Ba_hat) + 3*sd(Ba_hat)]

## Plot histogram WITHOUT outliers
ggplot(Ba_hat_hist, aes(x=Ba_hat)) +
  geom_histogram(bins = 20) +
  geom_vline(aes(xintercept=mean(Ba_hat)),
             color="red", linetype="dashed", size=1, alpha=.5) +
  ggtitle("Histogram of Ba_hat (no outliers)")