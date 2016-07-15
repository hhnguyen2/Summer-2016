#### RUN METHOD
# Initialize variables
init_gamma <- rep(0,16)
init_alpha <- rep(0,16)
init_B <- rep(0,2)


# Set loop conditions
iter <- 10
Ba_hat_est <- rep(NA,iter)
N <- nrow(card.data)

for(i in 1:iter){
  ## Set seed
  set.seed(140000 + i) # Seed ranges from 130001-131000
  
  ## Bootstrap
  index = sample(1:N,replace=TRUE)
  data.btstrap=card.data[index,]
  
  # Approximate gamma_hat
  gamma_hat <- as.numeric(glm(zi~x.1+x.2+x.3+x.4+x.5+x.6+x.7+x.8+x.9+x.10+x.11+x.12+x.13+x.14+x.15
  ,data=data.btstrap,family=binomial)$coefficients)
  #gamma_hat <- newtonRaphson(init_gamma,data.btstrap,usem1=FALSE)
  #gamma_hat <- optim(init_gamma,opt_fr,data=data.btstrap)$par
  
  # Generate f.zx & W into preallocated spot in data.btstrap
  data.btstrap$f.hat.zx <- gen_f.hat.zx(gamma_hat,data.btstrap)
  data.btstrap$W <- gen_W(data.btstrap)
  
  ## Fit E[W|X] = {exp(alpha'x) - 1} / {exp(alpha'x) + 1} for each i
  # Approximate alpha_hat
  #alpha_hat <- newtonRaphson(init_alpha,data.btstrap,usem1 = TRUE)
  #alpha_hat <- optim(init_alpha,opt_grr,data=data.btstrap)$par
  alpha_hat <- rep(0,16)
  alpha_hat[1] <- alpha_hat[1] + 1
  
  # Generate E[W|X],R, M for each person i
  data.btstrap$E.Wx <- gen_E.Wx(alpha_hat,data.btstrap)
  data.btstrap$R <- data.btstrap$Yi / data.btstrap$E.Wx
  #data.btstrap$M <- 1 / (data.btstrap$f.hat.zx)
  
  # Define Z_bar, M_bar, R_bar
  #Z_bar <- cbind(data.btstrap$ones,data.btstrap$zi)
  #M_bar <- diag(data.btstrap$M)
  #R_bar <- data.btstrap$R
  
  # Compute B_bar
  #B_hat <- solve(t(Z_bar) %*% M_bar %*% Z_bar) %*% (t(Z_bar) %*% M_bar %*% R_bar)
  #B_hat <- optim(init_B,opt_hrrr,data=data.btstrap)$par
  B_hat <- lm(R~1+zi, weight = 1/f.hat.zx, data = data.btstrap)$coefficients
  Ba_hat_est[i] <- B_hat[2]

  cat(100*i/iter, "% done", "\n")
}

hist(Ba_hat_est, main="Histogram of Beta_a using glm,NR,lm", xlab="Beta_a")
Ba_hat_hist_full <- data.frame(Ba_hat=Ba_hat_est)
ggplot(Ba_hat_hist_full, aes(x=Ba_hat)) +
  geom_histogram(bins = 20) +
  geom_vline(aes(xintercept=mean(Ba_hat)),
             color="red", linetype="dashed", size=1, alpha=.5) +
  ggtitle(expression("Histogram of "*hat(beta[a])))