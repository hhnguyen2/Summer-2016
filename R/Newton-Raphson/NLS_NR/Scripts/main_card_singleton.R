## Run estimator
# Set initial conditions
init_gamma <- rep(0,16)
init_alpha <- 0
init_B <- rep(0,2)

## Approximate alpha_hat by fitting logit Pr(z|X) = alpha'x
gamma_glm <- as.numeric(glm(zi~x.1+x.2+x.3+x.4+x.5+x.6+x.7+x.8+x.9+x.10+x.11+x.12+x.13+x.14+x.15
                            ,data=card.data,family=binomial)$coefficients)
gamma_hat <- gamma_glm     

#gamma_opt
#gamma_nr
gamma_glm

# Generate f.zx & W into preallocated spot in card.data
card.data$f.hat.zx <- gen_f.hat.zx(gamma_hat,card.data)
card.data$W <- gen_W(card.data)

## Fit E[W|X] = {exp(alpha'x) - 1} / {exp(alpha'x) + 1} for each i
# Approximate alpha_hat
alpha_opt <- optim(init_alpha,opt_grr_int,data=card.data,method="Brent",lower=-10,upper=10)$par
#alpha_nr <- newtonRaphson(init_alpha,card.data,usem1 = TRUE)
alpha_hat <- alpha_opt

alpha_opt
#alpha_nr

# Generate E[W|X],R, M for each person i
card.data$E.Wx <- gen_E.Wx(alpha_hat,card.data)
card.data$R <- card.data$Yi / card.data$E.Wx
#card.data$M <- 1 / (card.data$f.hat.zx)

# Define Z_bar, M_bar, R_bar
#Z_bar <- cbind(card.data$ones,card.data$zi)
#M_bar <- diag(card.data$M)
#R_bar <- card.data$R

# Compute B_bar
#B_closed <- solve(t(Z_bar) %*% M_bar %*% Z_bar) %*% (t(Z_bar) %*% M_bar %*% R_bar)
B_lm <- lm(R~1+zi, weight = 1/f.hat.zx, data = card.data)$coefficients
#B_opt <- optim(init_B,opt_hrrr,data=card.data)$par

#B_closed
B_lm
#B_opt