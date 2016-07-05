####################
## BEGIN SIMULATION

### Set initial seed
set.seed(-23018349)

## Generate simulation data
sim.data <- gen_sim.data()

## Approximate alpha_hat by fitting logit Pr(z|X) = alpha'x
gamma_hat <- newtonRaphson(initial_gamma,sim.data,usem1 = FALSE)

## Generate f.zx & W into preallocated spot in sim.data
sim.data$f.zx <- gen_f.zx(gamma_hat,sim.data)
sim.data$W <- gen_W(sim.data)

## Fit E[W|X] = {exp(alpha'x) - 1} / {exp(alpha'x) + 1} for each i
# Approximate alpha_hat
alpha_hat <- newtonRaphson(initial_alpha,sim.data,usem1 = TRUE)

## Generate E[W|X],R, M for each person i
sim.data$E.Wx <- gen_E.Wx(alpha_hat,sim.data)
sim.data$R <- sim.data$Yi / sim.data$E.Wx
sim.data$M <- 1 / (sim.data$f.zx)

## Define Z_bar, M_bar, R_bar
Z_bar <- cbind(sim.data$ones,sim.data$zi)
M_bar <- diag(sim.data$M)
R_bar <- sim.data$R

## Compute B_bar
B_hat <- solve(t(Z_bar) %*% M_bar %*% Z_bar) %*% (t(Z_bar) %*% M_bar %*% R_bar)

## END SIMULATION
###################