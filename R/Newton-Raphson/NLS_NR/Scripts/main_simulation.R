## CONDUCT SIMULATION
# Set initial conditions and initialize variables
init_gamma <- rep(0,5) 
init_alpha <- rep(0,5)
init_B <- rep(0,2)
# Set truth
my.gamma <- c(0.25,-0.25,0.10,-0.45,0.75)  # my.gamma: logistic parameters for Z
alph <- 0.35*c(.1,.1,.1,.1,-0.5)                       # alph: logistic parameters for A|Z,X
mu <- c(1.62,0,0.58,-1.54,-0.65,0.55,1.26) # mu: logistic paramters for A, delta, x
psi <- 0.6
#psi <- c(0.25,0.17,0.13,0.38,-0.16)        # psi: used to coerce Pr.A.zx between (0,1)

####################
## BEGIN SIMULATION

### Set initial seed
i <- 550
set.seed(-23018349 - i)

## Generate simulation data
sim.data <- gen_sim.data(my.gamma,alph,mu,psi)
xi <- extract_xi(sim.data)

## Approximate alpha_hat by fitting logit Pr(z|X) = alpha'x
gamma_opt <- optim(init_gamma,opt_fr,data=sim.data)$par
gamma_nr <- newtonRaphson(init_gamma,sim.data,usem1 = FALSE)
gamma_glm <- as.numeric(glm(zi~x.1+x.2+x.3+x.4,data=sim.data,family=binomial)$coefficients)

gamma_hat <- gamma_nr
## Generate f.zx & W into preallocated spot in sim.data
sim.data$f.hat.zx <- gen_f.hat.zx(gamma_hat,sim.data)
sim.data$W <- gen_W(sim.data)

## Fit E[W|X] = {exp(alpha'x) - 1} / {exp(alpha'x) + 1} for each i
# Approximate alpha_hat
alpha_opt <- optim(init_alpha,opt_grr,data=sim.data)$par
alpha_nr <- newtonRaphson(init_alpha,sim.data,usem1 = TRUE)

alpha_hat <- alpha_nr

## Generate E[W|X],R, M for each person i
sim.data$E.Wx <- gen_E.Wx(alpha_hat,sim.data)
sim.data$R <- sim.data$Yi / sim.data$E.Wx
sim.data$M <- 1 / (sim.data$f.hat.zx)

## Define Z_bar, M_bar, R_bar
Z_bar <- cbind(sim.data$ones,sim.data$zi)
M_bar <- diag(sim.data$M)
R_bar <- sim.data$R

## Compute B_bar
B_closed <- solve(t(Z_bar) %*% M_bar %*% Z_bar) %*% (t(Z_bar) %*% M_bar %*% R_bar)
B_lm <- lm(R~1+zi, weight = 1/f.hat.zx, data = sim.data)$coefficients
B_opt <- optim(init_B,opt_hrrr,data=sim.data)$par

## Output
## Diag
my.gamma
gamma_nr
gamma_glm
gamma_opt
alph
alpha_nr
alpha_opt
B_closed
B_lm
B_opt

mu


## END SIMULATION
###################