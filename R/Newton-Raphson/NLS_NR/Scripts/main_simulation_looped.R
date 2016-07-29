## CONDUCT SIMULATION
# Set initial conditions and initialize variables
init_gamma <- rep(0,5) 
init_alpha <- 0
init_B <- rep(0,2)
# Set truth
my.gamma <- c(0.25,-0.25,0.10,-0.45,0.75)  # my.gamma: logistic parameters for Z
alph <- 1
mu <- c(1.62,0.73,0.58,-1.54,-0.65,0.55,1.26) # mu: logistic paramters for A, delta, x
#psi <- c(0.1076,-0.00779,0.0044,-0.0054,0.608)
#psi <- 0.45
#psi <- c(0.57,-0.05,-0.01,-0.08,0.39)       # psi: used to coerce Pr.A.zx between (0,1)

#alph <- c(0.67,-0.34,-0.26,-0.73,0.28)     # alph: logistic parameters for A|Z,X

# Initialize loop conditions
n <- c(500,1000,5000)
Ba_hat <- mapply(rep,0,rep(1000,3))
iter <- 1000
####################
## BEGIN SIMULATION

for(j in 1:3){
  for(i in 1:iter){
    ### Set initial seed
    set.seed(-500000 - (i+j*1000)) # Seed list: (500 iterations)  -- -401001 to -401500
                                   #            (1000 iterations) -- -402001 to -403000
                                   #            (5000 iterations) -- -403001 to -408000
    
    ## Generate simulation data
    sim.data <- gen_sim.data(n[j],my.gamma,alph,mu,psi)
    xi <- extract_xi(sim.data)
    
    ## Approximate alpha_hat by fitting logit Pr(z|X) = alpha'x
    gamma_hat <- as.numeric(glm(zi~x.1+x.2+x.3+x.4,data=sim.data,family=binomial)$coefficients)
  
    ## Generate f.zx & W into preallocated spot in sim.data
    sim.data$f.hat.zx <- gen_f.hat.zx(gamma_hat,sim.data)
    sim.data$W <- gen_W(sim.data)
    
    ## Fit E[W|X] = {exp(alpha'x) - 1} / {exp(alpha'x) + 1} for each i
    # Approximate alpha_hat
    #alpha_hat <- newtonRaphson(init_alpha,sim.data,usem1 = TRUE)
    alpha_hat <- optim(init_alpha,opt_grr_int,data=sim.data,method="Brent",lower=-20,upper=20)$par
    
    ## Generate E[W|X],R, M for each person i
    sim.data$E.Wx <- gen_E.Wx(alpha_hat,sim.data)
    sim.data$R <- sim.data$Yi / sim.data$E.Wx
    
    ## Compute B_bar
    Ba_hat[i,j] <- lm(R~1+zi, weight = 1/f.hat.zx, data = sim.data)$coefficients[2]
    
    ## Check done
    cat(i/iter*100,"% done\n")
  }
}
## END SIMULATION
###################

results.summary.500 <- gen_results(Ba_hat[,1],mu)
results.summary.1000 <- gen_results(Ba_hat[,2],mu)
results.summary.5000 <- gen_results(Ba_hat[,3],mu)

results.summary.500
results.summary.1000
results.summary.5000

par(mar=c(4,5,4,4))
boxplot(Ba_hat, names=c("n = 500","n = 1000","n = 5000"),ylab=expression(hat(beta[a])),main="Boxplot of Simulation Study at n={500,1000,5000}")
abline(h=mu[1],col=2,lty=3)

lm(Yi~A+x.1+x.2+x.3+x.4+x.5+x.6+x.7+x.8+x.9+x.10+x.11+x.12+x.13+x.14+x.15
   ,data=card.data)