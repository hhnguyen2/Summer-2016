### IV_NR_simulation.R
## Desc: This is the code to simulate the binary IV data to test if method works

## Dependencies: Call graphing libraries for histograms
library(ggplot2)

## Define expit and dexpit
# expit: {exp(x)}/{1+exp(x)} = 1/{exp(-x)+1}
expit <- function(x){
  1/{exp(-x)+1}
}

# expitm1: (exp(x)-1)/(exp(x)+1) = 1 - 2/(exp(x)+1)
# Read: "expit minus one"
expitm1 <- function(x){
  1 - 2/{exp(x)+1}
}

# dexpit: expit'(x) = expit(x)(1-expit(x))
dexpit <- function(x){
  expit(x)*{1 - expit(x)}
}

# EXTRACT Xi
# Input: simulation data
# Output: n-by-p matrix: 1,x1,x2,...,xn
extract_xi <- function(sim.data){
  as.matrix(sim.data[,c(1,
                        grep("x.",colnames(sim.data)))])
}

# Linear combination
# Input: x -- n-by-p matrix of n-component vectors
#        eta -- p length vector of coefficients
# Output: Scalar: By-row sum of linear combination: 
#                 eta_1*x_1 + ... + eta_p*x_n
lcomb <- function(x,eta){
  # t(t(x) * eta) computes eta_1*x_1i eta_2*x_2i ... eta_p*x_2n for all i = 1,...,n
  rowSums(t(t(x) * eta))
}


## gen_sim.data: Generate simulation data
## Input:   eta: 5-component vector used to generate Z|x
##          omega: 5-component vector used to generate A|Z,x
##          mu: 7-component vector used to generate Y|A,Z,x
##             * mu[1] -- our target parameter
##             * mu[2] -- confounding on A
## Output: data.frame of simulation variables (see output for more details)
gen_sim.data <- function(){
  # Set initial conditions
  n <- 1000    # Number of rows
  p <- 4       # Number of x_i variables
  prob <- 0.5  # Probability of binary x_i values
  # Set truth
  eta <- c(0.25,-0.25,0.10,-0.45,0.75)  # eta: logistic parameters for Z
  omega <- c(-0.15,0.85,-0.30,0.75,0.60)        # omega: logistic parameters for A|Z,X
  mu <- c(1.62,0,0.58,-1.54,-0.65,0.55,1.26,-0.77) # mu: logistic paramters for A, delta, x
  my.gamma <- c(0.59,-0.16,0.06,-0.13,-0.23)    # gamma: used to coerce Pr.A.zx between (0,1)
  
  # Generate xi, zi, A, yi
  xi <- cbind(1,
              matrix(rbinom(n*p,1,prob),n,p))  # xi = (x1,x2,...,x_n) ~ Bernoulli(prob)
  zi <- rbinom(n,1,expit(lcomb(xi,eta)))       # zi ~ expit(eta'*x + eta_0)
  
  Pr.A.zx <- zi * expitm1(lcomb(xi,omega))   # Pr(A=1|Z,X) ~ expit1m(omega'x + omega_0)
  Pr.A.zx <- Pr.A.zx + lcomb(xi,my.gamma) # Pr(A=1|Z,x) + gamma'x

  # Stop simulation if we did not generate a probability. 
  if(min(Pr.A.zx) < 0 | max(Pr.A.zx) > 1){
    stop("Oops! Pr(A|Z,x) contains a value outside of [0,1]")
  }
  
  summary(Pr.A.zx)
  
  # Generate A|Z,x and Y|A,Z,x
  A <- rbinom(n,1,Pr.A.zx)   # A ~ Bernoulli(Pr(A=1|Z,X)), column vector
  Yi <- mu[1]*A + mu[2]*{A - Pr.A.zx} + lcomb(xi,mu[3:length(mu)]) #mu_{1}*A + mu_{2}*(A - Pr(A|Z,x)) + mu_{3}'(x')
  Yi <- Yi + rnorm(n,0,1)    # E[Y|A,z,x] + eps
  
  ###################
  ## CODE FOR Y CONTINUOUS: Y = mu_1 * a + mu_3 * x + mu_4 * z + eps
  ###################
  
  ###################
  ## CODE FOR Y BINARY: NOT SURE 
  ###################
  
  ## Output simulation data as data.frame
  ## Note: We preallocate f(z|x), W, E[W|x], R, N as NA values.
  data.frame(
    ones=xi[,1],       # 1:   1
    x=xi[,2:ncol(xi)], # 2-5: x1, x2, x3, x4
    zi=zi,             # 6:   z|x
    A=A,               # 7:   A|z,x
    Yi=Yi,             # 8:   E[Y|A,z,x]
    f.zx=rep(NA,n),    # 9:   f(z|x)
    W=rep(NA,n),       # 10:  W
    E.Wx=rep(NA,n),    # 11:  E[W|X]
    R=rep(NA,n),       # 12:  R
    M=rep(NA,n))       # 13:  N
}


## Solver

# Define function fr that evalutes f(eta)
fr <- function(eta,sim.data,usem1 = FALSE){
  # INPUTS: eta, with form (eta_0, eta_1, ..., eta_n)
  # OUTPUTS: f(eta) as column vector
  # ----------------
  # Extract xi
  xi <- extract_xi(sim.data)
  eta <- as.numeric(eta) # coerce data type to ensure operation works

  # Evaluate f(eta); output as column vector
  if(usem1){ # If usem1 == TRUE, then use expitm1, since this is A|Z,x
    col.minus.expit <- sim.data$W - expitm1(lcomb(xi,eta)) #  W - expitm1(eta'x)
    } else{  # Otherwise just use expit, since this is Z|x
    col.minus.expit <- sim.data$zi - expit(lcomb(xi,eta))  # zi - expit(eta'x)
  }
  
  # Output answer as column vector;
  as.matrix(colSums(xi * col.minus.expit))
}

# Define function grr that evalutes f'(eta)
grr <- function(eta,sim.data,usem1 = FALSE){
  # INPUTS: eta, with form (eta_0, eta_1, ..., eta_n)
  # OUTPUTS: Inverse Jacobian of f(eta) 
  # ----------------
  # Extract xi
  xi <- extract_xi(sim.data)
  # Allocate empty answer matrix
  ans <- matrix(0,length(eta),length(eta))
  
  # Iteratively generate Jacobian matrix and add to allocated matrix
  for (i in 1:nrow(sim.data)){
    # Analytic coeff is outer product of xi*t(xi), so take adv of this. 
    add.me <- -1 * dexpit(sum(xi[i,] * eta)) * {xi[i,] %*% t(xi[i,])}
    ans <- ans + add.me
  }
  
  # Evaluate inverse Jacobian. Return NA matrix if singular. 
  if(usem1){ # If using expitm1, dexpit should be scaled by 2
    solve(ans*2) 
  } else{    # Otherwise, just solve as is 
    solve(ans)
  }
  
}

# Newton Raphson to solve for eta
newtonRaphson <- function(eta_0,sim.data,usem1 = FALSE){
  # Initial setup
  count <- 0
  max_iterations <- 100
  tolerance <- 10^(-15)
  
  # Set eta_t using eta_0 input, and then update eta_t.plus.one. 
  eta_t <- eta_0
  eta_t.plus.one <- eta_t - (grr(eta_t,sim.data,usem1) %*% fr(eta_t,sim.data,usem1))
  
  # Loop until squared diff between eta_t and eta_t.plus.one is minimized
  #      or maximum alloted iterations is reached
  while(count < max_iterations &
        norm(matrix(eta_t.plus.one - eta_t), "F") > tolerance){
    # Increment counter
    count <- count + 1
    
    # Update eta_t
    eta_t <- eta_t.plus.one
    eta_t.plus.one <- eta_t - (grr(eta_t,sim.data,usem1) %*% fr(eta_t,sim.data,usem1))
  }
  
  # Return: eta_t.plus.one, aka eta_hat
  as.numeric(eta_t.plus.one)
}

#Generate f(z|x) from zi's and x's
gen_f.zx <- function(alpha_hat,sim.data){
  # Extract xi and zi
  xi <- extract_xi(sim.data)
  zi <- sim.data$zi
  # Allocate empty output matrix
  f.zx <- rep(NA,1000)
  
  # Vectorized: f(z|x) = expit(alpha_hat'x)     if z==1, 
  #                    = 1 - expit(alpha_hat'x) otherwise
  f.zx[which(zi==1)] <- expit(lcomb(xi[which(zi==1),],alpha_hat))
  f.zx[which(zi==0)] <- 1 - expit(lcomb(xi[which(zi==0),],alpha_hat))
  f.zx
}

# Generate W's from f(z|x)'s 
gen_W <- function(sim.data){
  # Extract f(z|x), A, and zi
  A <- sim.data$A
  f.zx <- sim.data$f.zx
  zi <- sim.data$zi
  
  # Compute W
  {A*{(-1)^(1-zi)}} / {f.zx}
}

gen_E.Wx <- function(omega_hat,sim.data){
  # Extract xi
  xi <- extract_xi(sim.data)
  
  # Generate E.Wx
  expitm1(lcomb(xi,omega_hat))
}

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