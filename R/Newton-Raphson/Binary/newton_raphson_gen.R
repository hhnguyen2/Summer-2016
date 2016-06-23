## Dependencies: Call graphing libraries for histograms
library(ggplot2)
library(grid)
library(gridExtra)

# Multiple plot function (Winston Chang, R Cookbook)
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

## Define expit and dexpit
# expit: (e^x)/(1+e^x) = 1/(exp(-x)+1)
expit <- function(x){
  1/{exp(-x)+1}
}

# expit1m: (e^(x)-1)/(1+e^x) = 1 - 2/(exp(x)+1)
expit1m <- function(x){
  1 - 2/{exp(x)+1}
}

# dexpit: expit'(x), or expit(x)(1-expit(x))
dexpit <- function(x){
  expit(x)*{1 - expit(x)}
}


## gen_sim.data: Simulate data
## Output: a list of 2 objects: 
##      - $eta: original eta values used to generate simulation data
##      - $sim.Matrix: n*p size matrix containing in order: zi, 1, xi
gen_sim.data <- function(eta = eta,omega = omega,mu = mu){
  n <- 1000
  p <- 4
  prob <- 0.5

  # Generate xi, zi, A, yi
  xi <- matrix(rbinom(n*p,1,prob),n,p) # x1, x2,...,xn
  zi <- rbinom(n,1,expit(rowSums(t(t(xi) * eta[-1])) + eta[1]))  # zi ~ expit(eta'*x + eta_0)
  Pr.A.zx <- zi * expit1m(rowSums(t(t(xi) * omega[-1])) + omega[1]) # Pr(A=1|Z,X)
  #my.gamma <- as.numeric(glm(0.1 * Pr.A.zx ~ xi)$coefficients)
  my.gamma <- c(0.01,0.05,0.02,0.01,0.03)
  Pr.A.zx <- Pr.A.zx + rowSums(t(t(xi) * my.gamma[-1])) + my.gamma[1]

  # Stop simulation if we did not generate a probability. 
  if(min(Pr.A.zx) < 0 | max(Pr.A.zx) > 1){
    stop("ABANDON SHIP!!!!~ WE GOTTA GET OUTTA HEREE....")
  }

  
  A <- rbinom(n,1,Pr.A.zx) # A ~ Bernoulli(Pr(A=1|Z,X)), column vector
  Yi <- mu[1]*A + mu[2]*{A - Pr.A.zx} + mu[3] + rowSums(t(t(xi)*mu[4:length(mu)])) #mu[3] = mu_3*1 (+ mu_4*x1...)
  Yi <- Yi + rnorm(n,0,0.15) # E[Y|A,z,x] + eps
  
  ## Output simulation data as matrix
  data.frame(
    ones=rep(1,n),  # 1:   1
    xi=xi,          # 2-5: x1, x2, x3, x4
    zi=zi,          # 6:   z|x
    A=A,            # 7:   A|z,x
    Yi=Yi,          # 8:   E[Y|A,z,x]
    f.zx=rep(NA,n), # 9:   f(z|x)
    W=rep(NA,n),    # 10:  W
    E.Wx=rep(NA,n), # 11:  E[W|X]
    R=rep(NA,n),    # 12:  R
    M=rep(NA,n))    # 13:  N
}


## Solver

# Define function fr that evalutes f(eta)
fr <- function(eta,sim.Matrix,use1m = FALSE){
  # INPUTS: eta, with form (eta_0, eta_1, ..., eta_n)
  # OUTPUTS: f(eta) as column vector
  # ----------------
  # Extract xi
  xi <- as.matrix(sim.Matrix[,c(1,
                                  grep("xi",colnames(sim.Matrix)))])
  eta <- as.numeric(eta) # coerce data type to ensure operation works

  # Evaluate f(eta); eta[1] is eta_0; eta[-1] is eta_1,...,eta_4
  if(use1m){ # If use1m == TRUE, then use expit1m
    zi.minus.pi <- sim.Matrix$W - expit1m(rowSums(t(t(xi) * eta))) # zi - xi*eta
  } else{    # Otherwise just use e
    zi.minus.pi <- sim.Matrix$zi - expit(rowSums(t(t(xi) * eta))) # zi - xi*eta
  }
  
  # Output answer as column vector; [,1:5] is 1,x1,...,x4
  as.matrix(colSums(xi * zi.minus.pi))
}

# Define function grr that evalutes f'(eta)
grr <- function(eta,sim.Matrix,use1m = FALSE){
  # INPUTS: eta, with form (eta_0, eta_1, ..., eta_n)
  # OUTPUTS: Inverse Jacobian of f(eta) 
  # ----------------
  # Extract xi
  xi <- as.matrix(sim.Matrix[,c(1,
                                  grep("xi",colnames(sim.Matrix)))])
  # Allocate empty answer matrix
  ans <- matrix(0,length(eta),length(eta))
  
  # Iteratively generate Jacobian matrix and add to allocated matrix
  for (i in 1:nrow(sim.Matrix)){
    # Analytic coeff is outer product of xi*t(xi), so take adv of this. 
    add.me <- -1 * dexpit(sum(xi[i,] * eta)) * {xi[i,] %*% t(xi[i,])}
    ans <- ans + add.me
  }
  
  # Evaluate inverse Jacobian. Return NA matrix if singular. 
  if(use1m){ # If using expit1m, dexpit should be scaled by 2
    solve(ans*2) 
  } else{    # Otherwise, just solve as is 
    solve(ans)
  }
  
}

# Newton Raphson to solve for eta
newtonRaphson <- function(eta_0,sim.Matrix,use1m = FALSE){
  # Initial setup
  count <- 0
  max_iterations <- 100
  diverged <- FALSE 
  tolerance <- 10^(-15)
  
  # Set eta_t using eta_0 input, and then update eta_t.plus.one. 
  eta_t <- eta_0
  eta_t.plus.one <- eta_t - (grr(eta_t,sim.Matrix,use1m) %*% fr(eta_t,sim.Matrix,use1m))
  
  # Loop until squared diff between eta_t and eta_t.plus.one is minimized
  #      or maximum alloted iterations is reached
  while(count < max_iterations &
        norm(matrix(eta_t.plus.one - eta_t), "F") > tolerance){
    # Increment counter
    count <- count + 1
    
    # Update eta_t
    eta_t <- eta_t.plus.one
    eta_t.plus.one <- eta_t - (grr(eta_t,sim.Matrix,use1m) %*% fr(eta_t,sim.Matrix,use1m))
  }
  
  # Flag as diverged if solver reaches max iterations
  if (count >= max_iterations){
    diverged <- TRUE
  }
  
  # Return: eta_t.plus.one, converge or diverge? 
  #data.frame(eta_hat=matrix(eta_t.plus.one,1,length(eta_t.plus.one)),diverged=diverged,iterations=count)
  as.numeric(eta_t.plus.one)
}

#Generate f(z|x) from zi's and x's
gen_f.zx <- function(alpha,sim.Matrix){
  # Extract xi and zi
  xi <- sim.Matrix[,c(1,
                      grep("xi",colnames(sim.Matrix)))] # 1,x1, ..., xn
  zi <- sim.Matrix$zi
  # Allocate empty answer matrix
  f.zx <- rep(NA,1000)
  
  # Vectorized: f(z|x) = expit(alpha'x) if z==1, 1 - expit(alpha'x) otherwise
  f.zx[which(zi==1)] <- expit(rowSums(t(t(xi[which(zi==1),])*alpha)))
  f.zx[which(zi==0)] <- 1 - expit(rowSums(t(t(xi[which(zi==0),])*alpha)))
  f.zx
}

# Generate W's from f(z|x)'s 
gen_W <- function(sim.Matrix){
  # Extract f(z|x), A, and zi
  A <- sim.Matrix$A
  f.zx <- sim.Matrix$f.zx
  zi <- sim.Matrix$zi
  
  # Compute W
  {A*{(-1)^(1-zi)}} / {f.zx}
}

gen_E.Wx <- function(omega,sim.Matrix){
  # Extract xi
  xi <- sim.Matrix[,c(1,
                      grep("xi",colnames(sim.Matrix)))] # 1,x1, ..., xn
  
  # Generate E.Wx
  expit1m(rowSums(t(t(xi)*omega)))
}

## CONDUCT SIMULATION
# Set truth
eta <- c(0.25,-0.25,0.10,-0.45,0.75)  # eta: logistic parameters for Z
omega <- c(0.3,0.6,0.75,0.80,0.25) # omega : logistic parameters for A|Z,X
mu <- c(1,0.30,0.55,-0.95,-0.29,-0.7,-0.97) # mu: logistic paramters for A + delta + x

# Set initial conditions
iter <- 1000
initial_eta <- rep(0,length(eta)) 
initial_omega <- rep(0,length(omega))
B_hat <- matrix(NA,iter,2) # Bo, Ba 
#set.seed(-5664498)
#set.seed(1199449)
set.seed(-23018349)

#################
## BEGIN SIMULATION

for(i in 1:iter){
  ## Fit logit Pr(z|X) = alpha'x
  sim.data <- gen_sim.data(eta,omega,mu)
  
  # Approximate alpha_hat
  alpha_hat <- newtonRaphson(initial_eta,sim.data,use1m = FALSE)
  
  # Generate f.zx & W into preallocated spot in sim.data
  sim.data$f.zx <- gen_f.zx(alpha_hat,sim.data)
  sim.data$W <- gen_W(sim.data)
  
  ## Fit E[W|X] = {exp(omega'x) - 1} / {exp(omega'x) + 1} for each i
  # Approximate omega_hat
  omega_hat <- newtonRaphson(initial_omega,sim.data,use1m = TRUE)
  # Generate E[W|X],R, M for each person i
  sim.data$E.Wx <- gen_E.Wx(omega_hat,sim.data)
  sim.data$R <- sim.data$Yi / sim.data$E.Wx
  sim.data$M <- 1 / (sim.data$f.zx)
  
  # Define Z_bar, M_bar, R_bar
  Z_bar <- cbind(sim.data$ones,sim.data$zi)
  M_bar <- diag(sim.data$M)
  R_bar <- sim.data$R
  
  # Compute B_bar
  B_hat[i,] <- solve(t(Z_bar) %*% M_bar %*% Z_bar) %*% (t(Z_bar) %*% M_bar %*% R_bar)

  # Update progress
  cat(100*i/iter, "% done", "\n")
}

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

## Make histogram, ah1p means B_hat.1 plot, etc.
Ba_hat_hist <- data.frame(Ba_hat=Ba_hat[Ba_hat > Ba_hat.avg - 3*sd(Ba_hat) & Ba_hat < Ba_hat.avg + 3*sd(Ba_hat)])
Ba_hat_outliers <- Ba_hat[Ba_hat < Ba_hat.avg - 3*sd(Ba_hat) & Ba_hat > Ba_hat.avg + 3*sd(Ba_hat)]

Ba_hat_hist_full <- data.frame(Ba_hat=Ba_hat)

ggplot(Ba_hat_hist, aes(x=Ba_hat)) +
  geom_histogram(bins = 20) +
  geom_vline(aes(xintercept=mean(Ba_hat)),
             color="red", linetype="dashed", size=1, alpha=.5) +
  ggtitle("Histogram of Ba_hat (no outliers)")

ggplot(Ba_hat_hist_full, aes(x=Ba_hat)) +
  geom_histogram(bins = 20) +
  geom_vline(aes(xintercept=mean(Ba_hat)),
             color="red", linetype="dashed", size=1, alpha=.5) +
  ggtitle("Histogram of Ba_hat (all values)")