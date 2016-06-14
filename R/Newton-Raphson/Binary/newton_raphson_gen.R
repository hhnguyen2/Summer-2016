## Dependencies: Call graphing libraries for histograms
library(ggplot2)
library(grid)
library(gridExtra)

# Multiple plot function (Winston Chang)
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
gen_sim.data <- function(){
  n <- 1000
  p <- 4
  prob <- 0.5
  eta <- c(0.25,-0.25,0.10,-0.45,0.75)  # eta: logistic parameters for Z
  omega <- c(-0.3,0.6,0.75,0.80,-0.25,0.33) # omega : logistic parameters for A|Z,X
  mu <- c(0.14,-0.50,0.18,0.16,-0.87,0.90,0.20) # mu: logistic paramters for A + delta + x
  
  xi <- matrix(rbinom(n*p,1,prob),n,p) # x1, x2,...,xn
  zi <- rbinom(n,1,expit(rowSums(t(t(xi) * eta[-1])) + eta[1]))  # zi ~ expit(eta'*x + eta_0)
  Pr.A.zx <- expit(rowSums(t(t(cbind(zi,xi)) * omega[-1])) + omega[1]) # Pr(A=1|Z,X)
  A <- rbinom(n,1,Pr.A.zx) # A ~ Bernoulli(Pr(A=1|Z,X)), column vector
  Yi <- mu[1]*A + mu[2]*{A - Pr.A.zx} + mu[3] + rowSums(t(mu[4:length(mu)]*t(xi))) #mu[3] = mu_3*1 (+ mu_4*x1...)
  Yi <- Yi + rnorm(n,0,0.1) # E[Y|A,z,x] + eps
  sim.Matrix <- cbind(rep(1,n), # 1
                      xi,       # x1, x2, x3, x4
                      zi,       # z|x
                      A,        # A|z,x
                      Yi)       # E[Y|A,z,x]
  
  return(list(eta=eta,sim.Matrix=sim.Matrix))
}


## Solver
# Define function fr that evalutes f(eta)
fr <- function(eta, sim.Matrix){
  # INPUTS: eta, with form (eta_0, eta_1, ..., eta_n)
  # OUTPUTS: f(eta) as column vector
  # ----------------
  # Extract x1 to xn
  xi <- sim.Matrix[,3:ncol(sim.Matrix)] # x1, ..., x4
  
  # Evaluate f(eta); eta[1] is eta_0; eta[-1] is eta_1,...,eta_4
  zi.minus.pi <- sim.Matrix[,1] - expit1m(eta[1] + rowSums(t(t(xi) * eta[-1]))) # zi - xi*eta
  
  # Output answer as column vector
  as.matrix(colSums(sim.Matrix[,2:6] * zi.minus.pi))
}

# Define function grr that evalutes f'(eta)
grr <- function(eta, sim.Matrix){
  # INPUTS: eta, with form (eta_0, eta_1, ..., eta_n)
  # OUTPUTS: Inverse Jacobian of f(eta) 
  # ----------------
  # Allocate empty answer matrix, should be 5x5
  ans <- matrix(0,length(eta),length(eta))
  
  # Iteratively generate Jacobian matrix and add to allocated matrix
  for (i in 1:nrow(sim.Matrix)){
    # Analytic coeff is outer product of xi*t(xi), so take adv of this. 
    add.me <- sim.Matrix[i,2:6] %*% t(sim.Matrix[i,2:6])
    add.me <- -2 * add.me * dexpit(eta[1] + sum(sim.Matrix[i,3:6] * eta[2:5]))
    
    ans <- ans + add.me
  }
  
  # Evaluate inverse Jacobian. Return NA matrix if singular. 
  solve(ans)
}

# Newton Raphson to solve for eta
newtonRaphson <- function(eta_0,real_eta,sim.Matrix){
  # Initial setup
  count <- 0
  max_iterations <- 100
  diverged <- FALSE 
  
  # Set eta_t using eta_0 input, and then update eta_t.plus.one. 
  eta_t <- eta_0
  eta_t.plus.one <- eta_t - (grr(eta_t,sim.Matrix) %*% fr(eta_t,sim.Matrix))
  
  # Loop until squared diff between eta_t and eta_t.plus.one is minimized
  #      or maximum alloted iterations is reached
  while(count < max_iterations &&
        norm(matrix(eta_t.plus.one - eta_t), "F") > .Machine$double.eps){
    # Increment counter
    count <- count + 1
    
    # Update eta_t
    eta_t <- eta_t.plus.one
    eta_t.plus.one <- eta_t - (grr(eta_t,sim.Matrix) %*% fr(eta_t,sim.Matrix))
  }
  
  # Flag as diverged if solver reaches max iterations
  if (count >= max_iterations){
    diverged <- TRUE
  }
  
  # Return: eta_t.plus.one, converge or diverge? 
  data.frame(eta_hat=matrix(eta_t.plus.one,1,5),diverged=diverged,iterations=count)
}

## CONDUCT SIMULATION
# Set initial conditions
iter <- 1000
initial_eta <- c(0,0,0,0,0)
results <- data.frame(eta_hat=matrix(0,iter,5),diverged=rep(NA,iter),iterations=rep(NA,iter))
#set.seed(187536841) # mean(runif(1000,-10103203,390439843))
#set.seed(146150199) # mean(runif(1000,-101032023,390439843))
#set.seed(-316198416) # mean(runif(1000,-1010320223,390439843))
set.seed(982839)

# Run simulation 
for (i in 1:iter){
  sim.data <- gen_sim.data()
  results[i,] <- newtonRaphson(initial_eta,sim.data$eta,sim.data$sim.Matrix)
  
  cat(100*i/iter, "% done", "\n")
}

## Results: avg. parameter, bias, variance, mean^2 error
true_eta <- sim.data$eta
divergers <- which(results$diverged == TRUE)
eta_hat <- results[,1:5]
eta_hat.avg <- apply(eta_hat,2,mean)
eta_hat.bias <- (eta_hat.avg - true_eta)
eta_hat.pcbias <- (eta_hat.bias / true_eta) * 100
eta_hat.mean2 <- eta_hat.bias^2 + apply(eta_hat,2,var)

## Tabulated results
results.summary <- data.frame(true_eta=true_eta,
                              avg_value=eta_hat.avg,
                              bias=eta_hat.bias,
                              pc_bias=eta_hat.pcbias,
                              mean2_err=eta_hat.mean2)

## Make histogram, eh1p means eta_hat.1 plot, etc.
eh1p <- ggplot(eta_hat, aes(x=eta_hat.1)) +
  geom_histogram(bins = 30) +
  geom_vline(aes(xintercept=mean(eta_hat.1)),
             color="red", linetype="dashed", size=1, eta=.5)
eh2p <- ggplot(eta_hat, aes(x=eta_hat.2)) +
  geom_histogram(bins = 30) +
  geom_vline(aes(xintercept=mean(eta_hat.2)),
             color="red", linetype="dashed", size=1, eta=.5)
eh3p <- ggplot(eta_hat, aes(x=eta_hat.3)) +
  geom_histogram(bins = 30) +
  geom_vline(aes(xintercept=mean(eta_hat.3)),
             color="red", linetype="dashed", size=1, eta=.5)
eh4p <- ggplot(eta_hat, aes(x=eta_hat.4)) +
  geom_histogram(bins = 30) +
  geom_vline(aes(xintercept=mean(eta_hat.4)),
             color="red", linetype="dashed", size=1, eta=.5)
eh5p <- ggplot(eta_hat, aes(x=eta_hat.5)) +
  geom_histogram(bins = 30) +
  geom_vline(aes(xintercept=mean(eta_hat.5)),
             color="red", linetype="dashed", size=1, eta=.5)

multiplot(eh1p,eh2p,eh3p,eh4p,eh5p, cols = 3) # plot all on one graph