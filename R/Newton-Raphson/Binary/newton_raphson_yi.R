## Dependencies: Call libraries
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
##      - $alpha: original alpha values used to generate simulation data
##      - $sim.Matrix: n*p size matrix containing in order: yi, 1, xi
gen_sim.data <- function(){
  n <- 1000
  p <- 4
  prob <- 0.5
  alpha <- c(0.25,-0.25,0.10,-0.45,0.75)  # alpha: The target truth/parameter
  
  randMatrix <- matrix(rbinom(n*p,1,prob),n,p) # x1, x2,...,xn
  # Note: double transpose operation is just multiplying each column by corresponding elements of vector... but fast
  sd <- 0.25 # SD of gaussian noise, epsilon
  yi <- expit1m(rowSums(t(t(randMatrix) * alpha[-1])) + alpha[1]) + rnorm(n,0,sd)
  # Ensure that -1 <= y <= 1
  while(sum(abs(yi) > 1) > 0){
    yi <- expit1m(rowSums(t(t(randMatrix) * alpha[-1])) + alpha[1]) + rnorm(n,0,sd)
  }
  
  sim.Matrix <- cbind(yi,           # yi's
                      rep(1,n),     # x0
                      randMatrix)   # x1,x2,...,xn
  
  return(list(alpha=alpha,sim.Matrix=sim.Matrix))
}


## Solver
# Define function fr that evalutes f(alpha)
fr <- function(alpha, sim.Matrix){
  # INPUTS: alpha, with form (alpha_0, alpha_1, ..., alpha_n)
  # OUTPUTS: f(alpha) as column vector
  # ----------------
  # Extract x1 to xn
  xi <- sim.Matrix[,3:ncol(sim.Matrix)] # x1, ..., x4
  
  # Evaluate f(alpha); alpha[1] is alpha_0; alpha[-1] is alpha_1...
  yi.minus.pi <- sim.Matrix[,1] - expit1m(alpha[1] + rowSums(t(t(xi) * alpha[-1]))) # yi - xi*alpha
  
  # Output answer as column vector
  as.matrix(colSums(sim.Matrix[,2:6] * yi.minus.pi))
}

# Define function grr that evalutes f'(alpha)
grr <- function(alpha, sim.Matrix){
  # INPUTS: alpha, with form (alpha_0, alpha_1, ..., alpha_n)
  # OUTPUTS: Inverse Jacobian of f(alpha) 
  # ----------------
  # Allocate empty answer matrix, should be 5x5
  ans <- matrix(0,length(alpha),length(alpha))

  # Iteratively generate Jacobian matrix and add to allocated matrix
  for (i in 1:nrow(sim.Matrix)){
    # Analytic coeff is outer product of xi*t(xi), so take adv of this. 
    add.me <- sim.Matrix[i,2:6] %*% t(sim.Matrix[i,2:6])
    add.me <- -2 * add.me * dexpit(alpha[1] + sum(sim.Matrix[i,3:6] * alpha[2:5]))
    
    ans <- ans + add.me
  }
  
  # Evaluate inverse Jacobian. Return NA matrix if singular. 
  solve(ans)
}

# Newton Raphson to solve for alpha
newtonRaphson <- function(alpha_0,real_alpha,sim.Matrix){
  # Initial setup
  count <- 0
  max_iterations <- 100
  diverged <- FALSE 
  
  # Set alpha_t using alpha_0 input, and then update alpha_t.plus.one. 
  alpha_t <- alpha_0
  alpha_t.plus.one <- alpha_t - (grr(alpha_t,sim.Matrix) %*% fr(alpha_t,sim.Matrix))
  
  # Loop until squared diff between alpha_t and alpha_t.plus.one is minimized
  #      or maximum alloted iterations is reached
  while(count < max_iterations &&
        norm(matrix(alpha_t.plus.one - alpha_t), "F") > .Machine$double.eps){
    # Increment counter
    count <- count + 1
    
    # Update alpha_t
    alpha_t <- alpha_t.plus.one
    alpha_t.plus.one <- alpha_t - (grr(alpha_t,sim.Matrix) %*% fr(alpha_t,sim.Matrix))
  }
  
  # Flag as diverged if solver reaches max iterations
  if (count >= max_iterations){
    diverged <- TRUE
  }
  
  # Return: alpha_t.plus.one, converge or diverge? 
  data.frame(alpha_hat=matrix(alpha_t.plus.one,1,5),diverged=diverged,iterations=count)
}

## CONDUCT SIMULATION
# Set initial conditions
iter <- 100
initial_alpha <- c(0,0,0,0,0)
results <- data.frame(alpha_hat=matrix(0,iter,5),diverged=rep(NA,iter),iterations=rep(NA,iter))
#set.seed(187536841) # mean(runif(1000,-10103203,390439843))
#set.seed(146150199) # mean(runif(1000,-101032023,390439843))
#set.seed(-316198416) # mean(runif(1000,-1010320223,390439843))

# Run simulation 
for (i in 1:iter){
  sim.data <- gen_sim.data()
  results[i,] <- newtonRaphson(initial_alpha,sim.data$alpha,sim.data$sim.Matrix)
  
  cat(100*i/iter, "% done", "\n")
}

## Results: avg. parameter, bias, variance, mean^2 error
true_alpha <- sim.data$alpha
divergers <- which(results$diverged == TRUE)
alpha_hat <- results[,1:5]
alpha_hat.avg <- apply(alpha_hat,2,mean)
alpha_hat.bias <- (alpha_hat.avg - true_alpha)
alpha_hat.pcbias <- (alpha_hat.bias / true_alpha) * 100
alpha_hat.mean2 <- alpha_hat.bias^2 + apply(alpha_hat,2,var)

## Tabulated results
results.summary <- data.frame(true_alpha=true_alpha,
                              avg_value=alpha_hat.avg,
                              bias=alpha_hat.bias,
                              pc_bias=alpha_hat.pcbias,
                              mean2_err=alpha_hat.mean2)

## Make histogram, ah1p means alpha_hat.1 plot, etc.
ah1p <- ggplot(alpha_hat, aes(x=alpha_hat.1)) +
  geom_histogram(bins = 30) +
  geom_vline(aes(xintercept=mean(alpha_hat.1)),
             color="red", linetype="dashed", size=1, alpha=.5)
ah2p <- ggplot(alpha_hat, aes(x=alpha_hat.2)) +
  geom_histogram(bins = 30) +
  geom_vline(aes(xintercept=mean(alpha_hat.2)),
             color="red", linetype="dashed", size=1, alpha=.5)
ah3p <- ggplot(alpha_hat, aes(x=alpha_hat.3)) +
  geom_histogram(bins = 30) +
  geom_vline(aes(xintercept=mean(alpha_hat.3)),
             color="red", linetype="dashed", size=1, alpha=.5)
ah4p <- ggplot(alpha_hat, aes(x=alpha_hat.4)) +
  geom_histogram(bins = 30) +
  geom_vline(aes(xintercept=mean(alpha_hat.4)),
             color="red", linetype="dashed", size=1, alpha=.5)
ah5p <- ggplot(alpha_hat, aes(x=alpha_hat.5)) +
  geom_histogram(bins = 30) +
  geom_vline(aes(xintercept=mean(alpha_hat.5)),
             color="red", linetype="dashed", size=1, alpha=.5)

multiplot(ah1p,ah2p,ah3p,ah4p,ah5p, cols = 3) # plot all on one graph