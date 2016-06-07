## Define expit and dexpit
# expit: e^x/(1+e^x)
expit <- function(x){
  exp(x)/{1+exp(x)}
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
    eta <- c(0.25,-0.25,0.10,-0.45,0.75)  # eta: The target truth/parameter
    
    randMatrix <- matrix(rbinom(n*p,1,prob),n,p) # x1, x2,...,xn
    # Note: double transpose operation is just multiplying each column by corresponding elements of vector... but fast
    zi <- rbinom(n,1,expit(rowSums(t(t(randMatrix) * eta[2:length(eta)])) + eta[1])) # zi ~ expit(eta'*x + eta_0)
    sim.Matrix <- cbind(zi,        # zi's
                     rep(1,n),     # x0
                     randMatrix)   # x1,x2,...,xn
    
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

    # Evaluate f(eta); eta[1] is eta_0
    zi.minus.pi <- sim.Matrix[,1] - expit(eta[1] + rowSums(t(t(xi) * eta[2:length(eta)]))) # zi - xi*eta
    
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
    #ones.matrix <- matrix(1,length(eta),length(eta))
    
    # Iteratively generate Jacobian matrix and add to allocated matrix
    for (i in 1:nrow(sim.Matrix)){
        # Analytic coeff is outer product of xi*t(xi), so take adv of this. 
        add.me <- sim.Matrix[i,2:6] %*% t(sim.Matrix[i,2:6])
        add.me <- -1 * add.me * dexpit(eta[1] + sum(sim.Matrix[i,3:6] * eta[2:5]))
        
        ans <- ans + add.me
    }
    
    # Evaluate inverse Jacobian. Return NA matrix if singular. 
    #tryCatch(solve(ans), error=function(e) matrix(NA,5,5))
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
    
    if (count >= max_iterations){
      diverged <- TRUE
      print("Another one bites the dust.")
    }
  
    # Return: eta_t.plus.one, converge or diverge? 
   data.frame(eta_hat=matrix(eta_t.plus.one,1,5),diverged=diverged,iterations=count)
}

#newtonRaphson(myEta,c(eta_0,eta),sim.Matrix)
iter <- 1000
initial_eta <- c(0,0,0,0,0)
results <- data.frame(eta_hat=matrix(0,iter,5),diverged=rep(NA,iter),iterations=rep(NA,iter))
set.seed(187536841) # mean(runif(1000,-10103203,390439843))

ptm <- proc.time() #let's time this slow code...

for (i in 1:iter){
  
    sim.data <- gen_sim.data()
    results[i,] <- newtonRaphson(initial_eta,sim.data$eta,sim.data$sim.Matrix)
    
    cat(100*i/iter, "% done", "\n")
}

## Results: avg. parameter, bias, variance, mean^2 error
sim.data$eta
divergers <- which(results$diverged == TRUE)
eta_hat <- results[,1:5]
eta_hat.avg <- apply(eta_hat,2,mean)
eta_hat.bias <- (eta_hat.avg - sim.data$eta)
eta_hat.mean2 <- eta_hat.bias^2 + apply(eta_hat,2,var)

## Tabulated results
data.frame(true_eta=sim.data$eta,eta_hat.avg,eta_hat.bias,eta_hat.mean2)

roc.time() - ptm