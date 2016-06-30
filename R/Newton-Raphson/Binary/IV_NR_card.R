### IV_NR_card.R
## Desc: Runs IV method on Card data with binary Z, X, Y

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
  omega <- c(0.3,0.6,0.75,0.80,0.25)    # omega : logistic parameters for A|Z,X
  mu <- c(1,0.30,0.55,-0.95,-0.29,-0.7,-0.97) # mu: logistic paramters for A, delta, x
  my.gamma <- c(0.01,0.05,0.02,0.01,0.03) # gamma: used to coerce Pr.A.zx between (0,1)
  
  # Generate xi, zi, A, yi
  # Note: t(t(xi) * eta) = sweep(xi,2,eta,'*') = x_1*eta_1 + x_2*eta_2 + ... 
  xi <- matrix(rbinom(n*p,1,prob),n,p) # x1, x2,...,xn
  zi <- rbinom(n,1,expit(rowSums(t(t(xi) * eta[-1])) + eta[1]))       # zi ~ expit(eta'*x + eta_0)
  
  Pr.A.zx <- zi * expitm1(rowSums(t(t(xi) * omega[-1])) + omega[1])   # Pr(A=1|Z,X) ~ expit1m(omega'x + omega_0)
  Pr.A.zx <- Pr.A.zx + rowSums(t(t(xi) * my.gamma[-1])) + my.gamma[1] # Pr(A=1|Z,x) + gamma'x

  # Stop simulation if we did not generate a probability. 
  if(min(Pr.A.zx) < 0 | max(Pr.A.zx) > 1){
    stop("Oops! Pr(A|Z,x) contains a value outside of [0,1]")
  }
  
  # Generate A|Z,x and Y|A,Z,x
  A <- rbinom(n,1,Pr.A.zx)   # A ~ Bernoulli(Pr(A=1|Z,X)), column vector
  Yi <- mu[1]*A + mu[2]*{A - Pr.A.zx} + mu[3] + rowSums(t(t(xi)*mu[4:length(mu)])) #mu_{1}*A + mu_{2}*(A - Pr(A|Z,x)) + mu_{3}'(x')
  Yi <- Yi + rnorm(n,0,1)    # E[Y|A,z,x] + eps
  
  ## Output simulation data as data.frame
  ## Note: We preallocate f(z|x), W, E[W|x], R, N as NA values.
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
fr <- function(eta,sim.data,usem1 = FALSE){
  # INPUTS: eta, with form (eta_0, eta_1, ..., eta_n)
  # OUTPUTS: f(eta) as column vector
  # ----------------
  # Extract xi
  xi <- as.matrix(sim.data[,c(1,
                                grep("x.",colnames(sim.data)))])
  eta <- as.numeric(eta) # coerce data type to ensure operation works

  # Evaluate f(eta); output as column vector
  if(usem1){ # If usem1 == TRUE, then use expitm1, since this is A|Z,x
    col.minus.expit <- sim.data$W - expitm1(rowSums(t(t(xi) * eta))) #  W - expitm1(eta'x)
    } else{  # Otherwise just use expit, since this is Z|x
    col.minus.expit <- sim.data$zi - expit(rowSums(t(t(xi) * eta)))  # zi - expit(eta'x)
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
  xi <- as.matrix(sim.data[,c(1,
                                grep("x.",colnames(sim.data)))])
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
  xi <- sim.data[,c(1,
                      grep("x.",colnames(sim.data)))] # 1,x1,...,xn
  zi <- sim.data$zi
  # Allocate empty output matrix
  f.zx <- rep(NA,1000)
  
  # Vectorized: f(z|x) = expit(alpha_hat'x)     if z==1, 
  #                    = 1 - expit(alpha_hat'x) otherwise
  f.zx[which(zi==1)] <- expit(rowSums(t(t(xi[which(zi==1),])*alpha_hat)))
  f.zx[which(zi==0)] <- 1 - expit(rowSums(t(t(xi[which(zi==0),])*alpha_hat)))
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
  xi <- sim.data[,c(1,
                      grep("x.",colnames(sim.data)))] # 1,x1, ..., xn
  
  # Generate E.Wx
  expitm1(rowSums(t(t(xi)*omega_hat)))
}

#### RUN METHOD
## Read data 
setwd("C:/Users/Jimmy/Desktop/git/Summer-2016/R/Newton-Raphson/Binary")
data <- read.csv("Data/card.csv")

x.1 <- data$black
x.2 <- data$south
x.3 <- data$smsa
x.4 <- apply(data[,12:20],1,function(x) which(x==1))
x.5 <- data$smsa66
x.6 <- data$exper
x.7 <- data$fatheduc;    x.8 <- as.numeric(is.na(x.7))
x.7[is.na(x.7)] <- mean(x.7,na.rm=TRUE)     # Card(1995)'s imputation
x.9 <- data$motheduc;     x.10 <- as.numeric(is.na(x.9))
x.9[is.na(x.9)] <- mean(x.9,na.rm=TRUE)
x.11 <- data$momdad14
x.12 <- data$sinmom14
x.13 <- data$step14
x.14 <- data$kww;    x.15 <- as.numeric(is.na(x.14))
x.14[is.na(x.14)] <- mean(x.14,na.rm=TRUE)

my.data <- data.frame(ones=rep(1,nrow(data)),
                      x.1=x.1,
                      x.2=x.2,
                      x.3=x.3,
                      x.4=x.4,
                      x.5=x.5,
                      x.6=x.6,
                      x.7=x.7,
                      x.8=x.8,
                      x.9=x.9,
                      x.10=x.10,
                      x.11=x.11,
                      x.12=x.12,
                      x.13=x.13,
                      x.14=x.14,
                      x.15=x.15,
                      zi=data$nearc4,
                      A=as.numeric(data$educ > 12),
                      Yi=as.numeric(data$wage > median(data$wage)),
                      f.zx=rep(NA,nrow(data)), # f(z|x)
                      W=rep(NA,nrow(data)),    # W
                      E.Wx=rep(NA,nrow(data)), # E[W|X]
                      R=rep(NA,nrow(data)),    # R
                      M=rep(NA,nrow(data)))    # M

## my estimators
iter <- 1000
N <- nrow(my.data)
Ba_hat_est <- rep(NA,iter)


for(i in 1:iter){
  ## Bootstrap
  index = sample(1:N,replace=TRUE)
  data.btstrap=my.data[index,]
  
  ## Run estimator
  # Set initial conditions
  initial_eta <- rep(0,16)
  initial_theta<- rep(0,16)
  
  # Approximate alpha_hat
  alpha_hat <- newtonRaphson(initial_eta,data.btstrap,usem1 = FALSE)
  
  # Generate f.zx & W into preallocated spot in data.btstrap
  data.btstrap$f.zx <- gen_f.zx(alpha_hat,data.btstrap)
  data.btstrap$W <- gen_W(data.btstrap)
  
  ## Fit E[W|X] = {exp(theta'x) - 1} / {exp(theta'x) + 1} for each i
  # Approximate theta_hat
  theta_hat <- newtonRaphson(initial_theta,data.btstrap,usem1 = TRUE)
  # Generate E[W|X],R, M for each person i
  data.btstrap$E.Wx <- gen_E.Wx(theta_hat,data.btstrap)
  data.btstrap$R <- data.btstrap$Yi / data.btstrap$E.Wx
  data.btstrap$M <- 1 / (data.btstrap$f.zx)
  
  # Define Z_bar, M_bar, R_bar
  Z_bar <- cbind(data.btstrap$ones,data.btstrap$zi)
  M_bar <- diag(data.btstrap$M)
  R_bar <- data.btstrap$R
  
  # Compute B_bar
  B_hat <- solve(t(Z_bar) %*% M_bar %*% Z_bar) %*% (t(Z_bar) %*% M_bar %*% R_bar)
  Ba_hat_est[i] <- B_hat[2,1]
  
  cat(100*i/iter, "% done", "\n")
}