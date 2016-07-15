## gen_sim.data: Generate simulation data
## Input:   my.gamma: 5-component vector used to generate Z|x
##          alph: 5-component vector used to generate A|Z,x
##          mu: 7-component vector used to generate Y|A,Z,x
##             * mu[1] -- our target parameter
##             * mu[2] -- confounding on A
## Output: data.frame of simulation variables (see output for more details)
gen_sim.data <- function(my.gamma,alpha,mu,psi){
  # Set initial conditions
  n <- 1000    # Number of rows
  p <- 4       # Number of x_i variables
  prob <- 0.5  # Probability of binary x_i values

  # Generate xi, zi, A, yi
  xi <- cbind(1,
              matrix(rbinom(n*p,1,prob),n,p))  # xi = (x1,x2,...,x_n) ~ Bernoulli(prob)
  zi <- rbinom(n,1,expit(xi %*% my.gamma))     # zi ~ expit(my.gamma'*x + gamma_0)
  
  #Pr.A.zx <- zi * expitm1(xi %*% alph)   # Pr(A=1|Z,X) ~ expit1m(alph'x + omega_0)
  Pr.A.zx <- zi * tanh(xi %*% alph)
  #Pr.A.zx <- Pr.A.zx + psi
  Pr.A.zx <- Pr.A.zx + xi %*% psi        # Pr(A=1|Z,x) + my.gamma'x
  summary(Pr.A.zx)
  # Stop simulation if we did not generate a probability. 
  if(min(Pr.A.zx) < 0 | max(Pr.A.zx) > 1){
    stop("Oops! Pr(A|Z,x) contains a value outside of [0,1]")
  }

  # Generate A|Z,x and Y|A,Z,x
  A <- rbinom(n,1,Pr.A.zx)   # A ~ Bernoulli(Pr(A=1|Z,X)), column vector
  Yi <- mu[1]*A + mu[2]*{A - Pr.A.zx} + xi %*% mu[3:length(mu)] #mu_{1}*A + mu_{2}*(A - Pr(A|Z,x)) + mu_{3}'(x')
  Yi <- Yi + rnorm(n,0,1)    # E[Y|A,z,x] + eps
  
  ## Output simulation data as data.frame
  ## Note: We preallocate f(z|x), W, E[W|x], R, N as NA values.
  data.frame(
    ones=xi[,1],       # 1:   1
    x=xi[,2:ncol(xi)], # 2-5: x1, x2, x3, x4
    zi=zi,             # 6:   z|x
    A=A,               # 7:   A|z,x
    Yi=Yi,             # 8:   E[Y|A,z,x]
    f.hat.zx=rep(NA,n),    # 9:   f(z|x)
    W=rep(NA,n),       # 10:  W
    E.Wx=rep(NA,n),    # 11:  E[W|X]
    R=rep(NA,n),       # 12:  R
    M=rep(NA,n))       # 13:  N
}