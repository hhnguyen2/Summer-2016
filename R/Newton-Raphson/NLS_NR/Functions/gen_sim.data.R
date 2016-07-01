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
  mu <- c(1.62,0,0.58,-1.54,-0.65,0.55,1.26) # mu: logistic paramters for A, delta, x
  my.gamma <- c(0.59,-0.16,0.06,-0.13,-0.23)    # gamma: used to coerce Pr.A.zx between (0,1)
  
  # Generate xi, zi, A, yi
  xi <- cbind(1,
              matrix(rbinom(n*p,1,prob),n,p))  # xi = (x1,x2,...,x_n) ~ Bernoulli(prob)
  zi <- rbinom(n,1,expit(xi %*% eta))       # zi ~ expit(eta'*x + eta_0)
  
  Pr.A.zx <- zi * expitm1(xi %*% omega)   # Pr(A=1|Z,X) ~ expit1m(omega'x + omega_0)
  Pr.A.zx <- Pr.A.zx + xi %*% my.gamma # Pr(A=1|Z,x) + gamma'x
  
  # Stop simulation if we did not generate a probability. 
  if(min(Pr.A.zx) < 0 | max(Pr.A.zx) > 1){
    stop("Oops! Pr(A|Z,x) contains a value outside of [0,1]")
  }
  
  summary(Pr.A.zx)
  
  # Generate A|Z,x and Y|A,Z,x
  A <- rbinom(n,1,Pr.A.zx)   # A ~ Bernoulli(Pr(A=1|Z,X)), column vector
  Yi <- mu[1]*A + mu[2]*{A - Pr.A.zx} + xi %*% mu[3:length(mu)] #mu_{1}*A + mu_{2}*(A - Pr(A|Z,x)) + mu_{3}'(x')
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