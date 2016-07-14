n <- 1000    # Number of rows
p <- 4       # Number of x_i variables
prob <- 0.5  # Probability of binary x_i values
# Set truth
eta <<- c(0.25,-0.25,0.10,-0.45,0.75)  # eta: logistic parameters for Z
omega <- c(0.8,0.8,0.8,0.8,-0.3)    # omega : logistic parameters for A|Z,X
mu <<- c(1,0.30,0.55,-0.95,-0.29,-0.7,-0.97) # mu: logistic paramters for A, delta, x
#my.gamma <<- c(0.78,-0.11,0.14,0.26,0.15) # gamma: used to coerce Pr.A.zx between (0,1)
# 0.7827138 -0.1139980 -0.1378345 -0.2649706  0.1541452
# Generate xi, zi, A, yi
# Note: t(t(xi) * eta) = sweep(xi,2,eta,'*') = x_1*eta_1 + x_2*eta_2 + ... 
xi <- matrix(rbinom(n*p,1,prob),n,p) # x1, x2,...,xn
Pr.zx <<- expit(rowSums(t(t(xi) * eta[-1])) + eta[1])
zi <- rbinom(n,1,Pr.zx)           # zi ~ expit(eta'*x + eta_0)

#omega <- round(runif(5,-1,1),2)
Pr.A.zx <- zi * tanh(rowSums(t(t(xi) * omega[-1])) + omega[1])   # Pr(A=1|Z,X) ~ expit1m(omega'x + omega_0)
det_fact <- {abs(max(Pr.A.zx)) + abs(min(Pr.A.zx))} / 2
my.gamma <- as.numeric(glm((0.60 - Pr.A.zx) ~ xi)$coefficients)
#my.gamma <- c(0.05,0,0,0,0)
Pr.A.zx <- Pr.A.zx + rowSums(t(t(xi) * my.gamma[-1])) + my.gamma[1] #                               + gamma'x
omega
my.gamma
det_fact
summary(Pr.A.zx)