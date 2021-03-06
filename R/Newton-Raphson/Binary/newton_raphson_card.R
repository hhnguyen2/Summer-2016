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


## Solver

# Define function fr that evalutes f(eta)
fr <- function(eta,my.data,use1m = FALSE){
  # INPUTS: eta, with form (eta_0, eta_1, ..., eta_n)
  # OUTPUTS: f(eta) as column vector
  # ----------------
  # Extract xi
  xi <- as.matrix(my.data[,c(1,
                             grep("x.",colnames(my.data)))])
  eta <- as.numeric(eta) # coerce data type to ensure operation works
  
  # Evaluate f(eta); eta[1] is eta_0; eta[-1] is eta_1,...,eta_4
  if(use1m){ # If use1m == TRUE, then use expit1m
    zi.minus.pi <- my.data$W - expit1m(rowSums(t(t(xi) * eta))) # zi - xi*eta
  } else{    # Otherwise just use e
    zi.minus.pi <- my.data$zi - expit(rowSums(t(t(xi) * eta))) # zi - xi*eta
  }
  
  # Output answer as column vector; [,1:5] is 1,x1,...,x4
  as.matrix(colSums(xi * zi.minus.pi))
}

# Define function grr that evalutes f'(eta)
grr <- function(eta,my.data,use1m = FALSE){
  # INPUTS: eta, with form (eta_0, eta_1, ..., eta_n)
  # OUTPUTS: Inverse Jacobian of f(eta) 
  # ----------------
  # Extract xi
  xi <- as.matrix(my.data[,c(1,
                             grep("x.",colnames(my.data)))])
  # Allocate empty answer matrix
  ans <- matrix(0,length(eta),length(eta))
  
  # Iteratively generate Jacobian matrix and add to allocated matrix
  for (i in 1:nrow(my.data)){
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
newtonRaphson <- function(eta_0,my.data,use1m = FALSE){
  # Initial setup
  count <- 0
  max_iterations <- 100
  tolerance <- 10^(-15)
  
  # Set eta_t using eta_0 input, and then update eta_t.plus.one. 
  eta_t <- eta_0
  eta_t.plus.one <- eta_t - (grr(eta_t,my.data,use1m) %*% fr(eta_t,my.data,use1m))
  
  # Loop until squared diff between eta_t and eta_t.plus.one is minimized
  #      or maximum alloted iterations is reached
  while(count < max_iterations &
        norm(matrix(eta_t.plus.one - eta_t), "F") > tolerance){
    # Increment counter
    count <- count + 1
    
    # Update eta_t
    eta_t <- eta_t.plus.one
    eta_t.plus.one <- eta_t - (grr(eta_t,my.data,use1m) %*% fr(eta_t,my.data,use1m))
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
gen_f.zx <- function(alpha_hat,my.data){
  # Extract xi and zi
  xi <- my.data[,c(1,
                   grep("x.",colnames(my.data)))] # 1,x1, ..., xn
  zi <- my.data$zi
  # Allocate empty answer matrix
  f.zx <- rep(NA,nrow(my.data))
  
  # Vectorized: f(z|x) = expit(alpha_hat'x) if z==1, 1 - expit(alpha_hat'x) otherwise
  f.zx[which(zi==1)] <- expit(rowSums(t(t(xi[which(zi==1),])*alpha_hat)))
  f.zx[which(zi==0)] <- 1 - expit(rowSums(t(t(xi[which(zi==0),])*alpha_hat)))
  f.zx
}

# Generate W's from f(z|x)'s 
gen_W <- function(my.data){
  # Extract f(z|x), A, and zi
  A <- my.data$A
  f.zx <- my.data$f.zx
  zi <- my.data$zi
  
  # Compute W
  {A*{(-1)^(1-zi)}} / {f.zx}
}

gen_E.Wx <- function(theta,my.data){
  # Extract xi
  xi <- my.data[,c(1,
                   grep("x.",colnames(my.data)))] # 1,x1, ..., xn
  
  # Generate E.Wx
  expit1m(rowSums(t(t(xi)*theta)))
}

### Read data 

setwd("C:/Users/jimmy/Desktop/git/Summer-2016/R/Newton-Raphson/Binary/")
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

## Run estimator
# Set initial conditions
initial_eta <- rep(0,16)
initial_theta<- rep(0,16)

# Approximate alpha_hat
alpha_hat <- newtonRaphson(initial_eta,my.data,use1m = FALSE)

# Generate f.zx & W into preallocated spot in my.data
my.data$f.zx <- gen_f.zx(alpha_hat,my.data)
my.data$W <- gen_W(my.data)

## Fit E[W|X] = {exp(theta'x) - 1} / {exp(theta'x) + 1} for each i
# Approximate theta_hat
theta_hat <- newtonRaphson(initial_theta,my.data,use1m = TRUE)
# Generate E[W|X],R, M for each person i
my.data$E.Wx <- gen_E.Wx(theta_hat,my.data)
my.data$R <- my.data$Yi / my.data$E.Wx
my.data$M <- 1 / (my.data$f.zx)

# Define Z_bar, M_bar, R_bar
Z_bar <- cbind(my.data$ones,my.data$zi)
M_bar <- diag(my.data$M)
R_bar <- my.data$R

# Compute B_bar
B_hat <- solve(t(Z_bar) %*% M_bar %*% Z_bar) %*% (t(Z_bar) %*% M_bar %*% R_bar)
B_hat