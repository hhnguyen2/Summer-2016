# Newton Raphson to solve for eta, where eta is either my.gamma (usem1==FALSE) or alph (usem1==TRUE)
newtonRaphson <- function(eta_0,data,usem1 = FALSE){
  # Initial setup
  count <- 0
  max_iterations <- 100
  tolerance <- 10^(-12)
  
  # Set eta_t using eta_0 input, and then update eta_t.plus.one. 
  eta_t <- eta_0
  eta_t.plus.one <- eta_t - (grr(eta_t,data,usem1) %*% fr(eta_t,data,usem1))
  
  # Loop until squared diff between eta_t and eta_t.plus.one is minimized
  #      or maximum alloted iterations is reached
  while(count < max_iterations &
        norm(matrix(eta_t.plus.one - eta_t), "F") > tolerance){
    # Increment counter
    count <- count + 1
    
    # Update eta_t
    eta_t <- eta_t.plus.one
    eta_t.plus.one <- eta_t - (grr(eta_t,data,usem1) %*% fr(eta_t,data,usem1))
  }
  
  if(count==max_iterations){
    print("Max iterations reached.")
  }
  # Return: eta_t.plus.one, aka eta_hat
  as.numeric(eta_t.plus.one)
}