# Define function grr that evalutes f'(eta)
grr <- function(eta,sim.data,usem1 = FALSE){
  # INPUTS: eta, with form (eta_0, eta_1, ..., eta_n)
  # OUTPUTS: Inverse Jacobian of f(eta) 
  # ----------------
  # Extract xi
  xi <- extract_xi(sim.data)
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