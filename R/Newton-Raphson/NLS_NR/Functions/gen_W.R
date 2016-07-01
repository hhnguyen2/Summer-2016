# Generate W's from f(z|x)'s 
gen_W <- function(sim.data){
  # Extract f(z|x), A, and zi
  A <- sim.data$A
  f.zx <- sim.data$f.zx
  zi <- sim.data$zi
  
  # Compute W
  {A*{(-1)^(1-zi)}} / {f.zx}
}