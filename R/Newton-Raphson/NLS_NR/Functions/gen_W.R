# Generate W's from f(z|x)'s 
gen_W <- function(data){
  # Extract f(z|x), A, and zi
  A <- data$A
  f.hat.zx <- data$f.hat.zx
  zi <- data$zi
  
  # Compute W
  {A*{(-1)^(1-zi)}} / {f.hat.zx}
}