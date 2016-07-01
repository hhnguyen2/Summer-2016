# expit: evaluate {exp(x)}/{1+exp(x)} = 1/{exp(-x)+1}
expit <- function(x){
  1/{exp(-x)+1}
}