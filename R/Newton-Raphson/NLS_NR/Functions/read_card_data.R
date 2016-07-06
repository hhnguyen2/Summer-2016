## Read data 
read_card_data <- function(){
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
  
  data.frame(ones=rep(1,nrow(data)),
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
             f.hat.zx=rep(NA,nrow(data)), # f(z|x)
             W=rep(NA,nrow(data)),    # W
             E.Wx=rep(NA,nrow(data)), # E[W|X]
             R=rep(NA,nrow(data)),    # R
             M=rep(NA,nrow(data)))    # M
}
