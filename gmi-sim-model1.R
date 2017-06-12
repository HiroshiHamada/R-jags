#gmi-sim-model1.R: estimate only p. 
modelstring = "model{
for( i in 1 : N ) {
  y[i] ~ dlnorm(m, 1/(s^2))
}
  y0 <- y0.data
  n <- n.data
  b <- b.data
  m <- log(y0)+n*log(1-b)+log((1+b)/(1-b))*n*p
  s <- sqrt(n*p*(1-p)*(log((1+b)/(1-b)))^2)
#prior for hyper parameters
  p ~ dbeta(1,1)
}"
# close quote for modelstring
# b ~ dbeta(1,1)

parameters = c("p","m","s")