#gmi-sim-model0.R

#modelstring0はパラメータをn,y0,p,b,を全て推定
modelstring = "model{
for( i in 1 : N ) {
  y[i] ~ dlnorm(m, 1/(s^2))
}
m <- log(y0)+n*log(1-b)+log((1+b)/(1-b))*n*p
s <- sqrt(n*p*(1-p)*(log((1+b)/(1-b)))^2)
a ~ dunif(1,50)
n ~ dpois(a)
c ~ dunif(1,50) 
y0 ~ dpois(c)
p ~ dbeta(1,1)
b ~ dbeta(1,1)}"
# close quote for modelstring

parameters = c("p","b","n","a","m","s","c","y0")