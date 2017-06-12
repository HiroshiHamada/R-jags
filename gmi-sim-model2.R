# gmi-sim-model2.R: 事前分布を母数の近くに設定
# modelstring0はパラメータをn,y0,p,b,を全て推定
modelstring = "model{
for( i in 1 : N ) {
  y[i] ~ dlnorm(m, 1/(s^2))
}
m <- log(y0)+n*log(1-b)+log((1+b)/(1-b))*n*p
s <- sqrt(n*p*(1-p)*(log((1+b)/(1-b)))^2)
n <- 50

c ~ dunif(15,25) 
y0 ~ dpois(c)

p ~ dbeta(1,1)
b ~ dbeta(1,1)}"
# close quote for modelstring

parameters = c("p","b","m","s","y0")

#真のパラメータ設定  y1, b0, p0, n0
# n.player=5000 #データ数
# y0.data=20; b.data=0.1; p.data=0.6; n.data=50

# a ~ dunif(45,55)
# n ~ dpois(a)