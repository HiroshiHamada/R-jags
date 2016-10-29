
#------------------------------------------------------------------------------

# Jags-Ymet-Xnom1grp-MlogNormal-Script.R
# April 19, 2016. John K. Kruschke.
# Requires auxiliary R scripts that accompany the book,
# Kruschke, J. K. (2015). Doing Bayesian Data Analysis, Second Edition:
# A Tutorial with R, JAGS, and Stan. Academic Press / Elsevier.

# そもそもgmcをMCMC推定できているのかどうかを人工データでチェック

# graphics.off()
# rm(list=ls(all=TRUE))

source("DBDA2E-utilities.R")
fileNameRoot="jags-gmc-sim"
graphFileType="png"

#------------------------------------------------------------------------------
# THE DATA.
# Generate a person's income under the given condition: 
set.seed(47406)

#dataの分散が足りないので，p0を確率変数化する．

gmi<-function(y1, b0, p0, n0){
  #y1:初期資本, b0:投資割合, p0:成功確率, n0:反復回数
  y = y1
for(i in 1:n0){
  if(runif(1) < p0 + rnorm(1,0,0.1)){
  y <- y + y*(b0+rnorm(1,0,0.01))}
else{y <- y - y*(b0+rnorm(1,0,0.1))} }#for end
y
}


#真のパラメータ設定  y1, b0, p0, n0
m=2000; y1=20; b0=0.1; p0=0.5; n0=10
dat0<-c()

#generate income data for m people
for(i in 1:m){dat0[i]<-gmi(y1, b0, p0, n0)}
#dat0

#data<-data.frame(x=dat0)

head(dat0)#人工データの確認
hist(dat0,breaks = 20)

#jagsに渡すincomeデータをインポートする.
# Package the data for shipping to JAGS:
dataList = list(
y = dat0,
N = length(y) 
)

#,n=10#初期値定数n = 10 #初期値定数　ゲーム反復回数
head(dataList$y)
hist(dataList$y,breaks=50)

#------------------------------------------------------------------------------
#jagsにわたすモデルをswich関数で複数定義する．
#modelにあわせて推定するパラメータとアウトプット関数も変更

#modelstring0はパラメータをn,y0,p,b,を全て推定
modelstring0 = "
model{
for( i in 1 : N ) {
y[i] ~ dlnorm(m, 1/(s^2))
}
m <- log(y0)+n*log(1-b)+log((1+b)/(1-b))*n*p
s <- sqrt(n*p*(1-p)*(log((1+b)/(1-b)))^2)
n ~ dpois(lamn)#nをポアソン分布で近似
lamn ~ dunif(5,20)#事前情報として真のn0に近い情報を与える
y0 ~ dpois(lamy0)#
lamy0 ~ dunif(1,1000) #事前情報として真のy0に近い情報を与える
p ~ dbeta(1,1)#p,bの事前分布にβ分布
b ~ dbeta(1,1)}"
# close quote for modelstring

#modelstring1は理論的に固定すべきパラメータn,y0,bを定数にする
#m=2000; y1=10; b0=0.2; p0=0.5; n0=15

modelstring1 = "
model{
for( i in 1 : N ) {
y[i] ~ dlnorm(m, 1/(s^2))
}
y0 <- y1
n <- n0
b <- b0
m <- log(y0)+n*log(1-b)+log((1+b)/(1-b))*n*p
s <- sqrt(n*p*(1-p)*(log((1+b)/(1-b)))^2)
#prior for hyper parameters
p ~ dbeta(1,1)
}"
# close quote for modelstring
# b ~ dbeta(1,1)

#model選択switch
a<-"1"

modelstring<-switch(a,
"0"=modelstring0,
"1"=modelstring1)

#modelstring

writeLines(modelstring,con="model.txt")

#------------------------------------------------------------------------------
# INTIALIZE THE CHAINS.
# Let JAGS do it

#------------------------------------------------------------------------------
# RUN THE CHAINS

require(rjags)

parameters0 = c("p","b","n","lamn","m","s","lamy0","y0")
parameters1 = c("p","m","s")

#parameters1 = c("p","b","m","s")
parameters<-switch(a,
  "0"=parameters0,
  "1"=parameters1)

#parameters = c("p", "b","n","a","m","s","c","y0")
adaptSteps = 1000         # 初期値1000 Number of steps to "tune" the samplers.
burnInSteps = 1000        # 初期値1000 Number of steps to "burn-in" the samplers.
nChains = 3               # 初期値3 Number of chains to run.
numSavedSteps=20000       # 初期値20000 Total number of steps in chains to save.
thinSteps=1               # Number of steps to "thin" (1=keep every step).
nPerChain = ceiling( ( numSavedSteps * thinSteps ) / nChains ) # Steps per chain.
# Create, initialize, and adapt the model:
jagsModel = jags.model( "model.txt" , data=dataList ,
                        n.chains=nChains , n.adapt=adaptSteps )
# Burn-in:
cat( "Burning in the MCMC chain...\n" )
update( jagsModel , n.iter=burnInSteps )
# The saved MCMC chain:
cat( "Sampling final MCMC chain...\n" )
mcmcCoda = coda.samples( jagsModel , variable.names=parameters ,
                            n.iter=nPerChain , thin=thinSteps )

#列はパラメータ毎のchain head(mcmcCoda)

#------------------------------------------------------------------------------
# EXAMINE THE RESULTS

# Display diagnostics of chain, for specified parameters:
parameterNames = varnames(mcmcCoda) # get all parameter names

#個別にdisplayするかどうかを関数show_p(1)で決める
#show_p(0)ならdisplayなし

show_p <-function(x){if (x==1){
for ( parName in parameterNames ) {
  diagMCMC( codaObject=mcmcCoda , parName=parName ,
            saveName=fileNameRoot , saveType=graphFileType )
}}}#diagMCMCを全パラメータに適用するかどうか決める

show_p(1)

# Convert coda-object codaSamples to matrix object for easier handling.
mcmcChain = as.matrix( mcmcCoda )
chainLength = NROW(mcmcChain)

#最後の事後分布出力,openGraphはKruschke作関数
openGraph(width=20,height=12)#(width=10,height=6)
layout(matrix(1:6,nrow=2,byrow=TRUE))#6分割して表示
# posterior predictive
hist( dataList$y , xlab="y" , main="Data w. Post. Pred." , breaks=30 ,
      col="pink" , border="white" , prob=TRUE , cex.lab=1.5)
#dataのヒストグラム

pltIdx = floor(seq(1,chainLength,length=20))#推定パラメータを何個使うか
xComb = seq( min(dataList$y) , max(dataList$y) , length=51 )#dataのレンジを51分割xComb

for ( chnIdx in pltIdx ) {
  lines( xComb ,
         dlnorm( xComb, mcmcChain[chnIdx,"m"], mcmcChain[chnIdx,"s"] ),
         col="skyblue" )
}


# param's of log(y)
postInfo = plotPost( mcmcChain[,"p"] , xlab="success p" )
postInfo = plotPost( mcmcChain[,"b"] , xlab="interest b" )
# param's of y
postInfo = plotPost( mcmcChain[,"n"] , xlab="time" )
#postInfo = plotPost( mcmcChain[,"beta1"] , xlab="beta1")
postInfo = plotPost( mcmcChain[,"lamn"] , xlab="lambda of n")
postInfo = plotPost( mcmcChain[,"y0"] , xlab="y0" )
#postInfo = plotPost( mcmcChain[,"lamy0"] , xlab="lambda of y0")


#postInfo = plotPost( mcmcChain[,"y"] , xlab="mu of y")


#postInfo = plotPost( mcmcChain[,"s"] , xlab="sigma of y")
saveGraph(file=fileNameRoot,type=graphFileType)

#-------------------------------------------------------------------------------

#table(mcmcChain[, "n"])

