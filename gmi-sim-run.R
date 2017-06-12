
#------------------------------------------------------------------------------
# GMIに基づく離散dataがMCMC推定できるか確認
#------------------------------------------------------------------------------
# THE DATA.
# Generate a person's income under the given condition: 
set.seed(47406)

gmi<-function(y0.data, b.data, p.data, n.data){
  #y1:初期資本, b0:投資割合, p0:成功確率, n0:反復回数
  y = y0.data
for(i in 1:n.data){
  if(runif(1) < p.data){
  y <- y + (y*b.data)}
else{y <- y - (y*b.data)} }#for end
y
}

#真のパラメータ設定  y1, b0, p0, n0
n.player=5000 #データ数
y0.data=20; b.data=0.31; p.data=0.2; n.data=50
dat0<-c()

#generate income data for m people
for(i in 1:n.player){dat0[i]<-gmi(y0.data,b.data,p.data,n.data)}
#dat0, m-> n.player

#jagsに渡すincomeデータをインポートする.
# Package the data for shipping to JAGS:
dataList = list(
y = dat0,
N = length(y) 
)#head(dataList$y)

hist(dataList$y,breaks=30)

#THE MODEL.
# modelstrings,parametersの指定はsource("gmi-jags-model***.R")で読み込む
# switch関数の引数numberでモデルを切り替える

number="2" #この数字でjagsモデルを選択

switch(number,
       "0"=source("gmi-sim-model0.R"),#base: \mu,\sigmaを決定論的関数で推定．
       "1"=source("gmi-sim-model1.R"),#GMI　による "p","b","n","m","s"
       "2"=source("gmi-sim-model2.R"),#"p","b","n","m","s"事前分布修正
       "20"=source("gmi-jags-model20.R"),# model20: mを教育年数と年齢に回帰
       "21"=source("gmi-jags-model21.R"),#nj,pi,bi,b0,b1,b2,b3 年齢jでnを階層化．p,b逆ロジット教育年数
       "22"=source("gmi-jags-model22.R"),#ni,pi,bi,b1,b3:p,b逆ロジットで教育年数,n:age
       "3"=source("gmi-jags-model3.R"),#age,eduyが説明変数のあてはめglm
       "31"=source("gmi-jags-model3.R"),#age,eduyが説明変数のあてはめglm
       "40"=source("gmi-jags-model40.R"),#年齢別にp,b,n,giniを推定
       "4"=source("gmi-jags-model4.R"),#年齢別にp[45],b[45],n[45]を推定
       "5"=source("gmi-jags-model.R"))

writeLines(modelstring,con="model.txt")
require(rjags)

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

summary(mcmcCoda)
plot(mcmcCoda)
gelman.diag(mcmcCoda)

para.names = varnames(mcmcCoda) # get all parameter names

coda::gelman.plot(mcmcCoda[,para.names])#shrinkfactor
coda::traceplot(mcmcCoda[,c("p","b","y0")])#traceplot

source("DBDA2E-utilities.R")
# Display diagnostics of chain, for specified parameters:
parameterNames = varnames(mcmcCoda) # get all parameter names
for ( parName in parameterNames ) {
  diagMCMC2( codaObject=mcmcCoda , parName=parName)
}
# diagMCMC2 defined by Hamada does not open windows


graphics.off()



# HDI plot
DbdaDensPlot(mcmcCoda,para.names)
DbdaDensPlot(mcmcCoda,"b")
DbdaDensPlot(mcmcCoda,"y0")

# Autocorrelation
DbdaAcfPlot(mcmcCoda,"y0")


for ( x in para.names ) {
  DbdaDensPlot(mcmcCoda, x)
}

# true values
# y0.data=20; b.data=0.1; p.data=0.6; n.data=50
# EXAMINE THE RESULTS


# Convert coda-object codaSamples to matrix object for easier handling.
mcmcChain = as.matrix( mcmcCoda )
hist(mcmcChain[,"p"])
hist(mcmcChain[,"n"])
hist(mcmcChain[,"b"])
hist(mcmcChain[,"y0"])


coda::traceplot( codaObject[,c(parName)] , main="" , ylab="Param. Value" ,
                 col=DBDAplColors )
