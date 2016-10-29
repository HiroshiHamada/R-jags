
#------------------------------------------------------------------------------

# Jags-Ymet-Xnom1grp-MlogNormal-Script.R
# April 19, 2016. John K. Kruschke.
# Requires auxiliary R scripts that accompany the book,
# Kruschke, J. K. (2015). Doing Bayesian Data Analysis, Second Edition:
# A Tutorial with R, JAGS, and Stan. Academic Press / Elsevier.

# ��������gmc��MCMC����ł��Ă���̂��ǂ�����l�H�f�[�^�Ń`�F�b�N

# graphics.off()
# rm(list=ls(all=TRUE))

source("DBDA2E-utilities.R")
fileNameRoot="jags-gmc-sim"
graphFileType="png"

#------------------------------------------------------------------------------
# THE DATA.
# Generate a person's income under the given condition: 
set.seed(47406)

#data�̕��U������Ȃ��̂ŁCp0���m���ϐ�������D

gmi<-function(y1, b0, p0, n0){
  #y1:�������{, b0:��������, p0:�����m��, n0:������
  y = y1
for(i in 1:n0){
  if(runif(1) < p0 + rnorm(1,0,0.1)){
  y <- y + y*(b0+rnorm(1,0,0.01))}
else{y <- y - y*(b0+rnorm(1,0,0.1))} }#for end
y
}


#�^�̃p�����[�^�ݒ�  y1, b0, p0, n0
m=2000; y1=20; b0=0.1; p0=0.5; n0=10
dat0<-c()

#generate income data for m people
for(i in 1:m){dat0[i]<-gmi(y1, b0, p0, n0)}
#dat0

#data<-data.frame(x=dat0)

head(dat0)#�l�H�f�[�^�̊m�F
hist(dat0,breaks = 20)

#jags�ɓn��income�f�[�^���C���|�[�g����.
# Package the data for shipping to JAGS:
dataList = list(
y = dat0,
N = length(y) 
)

#,n=10#�����l�萔n = 10 #�����l�萔�@�Q�[��������
head(dataList$y)
hist(dataList$y,breaks=50)

#------------------------------------------------------------------------------
#jags�ɂ킽�����f����swich�֐��ŕ�����`����D
#model�ɂ��킹�Đ��肷��p�����[�^�ƃA�E�g�v�b�g�֐����ύX

#modelstring0�̓p�����[�^��n,y0,p,b,��S�Đ���
modelstring0 = "
model{
for( i in 1 : N ) {
y[i] ~ dlnorm(m, 1/(s^2))
}
m <- log(y0)+n*log(1-b)+log((1+b)/(1-b))*n*p
s <- sqrt(n*p*(1-p)*(log((1+b)/(1-b)))^2)
n ~ dpois(lamn)#n���|�A�\�����z�ŋߎ�
lamn ~ dunif(5,20)#���O���Ƃ��Đ^��n0�ɋ߂�����^����
y0 ~ dpois(lamy0)#
lamy0 ~ dunif(1,1000) #���O���Ƃ��Đ^��y0�ɋ߂�����^����
p ~ dbeta(1,1)#p,b�̎��O���z�Ƀ����z
b ~ dbeta(1,1)}"
# close quote for modelstring

#modelstring1�͗��_�I�ɌŒ肷�ׂ��p�����[�^n,y0,b��萔�ɂ���
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

#model�I��switch
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
adaptSteps = 1000         # �����l1000 Number of steps to "tune" the samplers.
burnInSteps = 1000        # �����l1000 Number of steps to "burn-in" the samplers.
nChains = 3               # �����l3 Number of chains to run.
numSavedSteps=20000       # �����l20000 Total number of steps in chains to save.
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

#��̓p�����[�^����chain head(mcmcCoda)

#------------------------------------------------------------------------------
# EXAMINE THE RESULTS

# Display diagnostics of chain, for specified parameters:
parameterNames = varnames(mcmcCoda) # get all parameter names

#�ʂ�display���邩�ǂ������֐�show_p(1)�Ō��߂�
#show_p(0)�Ȃ�display�Ȃ�

show_p <-function(x){if (x==1){
for ( parName in parameterNames ) {
  diagMCMC( codaObject=mcmcCoda , parName=parName ,
            saveName=fileNameRoot , saveType=graphFileType )
}}}#diagMCMC��S�p�����[�^�ɓK�p���邩�ǂ������߂�

show_p(1)

# Convert coda-object codaSamples to matrix object for easier handling.
mcmcChain = as.matrix( mcmcCoda )
chainLength = NROW(mcmcChain)

#�Ō�̎��㕪�z�o��,openGraph��Kruschke��֐�
openGraph(width=20,height=12)#(width=10,height=6)
layout(matrix(1:6,nrow=2,byrow=TRUE))#6�������ĕ\��
# posterior predictive
hist( dataList$y , xlab="y" , main="Data w. Post. Pred." , breaks=30 ,
      col="pink" , border="white" , prob=TRUE , cex.lab=1.5)
#data�̃q�X�g�O����

pltIdx = floor(seq(1,chainLength,length=20))#����p�����[�^�����g����
xComb = seq( min(dataList$y) , max(dataList$y) , length=51 )#data�̃����W��51����xComb

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
