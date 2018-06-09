#ToS_Garch直下で実行 scriptフォルダじゃない

#garch関数
library(tseries)
#sp <- read.csv("data/nky.csv",header=TRUE,skip=4)[,2] 
sp <- read.csv("data/SP500.csv",header=TRUE)[,2]
sp %>% as.character() %>% as.numeric()
sp.ret <- sp/lag(sp)-1
sp.ret <- sp.ret[-1]
sp.ts.garch <- garch(sp.ret) #デフォルトはGARCH(1,1)

#sigma: fitted values
sp.ts.sigma <- xts(sp.ts.garch$fitted.values[,1]*sqrt(250), order.by=as.Date(read.csv("data/nky.csv",header=TRUE,skip=4)[-1,1]))

#plot index and volatility
default.par <- par()
mai <- par()$mai
mai[4] <- mai[1]
par(mai = mai)
plot(sp, type="l", ylab="index", main="SP500", minor.ticks=F, auto.grid=F)
par(new=T)
plot(sp.ts.sigma, type="l", axes=F, ylab="", col=2, main="", minor.ticks=F, auto.grid=F)
axis(4)
mtext("volatility", side=4, line=3)
legend("topleft", c("index(LHS)", "volatility(RHS)"), lty=rep(1,2), col=c(1,2), bty="n")
par(default.par)

#----------------------------------------------------------------------------
# MLE - GARCH(1,1)
#----------------------------------------------------------------------------
sp.n <- length(sp.ret)
sp.garch.logLikelihood <- rep(0, sp.n)
sp.garch.sigma <- rep(0, sp.n) #実際にはsigma^2が格納される
sp.garch.l <- function(params){
  #1: const
  #2: beta(sigma)
  #3: alpha(epsilon)
  for(i in 1:sp.n){
    if(i == 1){
      sp.garch.sigma[i] <<- sum(sp.ret^2)/sp.n
    }else{
      sp.garch.sigma[i] <<- exp(params[1]) + exp(params[2]) * sp.garch.sigma[i-1] + exp(params[3]) * sp.ret[i-1]^2
    }
    sp.garch.logLikelihood[i] <<- log(dnorm(sp.ret[i], mean = 0, sd = sqrt(sp.garch.sigma[i])))
  }
  return(sum(sp.garch.logLikelihood))
}

sp.garch.mle <- optim(c(0, 0, 0), sp.garch.l, control=list(fnscale=-1))

#parameter comparison
sp.ts.garch$coef #tseries
exp(sp.garch.mle$par) #using likelihood function

sp.garch.sigma <- xts(sqrt(sp.garch.sigma*250), order.by=as.Date(read.csv("data/nky.csv",header=TRUE,skip=4)[-1,1]))

plot(sp.garch.sigma,  type="l", main="MLE GARCH volatility ", ylab="", minor.ticks=F)
lines(sp.ts.sigma[,1], col=2)
plot(sp.ts.sigma[,1],  type="l", main="MLE GARCH volatility ", ylab="", minor.ticks=F)
legend("topright", c("optim", "tseries"), lty=rep(1,2), col=c(1,2), bty="n")

plot(sp.garch.sigma-sp.ts.sigma[,1], type="l", main="Vol diff", minor.ticks=F)