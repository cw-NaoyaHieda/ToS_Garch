
# 局度変換を伴うsinh-arcsinh分布のパラメータ推定関数
mle.dt <- function(x, ini){
  # 対数尤度を計算する関数
  # (最適化の際にパラメータに制約があると望ましくないので指数変換している)
  # (Rの最適化関数が最小化を目指すので正負を入れ替えてある) 
  obj <- function(par){
    df <- exp(par[1])
    ncp <- exp(par[2])
    lik <-  -sum(log(dt(x=x,df=df,ncp=ncp)))
    lik
  }
  # optimで上記の関数が最小になるパラメータiniを探してくれる
  # BFGS法　興味のある人は自分で調べて
  out <- optim( ini, obj, hessian=TRUE)
  out$par2[1] <- exp(out$par[1])
  out$par2[2] <- exp(out$par[2])
  out
}


tmp<- mle.dt(df$log_x[-1],ini=c(0,0))

ggplot(df, aes(x=log_x))+
  geom_histogram(binwidth=.3, colour="black", fill="white",
                 aes(y=..density..))+
  #scale_fill_gradient("Count", low="#DCDCDC", high="#7C7C7C")+
  xlim(min(rt),max(rt))+
  stat_function(fun=dnorm,
                #colour="violetred",
                aes(colour ="Normal distribution"),
                size=1,
                n=401,
                args=list(mean=mean(rt), sd = sd(rt)))+
  stat_function(fun=dfas2,
                #colour="green",
                aes(colour ="Tos sinh-arcsinh distribution"),
                size=1,
                n=401,
                args=list(mu=fit$par2[1],
                          sigma = fit$par2[2],
                          lambda = fit$par2[3],
                          delta = fit$par2[4]))+
  stat_function(fun=dt,
                #colour="green",
                aes(colour ="T distribution"),
                size=1,
                n=401,
                args=list(df=tmp$par2[1],
                          ncp=tmp$par2[2]))+
  theme_bw()+
  theme(legend.position = "bottom")+
  theme(axis.title.x = element_text(size=25),axis.title.y = element_text(size=25))+
  theme(axis.text.x = element_text(size=25),axis.text.y = element_text(size=25)) +
  theme(legend.title = element_text(size=25),legend.text = element_text(size=25))



t_var <- qt(c(0.01,0.025,0.05),df=tmp$par2[1],ncp=tmp$par2[2])
f <-　function(x)  x*dt(x,df=theta[1], ncp=theta[2])
## Simple Monte Carlo Method VaR　(単純モンテカルロ法)
y_ <- rt(10000, df= tmp$par2[1], ncp =tmp$par2[2])

t_SMC_1000 <- pforeach::pforeach(n = 1:100, .combine = rbind,.cores = )({
  #for(n in 100:10000){
  y <- rt(1000, df= tmp$par2[1], ncp =tmp$par2[2])
  out1<- sapply(1000, function(x){
    VaR1 <- quantile(y[1:x], c(0.01,0.025))
    ## Simple Monte Carlo Method ES　(単純モンテカルロ法)
    Es1 <- c( mean(y[1:x][y[1:x]<VaR1[1]]),   
              mean(y[1:x][y[1:x]<VaR1[2]]))
    return(c(VaR1,Es1))
  })
  t(out1)
})


VaR.true <- qt(c(0.01,0.025,0.05), df=tmp$par2[1], ncp=tmp$par2[2])
f <- function(x)  x*dt(x,df=tmp$par2[1], ncp=tmp$par2[2]) 

ES.true <-sapply(c(0.01,0.025, 0.05), function(x){
  integrate(f, lower=-Inf, upper=qt(x,df=tmp$par2[1], ncp=tmp$par2[2])
  )$value/x} )

t_IS_1000 <- pforeach::pforeach(n = 1:100, .combine = rbind,.cores = 45)({
  
  #正規分布を提案分布とした重点サンプリング
  rfa1 <-rnorm(10000, mean = ES.true[1], sd =1)
  rfa25 <-rnorm(10000, mean = ES.true[2], sd = 1)
  
  w1 <-dt(rfa1,df=tmp$par2[1], ncp=tmp$par2[2]) / dnorm(rfa1,mean = ES.true[1], sd =1)
  w25 <-dt(rfa25,df=tmp$par2[1], ncp=tmp$par2[2]) / dnorm(rfa25,mean = ES.true[2], sd =1)
  
  N <-10000
  
  #99%点での計算 100~10000までサンプル数を増やして行う
  out1<- sapply(1000, function(x){
    # サンプルの対応するweightをくっつける
    out1<-cbind( rfa1[1:x],  w1[1:x]/x)
    # サンプルの小さい順にならべかえる
    A <- out1[sort.list(out1[,1]),]
    # weightの累積和を並べる
    A <- cbind(A, cumsum(A[,2]))
    # 累積和が0.01に一番近いサンプルが99%VaR
    # v1までのサンプルからES0.01の推定値を求める
    v1 <- A[which.min(abs(A[,3]-0.01)),1]
    es1<- sum(apply(A[1:which.min(abs(A[,3]-0.01)),1:2],1,prod))/0.01
    return(c(v1,es1))})
  
  #matplot(t(out1),type="l")
  #abline(h=VaR.true.FA[1],col=1)
  #abline(h=ES.true.FA[1],col=2)
  
  out25<-sapply(1000, function(x){
    out1<-cbind(rfa25[1:x],  w25[1:x]/x)
    A <- out1[sort.list(out1[,1]),]
    A <- cbind(A, cumsum(A[,2]))
    v1 <- A[which.min(abs(A[,3]-0.025)),1]
    es1<- sum(apply(A[1:which.min(abs(A[,3]-0.025)),1:2],1,prod))/0.025
    return(c(v1,es1))})
  
  
  
  c(out1,out25)
  
  #matplot(t(out5),type="l")
  #abline(h=VaR.true[3],col=1)
  #abline(h=ES.true[3],col=2)
})


colMeans(t_SMC_1000)
apply(t_SMC_1000,2,sd)
colMeans(t_IS_1000)
apply(t_IS_1000,2,sd)


t_SMC_1000_e <- colMeans(t_SMC_1000)
t_IS_1000_e <- colMeans(t_IS_1000)

var_matrix <- sapply(df$log_x[-1], function(x) x<c(t_SMC_1000_e[c(1,2)],t_IS_1000_e[c(1,3)])) %>% t()
colSums(var_matrix)/dim(var_matrix)[1]

Ec_99 <- length(df$log_x) * 0.01
Ec_97.5 <- length(df$log_x) * 0.025

sum((df$log_x[-1][which(var_matrix[,1] == TRUE)] / t_SMC_1000_e[3]  -1)/Ec_99)
sum((df$log_x[-1][which(var_matrix[,2] == TRUE)] / t_SMC_1000_e[4]  -1)/Ec_97.5)
sum((df$log_x[-1][which(var_matrix[,3] == TRUE)] / t_IS_1000_e[2]  -1)/Ec_99)
sum((df$log_x[-1][which(var_matrix[,4] == TRUE)] / t_IS_1000_e[4]  -1)/Ec_97.5)



SMC.t_pre <-function(theta){
  
  f <-　function(x)  x*dt(x,df=theta[1], ncp=theta[2])
  VaR.true <- qt(c(0.01,0.025,0.05), df=theta[1],
                 ncp=theta[2])
  ES.true  <- sapply( c(0.01,0.025, 0.05), function(x){
    integrate(f,
              lower=-Inf, upper=qt(x,df=theta[1],
                                   ncp=theta[2]))$value/x})
  ## Simple Monte Carlo Method VaR　(単純モンテカルロ法)
  y <- rt(10000, df= theta[1], ncp =theta[2])
  VaR1 <- quantile(y, c(0.01,0.025, 0.05))
  ## Simple Monte Carlo Method ES　(単純モンテカルロ法)
  Es1 <- c( mean(y[y<VaR1[1]]),   
            mean(y[y<VaR1[2]]),
            mean(y[y<VaR1[3]]))
  
  out <- cbind(t(VaR1),t(Es1))
  
  return(out)
} 

# 正規分布のIS
IS.t_pre<-function(theta)
{
  VaR.true <- qt(c(0.01,0.025,0.05), df=theta[1], ncp=theta[2])
  f <- function(x)  x*dt(x,df=theta[1], ncp=theta[2]) 
  
  ES.true <-sapply(c(0.01,0.025, 0.05), function(x){
    integrate(f, lower=-Inf, upper=qt(x,df=theta[1], ncp=theta[2])
    )$value/x} )
  
  th1 <-Expected.val.th(VaR=VaR.true[1], mu=theta[1], sd = theta[2])
  th25 <-Expected.val.th(VaR=VaR.true[2], mu=theta[1], sd = theta[2])
  th5 <-Expected.val.th(VaR=VaR.true[3], mu=theta[1], sd = theta[2])
  
  
  #--------------------------------
  IS1 <-rnorm(10000, mean = theta[1]+theta[2]^2*th1, sd = theta[2])
  IS25 <-rnorm(10000, mean = theta[1]+theta[2]^2*th25, sd = theta[2])
  IS5 <-rnorm(10000, mean = theta[1]+theta[2]^2*th5, sd = theta[2])
  
  M1 <- exp(theta[1]*th1+ theta[2]^2*th1^2/2)  
  w1 <- exp(-th1*IS1)*M1
  
  M25 <- exp(theta[1]*th25+ theta[2]^2*th25^2/2)  
  w25 <- exp(-th25*IS25)*M25
  
  M5 <- exp(theta[1]*th5+ theta[2]^2*th5^2/2)  
  w5 <- exp(-th5*IS5)*M5
  
  out1<-cbind( IS1,  w1/10000)
  A <- out1[sort.list(out1[,1]),]
  A <- cbind(A, cumsum(A[,2]))
  # VaR0.01
  v1 <- A[which.min(abs(A[,3]-0.01)),1]
  # ES0.01の推定値
  es1<-   sum(apply(A[1:which.min(abs(A[,3]-0.01)),1:2],1,prod))/0.01 
  out1_ <- c(v1, es1)
  
  out1<-cbind( IS25,  w25/10000)
  A <- out1[sort.list(out1[,1]),]
  A <- cbind(A, cumsum(A[,2]))
  # VaR0.01
  v1<-A[which.min(abs(A[,3]-0.025)),1]
  # ES0.01の推定値
  es1<-   sum(apply(A[1:which.min(abs(A[,3]-0.025)),1:2],1,prod))/0.025 
  out25 <- c(v1, es1)
  
  out1<-cbind( IS5,  w5/10000)
  A <- out1[sort.list(out1[,1]),]
  A <- cbind(A, cumsum(A[,2]))
  # VaR0.01
  v1 <- A[which.min(abs(A[,3]-0.05)),1]
  # ES0.01の推定値
  es1 <- sum(apply(A[1:which.min(abs(A[,3]-0.05)),1:2],1,prod))/0.05 
  out5 <- c(v1, es1)
  
  #matplot(t(out5),type="l")
  #abline(h=VaR.true[3],col=1)
  #abline(h=ES.true[3],col=2)
  return(out = cbind(t(out1_),t(out25),t(out5)))  }




result_tmp <- c()
for(i in 1:dim(result_para)[1]){
  rt <- df$log_x[(i+1):(i+250)]
  tmp<- mle.dt(df$log_return[-1],ini=c(0,0))
  
  SMC.norm.out <- SMC.t_pre(theta)
  IS.norm.out <- IS.t_pre(theta)
  result_tmp <- rbind(result_tmp,data.frame(SMC.norm.out,IS.norm.out))
  print(i)
}

