#局度変換を伴うsinh-arcsinh分布の単純モンテカルロ法によってVaR,ESを求める関数
SMC.fa_pre <-function(theta){
    ## 正規分布からの重点サンプリングで乱数を取得
  rand.fa<-rfa_SIR(n=20000, mu=theta[1], 
                          sigma=theta[2], 
                          lambda=theta[3],
                          delta=theta[4])
  y <- sample(rand.fa$q,10000)
    
  # 単純モンテカルロ法
  #VaRの計算
  VaR1 <- quantile(y, c(0.01,0.025, 0.05))
  # ESの計算
  Es1 <- c( mean(y[y < VaR1[1]]),   
              mean(y[y < VaR1[2]]),
              mean(y[y < VaR1[3]]))
    
  # 真値と単純モンテカルロ法の結果をまとめる
  out <- cbind(t(VaR1),t(Es1))
  return(out)
} 

# 重点サンプリングの関数(d.ISを提案分布に局度変換を伴うsinh-arcsinh分布のVaR,ESを求める)
IS.fa_pre <- function(){
  f <- function(x, theta, par){
    exp(theta*x)*dfas2(x, mu=par[1], sigma=par[2],
                       lambda=par[3], delta=par[4])
  }
  #weightを計算するためにMを計算する．99%,97.5%,95%をまとめて行う
  M1 <- integrate(f, -30, 30, theta=theta.val1, par=fit$par2)$value
  M25 <- integrate(f, -30, 30, theta=theta.val25, par=fit$par2)$value
  M5 <- integrate(f, -30, 30, theta=theta.val5, par=fit$par2)$value
  
  #weightを計算する
  w1 <- exp(-theta.val1*rfa1)*M1
  w25 <- exp(-theta.val25*rfa25)*M25
  w5 <- exp(-theta.val5*rfa5)*M5
  
  #99%点での計算 100~10000までサンプル数を増やして行う
  out1<-cbind( rfa1,  w1/10000)
  # サンプルの小さい順にならべかえる
  A <- out1[sort.list(out1[,1]),]
  # weightの累積和を並べる
  A <- cbind(A, cumsum(A[,2]))
  # 累積和が0.01に一番近いサンプルが99%VaR
  # v1までのサンプルからES0.01の推定値を求める
  v1 <- A[which.min(abs(A[,3]-0.01)),1]
  es1<- sum(apply(A[1:which.min(abs(A[,3]-0.01)),1:2],1,prod))/0.01
  out1 <- c(v1, es1)
  

  out25<-cbind(rfa25,  w25/10000)
  A <- out25[sort.list(out25[,1]),]
  A <- cbind(A, cumsum(A[,2]))
  v25 <- A[which.min(abs(A[,3]-0.025)),1]
  es25<- sum(apply(A[1:which.min(abs(A[,3]-0.025)),1:2],1,prod))/0.025
  out25 <- c(v25, es25)
  
  out5<-cbind( rfa5,  w5/10000)
  A <- out5[sort.list(out5[,1]),]
  A <- cbind(A, cumsum(A[,2]))
  v5 <- A[which.min(abs(A[,3]-0.05)),1]
  es5<- sum(apply(A[1:which.min(abs(A[,3]-0.05)),1:2],1,prod))/0.05
  out5 <- c(v5, es5)
  
  return(out = cbind(t(out1),t(out25),t(out5)))
  
}

# 正規分布の単純モンテカルロ法
SMC.norm_pre <-function(theta){
  
  f <-　function(x)  x*dnorm(x,mean=theta[1], sd=theta[2])
  VaR.true <- qnorm(c(0.01,0.025,0.05), mean=theta[1],
                    sd=theta[2])
  ES.true  <- sapply( c(0.01,0.025, 0.05), function(x){
    integrate(f,
              lower=-Inf, upper=qnorm(x,mean=theta[1],
                                      sd=theta[2]))$value/x})
  ## Simple Monte Carlo Method VaR　(単純モンテカルロ法)
  y <- rnorm(10000, mean= theta[1], sd =theta[2])
  VaR1 <- quantile(y, c(0.01,0.025, 0.05))
  ## Simple Monte Carlo Method ES　(単純モンテカルロ法)
  Es1 <- c( mean(y[y<VaR1[1]]),   
            mean(y[y<VaR1[2]]),
            mean(y[y<VaR1[3]]))
  
  out <- cbind(t(VaR1),t(Es1))
  
  return(out)
} 

# 正規分布のIS
IS.norm_pre<-function(theta)
{
  VaR.true <- qnorm(c(0.01,0.025,0.05), mean=theta[1], sd=theta[2])
  f <- function(x)  x*dnorm(x,mean=theta[1], sd=theta[2]) 
  
  ES.true <-sapply(c(0.01,0.025, 0.05), function(x){
    integrate(f, lower=-Inf, upper=qnorm(x,mean=theta[1], sd=theta[2])
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

# 重点分布(d.IS)からのサンプリング (fucntionsでの同じ関数は並列にしているので上書き)
rIS_SIR <- function(n, par, par2, theta){
  ## 正規分布を提案分布に
  q <- rnorm(n,mean=par2,sd=15)
  ## 重み
  w <- sapply(q, 
                 d.IS, theta =theta, par=par)/dnorm(q, mean=par2, sd=15)
  w <- w/sum(w)
  ## resample
  q.resample <- Resample1(q, weight=w, NofSample = n)
  list( q=q.resample, w=w)
}
# 重点分布(dfas2)からのサンプリング (fucntionsでの同じ関数は並列にしているので上書き)
rIS_SIR2 <- function(n, par, par2){
  ## 正規分布を提案分布に
  q <- rnorm(n, mean=par2[1], sd=par2[2])
  ## 重み
  w <- sapply(q, 
              dfas2, mu = par[1], sigma= par[2] , lambda=par[3], delta = par[4])/dnorm(q, mean=par2[1], sd=par2[2])
  w <- w/sum(w)
  ## resample
  q.resample <- Resample1(q, weight=w, NofSample = n)
  list( q=q.resample, w=w)
}





#局度変換を伴うsinh-arcsinh分布の単純モンテカルロ法によってVaR,ESを求める関数並列
SMC.fa_para_pre <-function(theta){
  ## 正規分布からの重点サンプリングで乱数を取得
  rand.fa<-rfa_SIR_para(n=20000, mu=theta[1], 
                   sigma=theta[2], 
                   lambda=theta[3],
                   delta=theta[4])
  y <- sample(rand.fa$q,10000)
  
  # 単純モンテカルロ法
  #VaRの計算
  VaR1 <- quantile(y, c(0.01,0.025, 0.05))
  # ESの計算
  Es1 <- c( mean(y[y < VaR1[1]]),   
            mean(y[y < VaR1[2]]),
            mean(y[y < VaR1[3]]))
  
  # 真値と単純モンテカルロ法の結果をまとめる
  out <- cbind(t(VaR1),t(Es1))
  return(out)
} 

# 重点分布(d.IS)からのサンプリング (fucntionsでの同じ関数は並列にしているので上書き)
rIS_SIR_para <- function(n, par, par2, theta){
  ## 正規分布を提案分布に
  q <- as.vector(parSapply(cl,1:cl_l,function(x) rnorm(n/cl_l,mean=par2,sd=15)))
  ## 重み
  w <- parSapply(cl,q, 
                 d.IS, theta =theta, par=par)/dnorm(q, mean=par2, sd=15)
  w <- w/sum(w)
  ## resample
  q.resample <- Resample_para(q, weight=w, NofSample = n)
  list( q=q.resample, w=w)
}

# 重点サンプリングの関数(d.ISを提案分布に局度変換を伴うsinh-arcsinh分布のVaR,ESを求める)
IS.fa_pre <- function(){
  f <- function(x, theta, par){
    exp(theta*x)*dfas2(x, mu=par[1], sigma=par[2],
                       lambda=par[3], delta=par[4])
  }
  #weightを計算するためにMを計算する．99%,97.5%,95%をまとめて行う
  M1 <- integrate(f, -30, 30, theta=theta.val1, par=fit$par2)$value
  M25 <- integrate(f, -30, 30, theta=theta.val25, par=fit$par2)$value
  M5 <- integrate(f, -30, 30, theta=theta.val5, par=fit$par2)$value
  
  #weightを計算する
  w1 <- exp(-theta.val1*rfa1)*M1
  w25 <- exp(-theta.val25*rfa25)*M25
  w5 <- exp(-theta.val5*rfa5)*M5
  
  #99%点での計算 100~10000までサンプル数を増やして行う
  out1<-cbind( rfa1,  w1/10000)
  # サンプルの小さい順にならべかえる
  A <- out1[sort.list(out1[,1]),]
  # weightの累積和を並べる
  A <- cbind(A, cumsum(A[,2]))
  # 累積和が0.01に一番近いサンプルが99%VaR
  # v1までのサンプルからES0.01の推定値を求める
  v1 <- A[which.min(abs(A[,3]-0.01)),1]
  es1<- sum(apply(A[1:which.min(abs(A[,3]-0.01)),1:2],1,prod))/0.01
  out1 <- c(v1, es1)
  
  
  out25<-cbind(rfa25,  w25/10000)
  A <- out25[sort.list(out25[,1]),]
  A <- cbind(A, cumsum(A[,2]))
  v25 <- A[which.min(abs(A[,3]-0.025)),1]
  es25<- sum(apply(A[1:which.min(abs(A[,3]-0.025)),1:2],1,prod))/0.025
  out25 <- c(v25, es25)
  
  out5<-cbind( rfa5,  w5/10000)
  A <- out5[sort.list(out5[,1]),]
  A <- cbind(A, cumsum(A[,2]))
  v5 <- A[which.min(abs(A[,3]-0.05)),1]
  es5<- sum(apply(A[1:which.min(abs(A[,3]-0.05)),1:2],1,prod))/0.05
  out5 <- c(v5, es5)
  
  return(out = cbind(t(out1),t(out25),t(out5)))
  
}

# 重点サンプリングの関数(d.ISを提案分布に局度変換を伴うsinh-arcsinh分布のVaR,ESを求める)
IS.fa_para_pre <- function(){
  f <- function(x, theta, par){
    exp(theta*x)*dfas2(x, mu=par[1], sigma=par[2],
                       lambda=par[3], delta=par[4])
  }
  #weightを計算するためにMを計算する．99%,97.5%,95%をまとめて行う
  M1 <- integrate(f, -30, 30, theta=theta.val1, par=fit$par2)$value
  M25 <- integrate(f, -30, 30, theta=theta.val25, par=fit$par2)$value
  M5 <- integrate(f, -30, 30, theta=theta.val5, par=fit$par2)$value
  
  #weightを計算する
  w1 <- exp(-theta.val1*rfa1)*M1
  w25 <- exp(-theta.val25*rfa25)*M25
  w5 <- exp(-theta.val5*rfa5)*M5
  
  #99%点での計算 100~10000までサンプル数を増やして行う
  out1<-cbind( rfa1,  w1/10000)
  # サンプルの小さい順にならべかえる
  A <- out1[sort.list(out1[,1]),]
  # weightの累積和を並べる
  A <- cbind(A, cumsum(A[,2]))
  # 累積和が0.01に一番近いサンプルが99%VaR
  # v1までのサンプルからES0.01の推定値を求める
  v1 <- A[which.min(abs(A[,3]-0.01)),1]
  es1<- sum(parApply(cl,A[1:which.min(abs(A[,3]-0.01)),1:2], 1, prod))/0.01
  out1 <- c(v1, es1)
  
  
  out25<-cbind(rfa25,  w25/10000)
  A <- out25[sort.list(out25[,1]),]
  A <- cbind(A, cumsum(A[,2]))
  v25 <- A[which.min(abs(A[,3]-0.025)),1]
  es25<- sum(parApply(cl, A[1:which.min(abs(A[,3]-0.025)),1:2], 1, prod))/0.025
  out25 <- c(v25, es25)
  
  out5<-cbind( rfa5,  w5/10000)
  A <- out5[sort.list(out5[,1]),]
  A <- cbind(A, cumsum(A[,2]))
  v5 <- A[which.min(abs(A[,3]-0.05)),1]
  es5<- sum(parApply(cl, A[1:which.min(abs(A[,3]-0.05)),1:2], 1, prod))/0.05
  out5 <- c(v5, es5)
  
  return(out = cbind(t(out1),t(out25),t(out5)))
  
}


# 局度変換を伴うsinh-arcsinh分布のパラメータ推定関数
mle.dfas2_weight <- function(x, weight, ini){
  # 対数尤度を計算する関数
  # (最適化の際にパラメータに制約があると望ましくないので指数変換している)
  # (Rの最適化関数が最小化を目指すので正負を入れ替えてある) 
  obj <- function(par){
    mu <- par[1]
    sigma <- exp(par[2])
    lambda <- par[3]
    delta <-par[4]
    lik <-  -sum(weight * log( dfas2(x, mu, sigma, lambda, delta)))
    lik
  }
  # optimで上記の関数が最小になるパラメータiniを探してくれる
  # BFGS法　興味のある人は自分で調べて
  out <- optim( ini, obj, hessian=TRUE)
  out$par2 <- out$par
  out$par2[2] <- exp(out$par[2])
  out
}

IS.fa_loacation_pre <- function(){
  par <- fit$par2
  #weightを計算する
  w1 <- dfas2(rfa1, mu=par[1], sigma=par[2],lambda=par[3], delta=par[4])/
    dfas2(rfa1, mu=ES1.fa, sigma=par[2],lambda=par[3], delta=par[4])
  w25 <- dfas2(rfa25, mu=par[1], sigma=par[2],lambda=par[3], delta=par[4])/
    dfas2(rfa25, mu=ES25.fa, sigma=par[2],lambda=par[3], delta=par[4])
  w5 <- dfas2(rfa5, mu=par[1], sigma=par[2],lambda=par[3], delta=par[4])/
    dfas2(rfa5, mu=ES5.fa, sigma=par[2],lambda=par[3], delta=par[4])
  
  #99%点での計算 100~10000までサンプル数を増やして行う
  out1<-cbind( rfa1,  w1/10000)
  # サンプルの小さい順にならべかえる
  A <- out1[sort.list(out1[,1]),]
  # weightの累積和を並べる
  A <- cbind(A, cumsum(A[,2]))
  # 累積和が0.01に一番近いサンプルが99%VaR
  # v1までのサンプルからES0.01の推定値を求める
  v1 <- A[which.min(abs(A[,3]-0.01)),1]
  es1<- sum(apply(A[1:which.min(abs(A[,3]-0.01)),1:2],1,prod))/0.01
  out1 <- c(v1, es1)
  
  
  out25<-cbind(rfa25,  w25/10000)
  A <- out25[sort.list(out25[,1]),]
  A <- cbind(A, cumsum(A[,2]))
  v25 <- A[which.min(abs(A[,3]-0.025)),1]
  es25<- sum(apply(A[1:which.min(abs(A[,3]-0.025)),1:2],1,prod))/0.025
  out25 <- c(v25, es25)
  
  out5<-cbind( rfa5,  w5/10000)
  A <- out5[sort.list(out5[,1]),]
  A <- cbind(A, cumsum(A[,2]))
  v5 <- A[which.min(abs(A[,3]-0.05)),1]
  es5<- sum(apply(A[1:which.min(abs(A[,3]-0.05)),1:2],1,prod))/0.05
  out5 <- c(v5, es5)
  
  return(out = cbind(t(out1),t(out25),t(out5)))
  
}



