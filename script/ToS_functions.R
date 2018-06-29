## 関数S(x)
Sf <- function(x, delta){
  sinh(delta*asinh(x))
}
## 関数C(x)
Cf <- function(x, delta){
  sqrt(1+Sf(x, delta)^2)
}
## sinh-arcsinh分布
dfas <- function(x, delta){
  Cf(x, delta)*exp(-Sf(x, delta)^2/2)*delta/sqrt(2*pi*(1+x^2))
}
## 関数s(x)
s <- function(x, lambda){
  x + (sqrt(1+lambda^2*x^2)-1)/lambda
}
## 関数s_inverse(x)
s_inverse = function(x, lambda){
  alambda=1-exp(-lambda^2)
  if(lambda !=0){ return(
    (lambda*x+ alambda-alambda*sqrt( (lambda*x+alambda)^2+1-alambda^2))/(lambda*(1-alambda^2)))
  } else if (lambda ==0){ return(x) }
}
## 局度変換を伴うsinh-arcsinh分布
dfas2 <- function(x, mu, sigma, lambda, delta){
  r <-s_inverse( (x-mu)/sigma, lambda)
  return(dfas(r, delta)/sigma)
}

# 局度変換を伴うsinh-arcsinh分布のパラメータ推定関数
mle.dfas2 <- function(x, ini){
  # 対数尤度を計算する関数
  # (最適化の際にパラメータに制約があると望ましくないので指数変換している)
  # (Rの最適化関数が最小化を目指すので正負を入れ替えてある) 
  obj <- function(par){
    mu <- par[1]
    sigma <- exp(par[2])
    lambda <- par[3]
    delta <-par[4]
    lik <-  -sum(log( dfas2(x, mu, sigma, lambda, delta)))
    lik
  }
  # optimで上記の関数が最小になるパラメータiniを探してくれる
  # BFGS法　興味のある人は自分で調べて
  out <- optim( ini, obj, hessian=TRUE)
  out$par2 <- out$par
  out$par2[2] <- exp(out$par[2])
  out
}
# 局度変換を伴うsinh-arcsinh(x)分布の平均，分散，歪度，尖度を計算する関数
fas2.moment<- function(par){
  f1 <- function(x) x*dfas2(x, mu=par[1], sigma=par[2], lambda=par[3], delta=par[4])
  m1 <- integrate( f1, lower=-Inf, upper=Inf)$value
  f2 <- function(x) (x-m1)^2*dfas2(x, mu=par[1], sigma=par[2], lambda=par[3], delta=par[4])
  v2 <- integrate( f2, lower=-Inf, upper=Inf)$value
  f3 <- function(x) (x-m1)^3*dfas2(x, mu=par[1], sigma=par[2], lambda=par[3], delta=par[4])/(v2^{3/2})
  b1 <- integrate( f3, lower=-Inf, upper=Inf)$value
  f4 <- function(x) (x-m1)^4*dfas2(x, mu=par[1], sigma=par[2], lambda=par[3], delta=par[4])/(v2^{2})
  b2 <- integrate( f4, lower=-Inf, upper=Inf)$value
  return(list(m1=m1, v2=v2,b1=b1, b2=b2)) 
}
# 局度変換を伴うsinh-arcsinh(x)分布の分散を計算する関数
fas2.var<- function(par){
  f1 <- function(x) x*dfas2(x, mu=par[1], sigma=par[2], lambda=par[3], delta=par[4])
  m1 <- integrate( f1, lower=-Inf, upper=Inf)$value
  f2 <- function(x) (x-m1)^2*dfas2(x, mu=par[1], sigma=par[2], lambda=par[3], delta=par[4])
  v2 <- integrate( f2, lower=-Inf, upper=Inf)$value
  return(v2=v2) 
}
#局度変換を伴うsinh-arcsinh分布の確率分布関数
pfa2 <- function(x, mu, sigma, lambda, delta){
  f <- function(y) dfas2(y, mu, sigma, lambda, delta)
  return( integrate( f, lower=-Inf, upper=x)$value)
}
#局度変換を伴うsinh-arcsinh分布の 分位点関数
qfas <-function(p, mu, sigma, lambda, delta){
  eps=0.001
  f <- function(x) return( pfa2(x,  mu, sigma, lambda, delta) -p)
  uniroot(f, interval=c(-10,10),extendInt="yes", trace=1)$root
}
#局度変換を伴うsinh-arcsinh分布のESを数値計算で求める関数
find.ES <- function(p, par){
  VaR<- qfas(p, mu=par[1], sigma=par[2], lambda=par[3], delta = par[4])
  #-無限からVaRまでの積分
  f <- function(x) x*dfas2(x,
                           mu=par[1], sigma=par[2], lambda=par[3], delta=par[4])
  #積分した値から期待値を計算
  return( integrate( f, lower=-Inf, upper=VaR)$value/p ) }

## 重点サンンプリングに用いるResampleの関数
Resample1 <- function(data, weight, NofSample){
  re_ind <- runif(NofSample)
  cmwt <- cumsum(weight)/sum(weight);
  st <- sapply(re_ind, function(x) sum(x>cmwt[-length(cmwt)]))
  newdata <- data[ (st+1) ]
  return(newdata)
}


#Fsaからの乱数 SIR並列処理
rfa_SIR_para <- function(n, mu, sigma, lambda, delta)
{
  ## 正規分布を提案分布に
  q <- as.vector(parSapply(cl,1:cl_l,
                           function(x) rnorm(n/cl_l,mean=mu,sd=5*sigma)))
  ## 重み
  w <-parSapply(cl, q, 
                dfas2, mu=mu, sigma=sigma, lambda=lambda, delta=delta)/
    parSapply(cl, q, dnorm, mean=mu, sd=5*sigma)
  ## 合計が1になるように重みを基準化
  w <- w/sum(w)
  ## 重みに従ってresample
  q.resample <- Resample_para(q, weight=w, NofSample = n)
  list(q,q=q.resample, w=w)
}

## 重点サンンプリングに用いるResampleの関数
Resample_para <- function(data, weight, NofSample){
  re_ind <- as.vector(parSapply(cl, 1:cl_l, function(x)runif(NofSample/cl_l) ))
  cmwt <- cumsum(weight)/sum(weight);
  st <- parSapply(cl,re_ind, function(x) sum(x>cmwt[-length(cmwt)]))
  newdata <- data[ (st+1) ]
  return(newdata)
}

#局度変換を伴うsinh-arcsinh分布の単純モンテカルロ法によってESを求める関数
SMC.fa_para <-function(theta){
  out <- c()
  for(i in 1:100){
    ## 正規分布からの重点サンプリングで乱数を取得
    rand.fa<-rfa_SIR_para(n=20000, mu=theta[1], 
                          sigma=theta[2], 
                          lambda=theta[3],
                          delta=theta[4])
    y <- sample(rand.fa$q,10000)
    
    # 単純モンテカルロ法
    #VaRの計算
    VaR1 <- parSapply(cl,100:10000, function(x) quantile(y[1:x], c(0.01,0.025, 0.05)))
    # ESの計算
    Es1 <- parSapply(cl,100:10000, function(x) {
      return(c( mean(y[1:x][y[1:x] < VaR1[1,x-99]]),   
                mean(y[1:x][y[1:x] < VaR1[2,x-99]]),
                mean(y[1:x][y[1:x] < VaR1[3,x-99]])
      ))})
    
    
    
    # 真値と単純モンテカルロ法の結果をまとめる
    out <- cbind(out,t(VaR1),t(Es1))
  }
  return(out)
} 


# 重点サンプリングの関数(d.ISを提案分布に局度変換を伴うsinh-arcsinh分布のESを求める)
IS.fa_single <- function(){
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
  
  N <-10000
  
  #99%点での計算 100~10000までサンプル数を増やして行う
  out1<- sapply(100:N, function(x){
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
  
  out25<-sapply(100:N, function(x){
    out1<-cbind(rfa25[1:x],  w25[1:x]/x)
    A <- out1[sort.list(out1[,1]),]
    A <- cbind(A, cumsum(A[,2]))
    v1 <- A[which.min(abs(A[,3]-0.025)),1]
    es1<- sum(apply(A[1:which.min(abs(A[,3]-0.025)),1:2],1,prod))/0.025
    return(c(v1,es1))})
  
  #matplot(t(out25),type="l")
  #abline(h=VaR.true.FA[2],col=1)
  #abline(h=ES.true.FA[2],col=2)
  
  out5<-sapply(100:N, function(x){
    out1<-cbind( rfa5[1:x],  w5[1:x]/x)
    A <- out1[sort.list(out1[,1]),]
    A <- cbind(A, cumsum(A[,2]))
    v1 <- A[which.min(abs(A[,3]-0.05)),1]
    es1<- sum(apply(A[1:which.min(abs(A[,3]-0.05)),1:2],1,prod))/0.05
    return(c(v1,es1))})
  
  #matplot(t(out5),type="l")
  #abline(h=VaR.true.FA[3],col=1)
  #abline(h=ES.true.FA[3],col=2)
  
  
  return(out = cbind(t(out1),t(out25),t(out5)))
  
}

find.theta <- function(p, par){
  VaR<- qfas(p, mu=par[1], sigma=par[2], lambda=par[3], delta=par[4])
  f <- function(x) return( mean.IS(x,  par) -VaR )
  out <- uniroot(f, interval=c(-10,10),extendInt="yes", trace=1)$root
  out}

# 重点分布の平均を計算する関数
mean.IS <- function(theta,par){
  f <- function(x, theta, par){
    exp(theta*x)*dfas2(x, mu=par[1], sigma=par[2],
                       lambda=par[3], delta=par[4])
  }
  #M <- integrate(f, -Inf, Inf, theta=theta, par=par)$value
  M <- integrate(f, -30, 30, theta=theta, par=par)$value
  m.f<- function(x, theta, par)  x*f(x,theta,par)/M
  # return( integrate(m.f, -Inf, Inf, theta=theta, par=par)$value )
  return( integrate(m.f, -30, 30, theta=theta, par=par)$value )
}

### 指数変換した重点分布の密度関数  
d.IS <- function(x, theta,par){
  # 重点分布の密度関数の分子
  f <- function(x, theta, par){
    exp(theta*x)*dfas2(x, mu=par[1], sigma=par[2],
                       lambda=par[3], delta=par[4])
  }
  # 分母
  M <- integrate(f, -30, 30, theta=theta, par=par)$value
  return( f(x, theta, par)/M )
}


find.theta <- function(p, par){
  VaR<- qfas(p, mu=par[1], sigma=par[2], lambda=par[3], delta=par[4])
  f <- function(x) return( mean.IS(x,  par) -VaR )
  out <- uniroot(f, interval=c(-10,10),extendInt="yes", trace=1)$root
  out}


### 重点分布からサンプリングをする 
### 指数変換した重点分布の密度関数  
d.IS <- function(x, theta,par){
  # 重点分布の密度関数の分子
  f <- function(x, theta, par){
    exp(theta*x)*dfas2(x, mu=par[1], sigma=par[2],
                       lambda=par[3], delta=par[4])
  }
  # 分母
  M <- integrate(f, -30, 30, theta=theta, par=par)$value
  return( f(x, theta, par)/M )
}



# 重点サンプリングの関数(d.ISを提案分布に局度変換を伴うsinh-arcsinh分布のVaR,ESを求める)
IS.fa <- function(){
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
  
  N <-10000
  
  #99%点での計算 100~10000までサンプル数を増やして行う
  out1<- sapply( 100:N, function(x){
    # サンプルの対応するweightをくっつける
    out1<-cbind( rfa1[1:x],  w1[1:x]/x)
    # サンプルの小さい順にならべかえる
    A <- out1[sort.list(out1[,1]),]
    # weightの累積和を並べる
    A <- cbind(A, cumsum(A[,2]))
    # 累積和が0.01に一番近いサンプルが99%VaR
    v1 <- A[which.min(abs(A[,3]-0.01)),1]
    # v1までのサンプルからES0.01の推定値を求める
    es1 <- sum(apply(A[1:which.min(abs(A[,3]-0.01)),1:2],1,prod))/0.01
    return(c(v1, es1))})
  
  #matplot(t(out1),type="l")
  #abline(h=VaR.true.FA[1],col=1)
  #abline(h=ES.true.FA[1],col=2)
  
  out25<-sapply(100:N, function(x){
    out1<-cbind( rfa25[1:x],  w25[1:x]/x)
    A <- out1[sort.list(out1[,1]),]
    A <- cbind(A, cumsum(A[,2]))
    v1<-A[which.min(abs(A[,3]-0.025)),1]
    es1<-   sum(apply(A[1:which.min(abs(A[,3]-0.025)),1:2],1,prod))/0.025
    return(c(v1, es1))})
  
  #matplot(t(out25),type="l")
  #abline(h=VaR.true.FA[2],col=1)
  #abline(h=ES.true.FA[2],col=2)
  
  out5<-sapply(100:N, function(x){
    out1<-cbind( rfa5[1:x],  w5[1:x]/x)
    A <- out1[sort.list(out1[,1]),]
    A <- cbind(A, cumsum(A[,2]))
    v1<-A[which.min(abs(A[,3]-0.05)),1]
    es1<- sum(apply(A[1:which.min(abs(A[,3]-0.05)),1:2],1,prod))/0.05
    return(c(v1, es1))})
  
  #matplot(t(out5),type="l")
  #abline(h=VaR.true.FA[3],col=1)
  #abline(h=ES.true.FA[3],col=2)
  
  
  return( list(out = cbind(t(out1),t(out25),t(out5)), VaR.true.FA, ES.true.FA ) ) }