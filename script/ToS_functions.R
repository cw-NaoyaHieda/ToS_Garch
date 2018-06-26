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
# 局度変換を伴うsinh-arcsinh(x)分布の平均，分散，歪度，尖度を計算する関数
fas2.var<- function(par){
  f1 <- function(x) x*dfas2(x, mu=par[1], sigma=par[2], lambda=par[3], delta=par[4])
  m1 <- integrate( f1, lower=-Inf, upper=Inf)$value
  f2 <- function(x) (x-m1)^2*dfas2(x, mu=par[1], sigma=par[2], lambda=par[3], delta=par[4])
  v2 <- integrate( f2, lower=-Inf, upper=Inf)$value
  return(v2=v2) 
}

#局度変換を伴うsinh-arcsinh分布の 分位点関数
qfas <-function(p, mu, sigma, lambda, delta){
  eps=0.001
  f <- function(x) return( pfa2(x,  mu, sigma, lambda, delta) -p)
  uniroot(f, interval=c(-10,10),extendInt="yes", trace=1)$root
}