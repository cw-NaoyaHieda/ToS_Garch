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