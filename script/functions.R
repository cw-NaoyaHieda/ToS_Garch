garch_f_opt <- function(data, order = c(1, 1), error = "norm", dis_para = c()){
  t <- length(data)
  sigma_t <- rep(0, t)
  garch_loglik <- rep(0, t)
  p <- order[1]
  q <- order[2]
  params <- rep(0,c(p + q + 1))
  # garchモデルの対数尤度を計算する関数
  garch_op <- function(params){
    #sigma_tの計算
    #計算できないところまでは、全体のボラティリティを与える
    for(i in 1:p){
      sigma_t[i] <- sum(data^2)/t
      garch_loglik[i] <- log(dnorm(data[i], mean = 0, sd = sqrt(sigma_t[i])))
    }
    for(i in (p+1):t){
      sigma_t[i] <- exp(params[1]) +
        exp(params[2:(p+1)]) * sigma_t[(i-p):(i-1)] +
        exp(params[(p+2):(p+q+1)]) * data[(i-q):(i-1)]^2
      garch_loglik[i] <- log(dnorm(data[i], mean = 0, sd = sqrt(sigma_t[i])))
    }
    return(-sum(garch_loglik))
    }
    
  garch_op_ToS <- function(params){
    #sigma_tの計算
    #計算できないところまでは、全体のボラティリティを与える
    for(i in 1:p){
      #sigma_t[i] <- sum(data^2)/t
      # ToS のsigma!=標準偏差なので
      sigma_t[i] <- tos_sigma^2
      garch_loglik[i] <- log(dfas2(data[i], delta=tos_delta,mu=tos_mu,sigma=sigma_t[i] %>% sqrt(),lambda=tos_lambda))
    }
    for(i in (p+1):t){
      sigma_t[i] <- exp(params[1]) +
        exp(params[2:(p+1)]) * sigma_t[(i-p):(i-1)] +
        exp(params[(p+2):(p+q+1)]) * data[(i-q):(i-1)]^2
      garch_loglik[i] <- log(dfas2(data[i], delta=tos_delta,mu=tos_mu,sigma=sigma_t[i] %>% sqrt(),lambda=tos_lambda))
    }
    return(-sum(garch_loglik))
    }
  
  if(grepl(error,"norm")){out <- optim(params, garch_op)}else{
    if(grepl(error,"ToS")){
      tos_mu <- dis_para["mu"]
      tos_sigma <- dis_para["sigma"]
      tos_delta <- dis_para["delta"]
      tos_lambda <- dis_para["lambda"]
      out <- optim(params, garch_op_ToS)}
  }
  out
}


# garchモデルのsigmaを計算する関数
garch_f <- function(data, params, order = c(1, 1), error = "norm", dis_para){
  t <- length(data)
  sigma_t <- rep(0, t)
  if(length(params) != order[1] + order[2] + 1 ){
    stop("The number of parameters is wrong")
  }
  p <- order[1]
  q <- order[2]
  #sigma_tの計算
  #計算できないところまでは、全体のボラティリティを与える
  if(grepl(error,"norm")){
    for(i in 1:p){
      sigma_t[i] <- sum(data^2)/t
    }
  }else{
    if(grepl(error,"ToS")){
      for(i in 1:p){
        sigma_t[i] <- dis_para["sigma"]^2
      }
    }
  }
  for(i in (p+1):t){
    sigma_t[i] <- exp(params[1]) +
      exp(params[2:(p+1)]) * sigma_t[(i-p):(i-1)] +
      exp(params[(p+2):(p+q+1)]) * data[(i-q):(i-1)]^2
  }
  
  return(sigma_t)
}