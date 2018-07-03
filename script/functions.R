# GARCHを適用する関数
garch_f_opt <- function(data, order = c(1, 1), error = "norm", dis_para = c()){
  t <- length(data)
  sigma_t <- rep(0, t)
  garch_loglik <- rep(0, t)
  #GARCHのbeta
  G_beta <- order[1]
  #GARCHのalpha
  G_alpha <- order[2]
  params <- rep(0,c(G_beta + G_alpha + 1))
  # garchモデルの対数尤度を計算する関数(関数によって変えてますが、多分もっとスマートに書けますね・・・)
  garch_op <- function(params){
    #sigma_tの計算
    #計算できないところまでは、全体のボラティリティを与える
    for(i in 1:G_beta){
      sigma_t[i] <- sum(data^2)/t
      garch_loglik[i] <- log(dnorm(data[i], mean = 0, sd = sqrt(sigma_t[i])))
    }
    for(i in (G_beta+1):t){
      sigma_t[i] <- exp(params[1]) +
        exp(params[2:(G_beta+1)]) * sigma_t[(i-G_beta):(i-1)] +
        exp(params[(G_beta+2):(G_beta+G_alpha+1)]) * data[(i-G_alpha):(i-1)]^2
      garch_loglik[i] <- log(dnorm(data[i], mean = 0, sd = sqrt(sigma_t[i])))
    }
    return(-sum(garch_loglik))
  }
  
  # garchモデルの対数尤度を計算する関数(FA分布版)  
  garch_op_ToS <- function(params){
    #sigma_tの計算
    #計算できないところまでは、仮に推定してるsigmaを与える
    for(i in 1:G_beta){
      #sigma_t[i] <- sum(data^2)/t
      # ToS のsigma!=標準偏差なので
      sigma_t[i] <- tos_sigma^2
      garch_loglik[i] <- log(dfas2(data[i], delta=tos_delta,mu=tos_mu,sigma=sigma_t[i] %>% sqrt(),lambda=tos_lambda))
    }
    for(i in (G_beta+1):t){
      sigma_t[i] <- exp(params[1]) +
        exp(params[2:(G_beta+1)]) * sigma_t[(i-G_beta):(i-1)] +
        exp(params[(G_beta+2):(G_beta+G_alpha+1)]) * data[(i-G_alpha):(i-1)]^2
      garch_loglik[i] <- log(dfas2(data[i], delta=tos_delta,mu=tos_mu,sigma=sigma_t[i] %>% sqrt(),lambda=tos_lambda))
    }
    return(-sum(garch_loglik))
  }
  
  # garchモデルの対数尤度を計算する関数(FA分布版)  
  garch_op_ToS_stand <- function(params){
    #sigma_tの計算
    #計算できないところまでは、仮に推定してるsigmaを与える
    for(i in 1:G_beta){
      #sigma_t[i] <- sum(data^2)/t
      # ToS のsigma!=標準偏差なので
      sigma_t[i] <- sum(data^2)/t
      garch_loglik[i] <- log(dfas2_stand_sd(data[i], delta=tos_delta,mu=tos_mu,sigma=tos_sigma,lambda=tos_lambda,sd=sigma_t[i]%>% sqrt()))
    }
    for(i in (G_beta+1):t){
      sigma_t[i] <- exp(params[1]) +
        exp(params[2:(G_beta+1)]) * sigma_t[(i-G_beta):(i-1)] +
        exp(params[(G_beta+2):(G_beta+G_alpha+1)]) * data[(i-G_alpha):(i-1)]^2
      garch_loglik[i] <- log(dfas2_stand_sd(data[i], delta=tos_delta,mu=tos_mu,sigma=tos_sigma,lambda=tos_lambda,sd=sigma_t[i]%>% sqrt()))
    }
    return(-sum(garch_loglik))
  }
  
  # garchモデルの対数尤度を計算する関数(FA分布版)  
  garch_op_ToS_stand_with_roll <- function(params){
    #sigma_tの計算
    #計算できないところまでは、仮に推定してるsigmaを与える
    for(i in 1:G_beta){
      #sigma_t[i] <- sum(data^2)/t
      # ToS のsigma!=標準偏差なので
      sigma_t[i] <- sum(data^2)/t
      garch_loglik[i] <- log(dfas2_stand_sd(data[i], delta=tos_delta[i],mu=tos_mu[i],
                                            sigma=tos_sigma[i],lambda=tos_lambda[i],sd=sigma_t[i]%>% sqrt()))
    }
    for(i in (G_beta+1):t){
      sigma_t[i] <- exp(params[1]) +
        exp(params[2:(G_beta+1)]) * sigma_t[(i-G_beta):(i-1)] +
        exp(params[(G_beta+2):(G_beta+G_alpha+1)]) * data[(i-G_alpha):(i-1)]^2
      garch_loglik[i] <- log(dfas2_stand_sd(data[i], delta=tos_delta[i],mu=tos_mu[i],
                                            sigma=tos_sigma[i],lambda=tos_lambda[i],sd=sigma_t[i]%>% sqrt()))
    }
    return(-sum(garch_loglik))
  }
  
  # garchモデルの対数尤度を計算する関数(FA分布版)  
  garch_op_ToS_with_roll <- function(params){
    #sigma_tの計算
    #計算できないところまでは、仮に推定してるsigmaを与える
    for(i in 1:G_beta){
      #sigma_t[i] <- sum(data^2)/t
      # ToS のsigma!=標準偏差なので
      sigma_t[i] <- sum(data^2)/t
      garch_loglik[i] <- log(dfas2_stand_sd(data[i], delta=tos_delta[i],mu=tos_mu[i],
                                   sigma=sigma_t[i] ,lambda=tos_lambda[i],sd=sigma_t[i]%>% sqrt()))
    }
    for(i in (G_beta+1):t){
      sigma_t[i] <- exp(params[1]) +
        exp(params[2:(G_beta+1)]) * sigma_t[(i-G_beta):(i-1)] +
        exp(params[(G_beta+2):(G_beta+G_alpha+1)]) * data[(i-G_alpha):(i-1)]^2
      garch_loglik[i] <- log(dfas2(data[i], delta=tos_delta[i],mu=tos_mu[i],
                                   sigma=sigma_t[i] ,lambda=tos_lambda[i]))
    }
    return(-sum(garch_loglik))
  }
  
  if(grepl(error,"norm")){out <- optim(params, garch_op)}else{
    if(grepl(error,"ToS") && is.null(dim(dis_para))){
      tos_mu <- dis_para["mu"]
      tos_sigma <- dis_para["sigma"]
      tos_delta <- dis_para["delta"]
      tos_lambda <- dis_para["lambda"]
      out <- optim(params, garch_op_ToS)}else{
        if(grepl(error,"ToS_stand")&& is.null(dim(dis_para))){
          tos_mu <- dis_para["mu"]
          tos_sigma <- dis_para["sigma"]
          tos_delta <- dis_para["delta"]
          tos_lambda <- dis_para["lambda"]
          out <- optim(params, garch_op_ToS_stand)
      }else{
        if(grepl(error,"ToS_stand")){
          tos_mu <- pull(dis_para,mu)
          tos_sigma <- pull(dis_para,sigma)
          tos_delta <- pull(dis_para,delta)
          tos_lambda <- pull(dis_para,lambda)
          out <- optim(params, garch_op_ToS_stand_with_roll)
        }else{
        if(grepl(error,"ToS")){
          tos_mu <- pull(dis_para,mu)
          tos_sigma <- pull(dis_para,sigma)
          tos_delta <- pull(dis_para,delta)
          tos_lambda <- pull(dis_para,lambda)
          out <- optim(params, garch_op_ToS_with_roll)
        }
        }
      }
      }
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
    if(grepl(error,"ToS") && is.null(dim(dis_para))){
      for(i in 1:p){
        sigma_t[i] <- dis_para["sigma"]^2
      }
    }else{
      if(grepl(error,"ToS_stand")){
        for(i in 1:p){
          sigma_t[i] <- sum(data^2)/t
        }
      }else{
      if(grepl(error,"ToS")){
        for(i in 1:p){
          sigma_t[i] <- pull(dis_para,sigma)[i]^2
        }
    }
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



# 流石にややこしいので、ToS分布のパラメータまで推定する分布は別版
garch_f_opt_with_para <- function(data, order = c(1, 1)){
  t <- length(data)
  sigma_t <- rep(0, t)
  garch_loglik <- rep(0, t)
  #GARCHのbeta
  G_beta <- order[1]
  #GARCHのalpha
  G_alpha <- order[2]
  params <- rep(0,c(G_beta + G_alpha + 1 + 4))
  params_len <- length(params)
  params[(params_len-3):params_len] <- c(0, log(0.2), -0.2, 0.5)
  
  # garchモデルの対数尤度を計算する関数(FA分布版)  
  garch_op_ToS_stand <- function(params){
    tos_mu <- params[params_len - 3]
    tos_sigma <- exp(params[params_len - 2])
    tos_lambda <- params[params_len - 1]
    tos_delta <- params[params_len]
    #sigma_tの計算
    #計算できないところまでは、仮に推定してるsigmaを与える
    for(i in 1:G_beta){
      #sigma_t[i] <- sum(data^2)/t
      # ToS のsigma!=標準偏差なので
      sigma_t[i] <- sum(data^2)/t
      garch_loglik[i] <- log(dfas2_stand_sd(data[i], delta=tos_delta,mu=tos_mu,sigma=tos_sigma,lambda=tos_lambda,sd=sigma_t[i]%>% sqrt()))
    }
    for(i in (G_beta+1):t){
      sigma_t[i] <- exp(params[1]) +
        exp(params[2:(G_beta+1)]) * sigma_t[(i-G_beta):(i-1)] +
        exp(params[(G_beta+2):(G_beta+G_alpha+1)]) * data[(i-G_alpha):(i-1)]^2
      garch_loglik[i] <- log(dfas2_stand_sd(data[i], delta=tos_delta,mu=tos_mu,sigma=tos_sigma,lambda=tos_lambda,sd=sigma_t[i]%>% sqrt()))
    }
    return(-sum(garch_loglik))
  }
 
  out <- optim(params, garch_op_ToS_stand)
     
  out
}


# 流石にややこしいので、ToS分布のパラメータまで推定する分布は別版
garch_f_opt_with_para2 <- function(data, order = c(1, 1), tos_mu, tos_sigma){
  t <- length(data)
  sigma_t <- rep(0, t)
  garch_loglik <- rep(0, t)
  #GARCHのbeta
  G_beta <- order[1]
  #GARCHのalpha
  G_alpha <- order[2]
  params <- rep(0,c(G_beta + G_alpha + 1 + 2))
  params_len <- length(params)
  params[c(params_len-1,params_len)] <- c(-0.2, 0.5)
  
  # garchモデルの対数尤度を計算する関数(FA分布版)  
  garch_op_ToS_stand <- function(params){
    tos_lambda <- params[params_len - 1]
    tos_delta <- params[params_len]
    #sigma_tの計算
    #計算できないところまでは、仮に推定してるsigmaを与える
    for(i in 1:G_beta){
      #sigma_t[i] <- sum(data^2)/t
      # ToS のsigma!=標準偏差なので
      sigma_t[i] <- sum(data^2)/t
      garch_loglik[i] <- log(dfas2_stand_sd(data[i], delta=tos_delta,mu=tos_mu,sigma=tos_sigma,lambda=tos_lambda,sd=sigma_t[i]%>% sqrt()))
    }
    for(i in (G_beta+1):t){
      sigma_t[i] <- exp(params[1]) +
        exp(params[2:(G_beta+1)]) * sigma_t[(i-G_beta):(i-1)] +
        exp(params[(G_beta+2):(G_beta+G_alpha+1)]) * data[(i-G_alpha):(i-1)]^2
      garch_loglik[i] <- log(dfas2_stand_sd(data[i], delta=tos_delta,mu=tos_mu,sigma=tos_sigma,lambda=tos_lambda,sd=sigma_t[i]%>% sqrt()))
    }
    return(-sum(garch_loglik))
  }
  
  out <- optim(params, garch_op_ToS_stand)
  
  out
}

# lambdaは0に固定してやる
garch_f_opt_with_para2_lambda0 <- function(data, order = c(1, 1), tos_mu, tos_sigma){
  t <- length(data)
  sigma_t <- rep(0, t)
  garch_loglik <- rep(0, t)
  #GARCHのbeta
  G_beta <- order[1]
  #GARCHのalpha
  G_alpha <- order[2]
  params <- rep(0,c(G_beta + G_alpha + 1 + 1))
  params_len <- length(params)
  params[c(params_len)] <- c(0.5)
  tos_lambda <- 0
  # garchモデルの対数尤度を計算する関数(FA分布版)  
  garch_op_ToS_stand <- function(params){
    tos_delta <- params[params_len]
    #sigma_tの計算
    #計算できないところまでは、仮に推定してるsigmaを与える
    for(i in 1:G_beta){
      #sigma_t[i] <- sum(data^2)/t
      # ToS のsigma!=標準偏差なので
      sigma_t[i] <- sum(data^2)/t
      garch_loglik[i] <- log(dfas2_stand_sd(data[i], delta=tos_delta,mu=tos_mu,sigma=tos_sigma,lambda=tos_lambda,sd=sigma_t[i]%>% sqrt()))
    }
    for(i in (G_beta+1):t){
      sigma_t[i] <- exp(params[1]) +
        exp(params[2:(G_beta+1)]) * sigma_t[(i-G_beta):(i-1)] +
        exp(params[(G_beta+2):(G_beta+G_alpha+1)]) * data[(i-G_alpha):(i-1)]^2
      garch_loglik[i] <- log(dfas2_stand_sd(data[i], delta=tos_delta,mu=tos_mu,sigma=tos_sigma,lambda=tos_lambda,sd=sigma_t[i]%>% sqrt()))
    }
    return(-sum(garch_loglik))
  }
  
  out <- optim(params, garch_op_ToS_stand)
  
  out
}

# mu+sigma*errorで行う
garch_f_opt_with_para2_muestimate <- function(data, order = c(1, 1), tos_mu, tos_sigma){
  t <- length(data)
  sigma_t <- rep(0, t)
  garch_loglik <- rep(0, t)
  #GARCHのbeta
  G_beta <- order[1]
  #GARCHのalpha
  G_alpha <- order[2]
  params <- rep(0,c(G_beta + G_alpha + 1 +1 + 2))
  params_len <- length(params)
  params[c(params_len-1,params_len)] <- c(-0.2, 0.5)
  params[c(params_len-2)] <- 0.05
  # garchモデルの対数尤度を計算する関数(FA分布版)  
  garch_op_ToS_stand <- function(params){
    tos_lambda <- params[params_len - 1]
    tos_delta <- params[params_len]
    #sigma_tの計算
    #計算できないところまでは、仮に推定してるsigmaを与える
    for(i in 1:G_beta){
      #sigma_t[i] <- sum(data^2)/t
      # ToS のsigma!=標準偏差なので
      sigma_t[i] <- sum(data^2)/t
      garch_loglik[i] <- log(dfas2_stand_sd((data[i] - params[params_len - 2]),
                                            delta=tos_delta,mu=tos_mu,sigma=tos_sigma,
                                            lambda=tos_lambda,sd=sigma_t[i]%>% sqrt()))
    }
    for(i in (G_beta+1):t){
      sigma_t[i] <- exp(params[1]) +
        exp(params[2:(G_beta+1)]) * sigma_t[(i-G_beta):(i-1)] +
        exp(params[(G_beta+2):(G_beta+G_alpha+1)]) * ((data[(i-G_alpha):(i-1)] - params[params_len - 2])^2)
      garch_loglik[i] <- log(dfas2_stand_sd( (data[i] - params[params_len - 2]), delta=tos_delta,mu=tos_mu,sigma=tos_sigma,
                                            lambda=tos_lambda,sd=sigma_t[i]%>% sqrt()))
    }
    return(-sum(garch_loglik))
  }
  
  out <- optim(params, garch_op_ToS_stand)
  
  out
}


# 正規分布で平均を考慮する
garch_f_opt_with_para_norm <- function(data, order = c(1, 1), tos_mu, tos_sigma){
  t <- length(data)
  sigma_t <- rep(0, t)
  garch_loglik <- rep(0, t)
  #GARCHのbeta
  G_beta <- order[1]
  #GARCHのalpha
  G_alpha <- order[2]
  params <- rep(0,c(G_beta + G_alpha + 1 + 1))
  params_len <- length(params)
  # garchモデルの対数尤度を計算する関数(FA分布版)  
  garch_op <- function(params){
    #sigma_tの計算
    #計算できないところまでは、全体のボラティリティを与える
    for(i in 1:G_beta){
      sigma_t[i] <- sum(data^2)/t
      garch_loglik[i] <- log(dnorm(data[i], mean = params[params_len], sd = sqrt(sigma_t[i])))
    }
    for(i in (G_beta+1):t){
      sigma_t[i] <- exp(params[1]) +
        exp(params[2:(G_beta+1)]) * sigma_t[(i-G_beta):(i-1)] +
        exp(params[(G_beta+2):(G_beta+G_alpha+1)]) * (data[(i-G_alpha):(i-1)]-params[params_len])^2
      garch_loglik[i] <- log(dnorm(data[i], mean = params[params_len], sd = sqrt(sigma_t[i])))
    }
    return(-sum(garch_loglik))
  }
  
  out <- optim(params, garch_op)
  
  out
}

# t分布で平均を考慮する
garch_f_opt_with_para_t <- function(data, order = c(1, 1), tos_mu, tos_sigma){
  t <- length(data)
  sigma_t <- rep(0, t)
  garch_loglik <- rep(0, t)
  #GARCHのbeta
  G_beta <- order[1]
  #GARCHのalpha
  G_alpha <- order[2]
  params <- rep(0,c(G_beta + G_alpha + 1 + 1))
  params_len <- length(params)
  # garchモデルの対数尤度を計算する関数(FA分布版)  
  garch_op <- function(params){
    #sigma_tの計算
    #計算できないところまでは、全体のボラティリティを与える
    for(i in 1:G_beta){
      sigma_t[i] <- sum(data^2)/t
      garch_loglik[i] <- log(dt(data[i], df = params[params_len])))
    }
    for(i in (G_beta+1):t){
      sigma_t[i] <- exp(params[1]) +
        exp(params[2:(G_beta+1)]) * sigma_t[(i-G_beta):(i-1)] +
        exp(params[(G_beta+2):(G_beta+G_alpha+1)]) * data[(i-G_alpha):(i-1)]^2
      garch_loglik[i] <- log(dnorm(data[i], mean = params[params_len], sd = sqrt(sigma_t[i])))
    }
    return(-sum(garch_loglik))
  }
  
  out <- optim(params, garch_op)
  
  out
}

garch_f_mu <- function(data, params,mu, order = c(1, 1), error = "norm", dis_para){
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
    if(grepl(error,"ToS") && is.null(dim(dis_para))){
      for(i in 1:p){
        sigma_t[i] <- dis_para["sigma"]^2
      }
    }else{
      if(grepl(error,"ToS_stand")){
        for(i in 1:p){
          sigma_t[i] <- sum(data^2)/t
        }
      }else{
        if(grepl(error,"ToS")){
          for(i in 1:p){
            sigma_t[i] <- pull(dis_para,sigma)[i]^2
          }
        }
      }
    }
  }
  for(i in (p+1):t){
    sigma_t[i] <- exp(params[1]) +
      exp(params[2:(p+1)]) * sigma_t[(i-p):(i-1)] +
      exp(params[(p+2):(p+q+1)]) * (data[(i-q):(i-1)]- mu)^2
  }
  
  return(sigma_t)
}
