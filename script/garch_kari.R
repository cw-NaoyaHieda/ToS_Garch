source("script/functions.R")
load(file = "data/GARCH/with_ToSpara2.Rdata")
my_fun_res2_with_para2 <- garch_f(df$log_x[-1], my_fun_res_with_para2$par[1:3],
                                  error = "ToS_stand")

load(file = "data/GARCH/with_ToSpara2_muestimate.Rdata")
#分散の計算
my_fun_res2_with_para2 <- garch_f(df$log_x[-1], my_fun_res_with_para2_mues$par[1:3],
                                  error = "ToS_stand")+ my_fun_res_with_para2_mues$par[4]

my_fun_res <- garch_f_opt(df$log_x[-1])
my_fun_res_mu <- garch_f_opt_with_para_norm(df$log_x[-1])
my_fun_res2 <- garch_f(df$log_x[-1], my_fun_res$par)
my_fun_res2_with_para2 <- garch_f(df$log_x[-1], my_fun_res_with_para2$par[1:3],
                                  error = "ToS_stand")
my_fun_res2_mu <- garch_f(df$log_x[-1], my_fun_res_mu$par[1:3]) + my_fun_res_mu$par[4]
my_fun_res2_with_para2_mu <- garch_f(df$log_x[-1], my_fun_res_with_para2_mues$par[1:3],
                                     error = "ToS_stand") + my_fun_res_with_para2_mues$par[4]

plot_d  <- data.frame(dt = df$dt[-1],
                      my_fun_res2 = my_fun_res2 %>% sqrt(),
                      my_fun_tos2_with_para2 = my_fun_res2_with_para2 %>% sqrt(),
                      my_fun_res2_mu = my_fun_res2_mu %>% sqrt(),
                      my_fun_tos2_with_para2_mu = my_fun_res2_with_para2_mu %>% sqrt())


# 局度変換を伴うsinh-arcsinh(x)分布の平均，分散，歪度，尖度を計算する関数
fas2.mu.sd<- function(par){
  f1 <- function(x) x*dfas2(x, mu=par[1], sigma=par[2], lambda=par[3], delta=par[4])
  m1 <- integrate( f1, lower=-Inf, upper=Inf)$value
  f2 <- function(x) (x-m1)^2*dfas2(x, mu=par[1], sigma=par[2], lambda=par[3], delta=par[4])
  v2 <- integrate( f2, lower=-Inf, upper=Inf)$value
  return(c(m1,sqrt(v2))) 
}

dfas2_stand <- function(x, mu, sigma, lambda, delta){
  mu.sd <- fas2.mu.sd(par = c(mu,sigma,lambda,delta))
  dfas2(mu.sd[1] + mu.sd[2]*x, mu, sigma, lambda, delta)*mu.sd[2]
}


##多分、元の分布の標準偏差をいじれる分布が必要
dfas2_stand_sd <- function(x, mu, sigma, lambda, delta, sd){
  mu.sd <- fas2.mu.sd(par = c(mu,sigma,lambda,delta))
  dfas2((mu.sd[1] + mu.sd[2]*x/sd), mu, sigma, lambda, delta)*mu.sd[2]/sd
}

garch_SMC_fsa_var <- pforeach::pforeach(i = 1:dim(plot_d)[1], .combine = rbind,.cores = 45)({
  #Fsastandからの乱数 SIR並列処理
  rfa_stand_SIR<- function(n, mu, sigma, lambda, delta, sd)
  {
    ## 正規分布を提案分布に
    q <- rnorm(n,mean=mu ,sd=5*sd)
    ## 重み
    w <- sapply(q, dfas2_stand_sd, mu=mu, sigma=sigma, lambda=lambda, delta=delta,sd = sd)/
      dnorm(q, mean=mu , sd=5*sd)
    ## 合計が1になるように重みを基準化
    w <- w/sum(w)
    ## 重みに従ってresample
    q.resample <- Resample1(q, weight=w, NofSample = n)
    list(q,q=q.resample, w=w)
  }
  ## 正規分布からの重点サンプリングで乱数を取得
  rand.fa<-rfa_stand_SIR(n=20000, mu=0, 
                     sigma=1, 
                     lambda=my_fun_res_with_para2$par[4],
                     delta=my_fun_res_with_para2$par[5],
                     sd = plot_d[i,3])
  y <- sample(rand.fa$q,10000)
    
    # 単純モンテカルロ法
    #VaRの計算
  VaR1 <- quantile(y, c(0.01,0.025, 0.05))
    # ESの計算
  Es1 <- c( mean(y[y < VaR1[1]]),   
              mean(y[y < VaR1[2]]),
              mean(y[y < VaR1[3]]))
    
    # 真値と単純モンテカルロ法の結果をまとめる
  cbind(t(VaR1),t(Es1))
  
})


garch_SMC_fsa_mu_var <- pforeach::pforeach(i = 1:dim(plot_d)[1], .combine = rbind,.cores = 45)({
  #Fsastandからの乱数 SIR並列処理
  rfa_stand_SIR<- function(n, mu, sigma, lambda, delta, sd)
  {
    ## 正規分布を提案分布に
    q <- rnorm(n,mean=mu ,sd=5*sd)
    ## 重み
    w <- sapply(q, dfas2_stand_sd, mu=mu, sigma=sigma, lambda=lambda, delta=delta,sd = sd)/
      dnorm(q, mean=mu , sd=5*sd)
    ## 合計が1になるように重みを基準化
    w <- w/sum(w)
    ## 重みに従ってresample
    q.resample <- Resample1(q, weight=w, NofSample = n)
    list(q,q=q.resample, w=w)
  }
  ## 正規分布からの重点サンプリングで乱数を取得
  rand.fa<-rfa_stand_SIR(n=20000, mu=0, 
                         sigma=1, 
                         lambda=my_fun_res_with_para2_mues$par[5],
                         delta=my_fun_res_with_para2_mues$par[6],
                         sd = plot_d[i,5])
  rand.fa$q <- rand.fa$q + my_fun_res_with_para2_mues$par[4]
  y <- sample(rand.fa$q,10000)
  
  # 単純モンテカルロ法
  #VaRの計算
  VaR1 <- quantile(y, c(0.01,0.025, 0.05))
  # ESの計算
  Es1 <- c( mean(y[y < VaR1[1]]),   
            mean(y[y < VaR1[2]]),
            mean(y[y < VaR1[3]]))
  
  # 真値と単純モンテカルロ法の結果をまとめる
  cbind(t(VaR1),t(Es1))
  
})

garch_SMC_norm_var <- pforeach::pforeach(i = 1:dim(plot_d)[1], .combine = rbind,.cores = 45)({
  
  ## 正規分布からの乱数を取得
  y<-rnorm(n=10000, mean=0,
                         sd=plot_d[i,2])
  
  # 単純モンテカルロ法
  #VaRの計算
  VaR1 <- quantile(y, c(0.01,0.025, 0.05))
  # ESの計算
  Es1 <- c( mean(y[y < VaR1[1]]),   
            mean(y[y < VaR1[2]]),
            mean(y[y < VaR1[3]]))
  
  # 真値と単純モンテカルロ法の結果をまとめる
  cbind(t(VaR1),t(Es1))
  
})

garch_SMC_norm_mu_var <- pforeach::pforeach(i = 1:dim(plot_d)[1], .combine = rbind,.cores = 45)({
  
  ## 正規分布からの乱数を取得
  y<-rnorm(n=10000, mean =my_fun_res_mu$par[4],
                 sd=plot_d[i,4])
  
  # 単純モンテカルロ法
  #VaRの計算
  VaR1 <- quantile(y, c(0.01,0.025, 0.05))
  # ESの計算
  Es1 <- c( mean(y[y < VaR1[1]]),   
            mean(y[y < VaR1[2]]),
            mean(y[y < VaR1[3]]))
  
  # 真値と単純モンテカルロ法の結果をまとめる
  cbind(t(VaR1),t(Es1))
  
})


garch_SMC <- cbind(garch_SMC_norm_var,garch_SMC_norm_mu_var,garch_SMC_fsa_var,garch_SMC_fsa_mu_var)

#save(list=c("garch_SMC"),file="data/garch_SMC.Rdata")
load(file="data/garch_SMC.Rdata")
garch_var99 <- cbind(df[-1,-2],garch_SMC[,c(1,7,13,19)])

check_VaR <- apply(garch_var99[,-1], 1,
                   function(x) x[1] < x[-1]) %>% t()
count_VaR <- colSums(check_VaR)
tmp <- count_VaR/ dim(check_VaR)[1]


eval_ES <- check_VaR * apply(garch_var99[,c(-1,-2)],2,
                             function(x) abs(df$log_x[-1] - x)/ -x)
tmp3 <- colMeans(eval_ES)



garch_var975 <- cbind(df[-1,-2],garch_SMC[,c(2,8,14,20)])

check_VaR <- apply(garch_var975[,-1], 1,
                   function(x) x[1] < x[-1]) %>% t()
count_VaR <- colSums(check_VaR)
tmp2 <- count_VaR/ dim(check_VaR)[1]

eval_ES <- check_VaR * apply(garch_var975[,c(-1,-2)],2,
                             function(x) abs(df$log_x[-1] - x)/ -x)
tmp4 <- colMeans(eval_ES)

tmp_matrix <- cbind(tmp,tmp2,tmp3,tmp4)

xtable(tmp_matrix,dig=4)

p1 <- ggplot() + 
  geom_line(data = garch_var99[,c(-2,-3,-5)] %>% gather(key=fun,value,-dt),
            aes(x = dt, y= value,color=fun)) +
  ylab("value") +
  theme_bw() +
  theme(axis.title.x = element_text(size=25),axis.title.y = element_text(size=25))+
  theme(axis.text.x = element_text(size=25),axis.text.y = element_text(size=25)) +
  theme(legend.title = element_text(size=25),legend.text = element_text(size=25)) +
  theme(legend.position = 'none')

p2 <- ggplot() + 
  geom_line(data = garch_var975[,c(-2,-3,-5)] %>% gather(key=fun,value,-dt),
            aes(x = dt, y= value,color=fun)) +
  ylab("value") +
  theme_bw() +
  theme(axis.title.x = element_text(size=25),axis.title.y = element_text(size=25))+
  theme(axis.text.x = element_text(size=25),axis.text.y = element_text(size=25)) +
  theme(legend.title = element_text(size=25),legend.text = element_text(size=25)) +
  theme(legend.position = 'bottom') +
  scale_color_hue(name = "", labels = c(expression(Norm),expression(ToS)))
grid.arrange(p1, p2,
             ncol = 1) 



sig<-function(x){(tanh(x)+1)/2}
sig_env<-function(y){(1/2)*log(y/(1-y))}

#シグモイド変換で
sig<-function(x){(tanh(x)+1)/2}
sig_env<-function(y){(1/2)*log(y/(1-y))}
#シグモイド変換で
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
        sig(params[2:(G_beta+1)]) * sigma_t[(i-G_beta):(i-1)] +
        sig(params[(G_beta+2):(G_beta+G_alpha+1)]) * data[(i-G_alpha):(i-1)]^2
      garch_loglik[i] <- log(dnorm(data[i], mean = params[params_len], sd = sqrt(sigma_t[i])))
    }
    return(-sum(garch_loglik))
  }
  
  out <- optim(params, garch_op)
  
  out
}
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
        sig(params[2:(G_beta+1)]) * sigma_t[(i-G_beta):(i-1)] +
        sig(params[(G_beta+2):(G_beta+G_alpha+1)]) * ((data[(i-G_alpha):(i-1)] - params[params_len - 2])^2)
      garch_loglik[i] <- log(dfas2_stand_sd( (data[i] - params[params_len - 2]), delta=tos_delta,mu=tos_mu,sigma=tos_sigma,
                                             lambda=tos_lambda,sd=sigma_t[i]%>% sqrt()))
    }
    return(-sum(garch_loglik))
  }
  
  out <- optim(params, garch_op_ToS_stand)
  
  out
}

# garchモデルのsigmaを計算する関数
garch_f_mu <- function(data, params, mu,order = c(1, 1), error = "norm", dis_para){
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
      sig(params[2:(p+1)]) * sigma_t[(i-p):(i-1)] +
      sig(params[(p+2):(p+q+1)]) * (data[(i-q):(i-1)]-mu)^2
  }
  
  return(sigma_t)
}


garch_norm_rolling <- pforeach::pforeach(i = 251:dim(df)[1], .combine = rbind,.cores = 45)({
  rt <-df$log_x[(i-249):i]
  my_fun_res_mu <- garch_f_opt_with_para_norm(rt)
  my_fun_res2_mu <- garch_f_mu(rt,params =  my_fun_res_mu$par[1:3],mu=my_fun_res_mu$par[4]) %>% sqrt() + my_fun_res_mu$par[4]
  c(my_fun_res_mu$par,my_fun_res2_mu[length(my_fun_res2_mu)])
})

garch_fsa_rolling <- pforeach::pforeach(i = 251:dim(df)[1], .combine = rbind,.cores = 45)({
  rt <-df$log_x[(i-249):i]
  my_fun_res_with_para2_mues <- garch_f_opt_with_para2_muestimate(
  rt,tos_mu = 0,tos_sigma = 1)
  my_fun_res2_with_para2_mu <- garch_f_mu(rt, my_fun_res_with_para2_mues$par[1:3],
                                       error = "ToS_stand",
                                       mu=my_fun_res_with_para2_mues$par[4]) %>% sqrt() + my_fun_res_with_para2_mues$par[4]
  c(my_fun_res_with_para2_mues$par,my_fun_res2_with_para2_mu[length(my_fun_res2_with_para2_mu)])
})

garch_rolling <- cbind(garch_norm_rolling,garch_fsa_rolling)

#save(list=c("garch_rolling"),file="data/garch_rolling.Rdata")
load(file="data/garch_rolling.Rdata")
garch_parameters <- data.frame(exp(garch_rolling[,c(1)]),
                               sig(garch_rolling[,c(2,3)]),
                               exp(garch_rolling[,c(6)]),
                               sig(garch_rolling[,c(7,8)]))
colnames(garch_parameters) <- c("norma0","normb1","norma1","tosa0","tosb1","tosa1")
ggplot() + geom_path(data=data.frame(dt=df$dt[251:length(df$dt)],garch_parameters) %>%
                       gather(key,value,-dt),
                     aes(x=dt,y=value,color=key))
sigma_garhc_rolling <- garch_rolling[,c(5,12)]
sigma_garhc_rolling <- sigma_garhc_rolling %>% data.frame()
colnames(sigma_garhc_rolling) <- c("N","ToS")
ggplot() + geom_path(data=data.frame(dt=df$dt[251:dim(df)[1]],sigma_garhc_rolling)%>% gather(key,value,-dt),
                     aes(x=dt,y=value,color=key))+theme_bw()


garch_rolling %>% head()
sigma_garhc_rolling %>% data.frame() %>% gather(key,value)


garch_SMC_fsa_mu_var <- pforeach::pforeach(i = 1:dim(garch_rolling)[1], .combine = rbind,.cores = 45)({
  #Fsastandからの乱数 SIR並列処理
  rfa_stand_SIR<- function(n, mu, sigma, lambda, delta, sd)
  {
    ## 正規分布を提案分布に
    q <- rnorm(n,mean=mu ,sd=5*sd)
    ## 重み
    w <- sapply(q, dfas2_stand_sd, mu=mu, sigma=sigma, lambda=lambda, delta=delta,sd = sd)/
      dnorm(q, mean=mu , sd=5*sd)
    ## 合計が1になるように重みを基準化
    w <- w/sum(w)
    ## 重みに従ってresample
    q.resample <- Resample1(q, weight=w, NofSample = n)
    list(q,q=q.resample, w=w)
  }
  ## 正規分布からの重点サンプリングで乱数を取得
  rand.fa<-rfa_stand_SIR(n=20000, mu=0, 
                         sigma=1, 
                         lambda=garch_rolling[i,10] %>% as.numeric(),
                         delta=garch_rolling[i,11]%>% as.numeric(),
                         sd = garch_rolling[i,12]%>% as.numeric())
  rand.fa$q <- rand.fa$q + garch_rolling[i,9]%>% as.numeric()
  y <- sample(rand.fa$q,10000)
  
  # 単純モンテカルロ法
  #VaRの計算
  VaR1 <- quantile(y, c(0.01,0.025, 0.05))
  # ESの計算
  Es1 <- c( mean(y[y < VaR1[1]]),   
            mean(y[y < VaR1[2]]),
            mean(y[y < VaR1[3]]))
  
  # 真値と単純モンテカルロ法の結果をまとめる
  cbind(t(VaR1),t(Es1))
  
})
garch_SMC_norm_mu_var <- pforeach::pforeach(i = 1:dim(garch_rolling)[1], .combine = rbind,.cores = 45)({
  
  ## 正規分布からの乱数を取得
  y<-rnorm(n=10000, mean =garch_rolling[i,4]%>% as.numeric(),
           sd=garch_rolling[i,5]%>% as.numeric())
  
  # 単純モンテカルロ法
  #VaRの計算
  VaR1 <- quantile(y, c(0.01,0.025, 0.05))
  # ESの計算
  Es1 <- c( mean(y[y < VaR1[1]]),   
            mean(y[y < VaR1[2]]),
            mean(y[y < VaR1[3]]))
  
  # 真値と単純モンテカルロ法の結果をまとめる
  cbind(t(VaR1),t(Es1))
  
})

garch_SMC_rolling <- cbind(garch_SMC_norm_mu_var,garch_SMC_fsa_mu_var)


#save(list=c("garch_SMC_rolling"),file="data/garch_SMC_rolling.Rdata")
load(file="data/garch_SMC_rolling.Rdata")
str(garch_SMC_rolling)
garch_var99 <- data.frame(logx=df$log_x[251:length(df$log_x)],
                          garch_SMC_rolling[,1],
                          garch_SMC_rolling[,7]) 

str(garch_var99)
check_VaR <- apply(garch_var99, 1,
                   function(x) x[1] < x[c(2,3)]) %>% t()
count_VaR <- colSums(check_VaR)
tmp <- count_VaR/ dim(check_VaR)[1]

es_99 <- data.frame(garch_SMC_rolling[,4],
                    garch_SMC_rolling[,10]) 
eval_ES <- check_VaR * apply(es_99,2,
                             function(x) abs(df$log_x[251:length(df$log_x)] - x)/ -x)
tmp3 <- colMeans(eval_ES)



garch_var975 <- data.frame(logx=df$log_x[251:length(df$log_x)],
                          garch_SMC_rolling[,2],
                          garch_SMC_rolling[,8]) 


check_VaR <- apply(garch_var975, 1,
                   function(x) x[1] < x[c(2,3)]) %>% t()
count_VaR <- colSums(check_VaR)
tmp2 <- count_VaR/ dim(check_VaR)[1]

es_975 <- data.frame(garch_SMC_rolling[,5],
                    garch_SMC_rolling[,11]) 
eval_ES <- check_VaR * apply(es_975,2,
                             function(x) abs(df$log_x[251:length(df$log_x)] - x)/ -x)
tmp4 <- colMeans(eval_ES)

tmp_matrix <- cbind(tmp,tmp2,tmp3,tmp4)

xtable(tmp_matrix,dig=4)
