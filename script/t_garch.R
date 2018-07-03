# 局度変換を伴うsinh-arcsinh(x)分布の平均，分散，歪度，尖度を計算する関数
ft.sd<- function(df){
  f2 <- function(x) x^2*dt(x, df=df)
  v2 <- integrate( f2, lower=-Inf, upper=Inf)$value
  return(sqrt(v2))
}

dt_stand <- function(x, df){
  mu.sd <- ft.sd(df = df)
  dt(mu.sd*x, df)*mu.sd
}

##多分、元の分布の標準偏差をいじれる分布が必要
dt_stand_sd <- function(x, df, sd){
  mu.sd <- ft.sd(df = df)
  dt(mu.sd*x/sd, df)*mu.sd/sd
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
  params <- rep(0,c(G_beta + G_alpha + 1 + 2))
  params_len <- length(params)
  params[params_len] <- -1
  # garchモデルの対数尤度を計算する関数(FA分布版)  
  garch_op <- function(params){
    df <- exp(params[params_len]) + 2
    #sigma_tの計算
    #計算できないところまでは、全体のボラティリティを与える
    for(i in 1:G_beta){
      sigma_t[i] <- sum(data^2)/t
      garch_loglik[i] <- log(dt_stand_sd(data[i], df = df,sd = sigma_t[i] %>% sqrt()))
    }
    for(i in (G_beta+1):t){
      sigma_t[i] <- exp(params[1]) +
        exp(params[2:(G_beta+1)]) * sigma_t[(i-G_beta):(i-1)] +
        exp(params[(G_beta+2):(G_beta+G_alpha+1)]) * (data[(i-G_alpha):(i-1)] - params[params_len-1])^2
      garch_loglik[i] <- log(dt_stand_sd((data[i] - params[params_len-1]),df = df,sd = sigma_t[i] %>% sqrt()))
    }
    return(-sum(garch_loglik))
  }
  
  out <- optim(params, garch_op)
  
  out
}

#自作のGARCH
#my_fun_res_t <- garch_f_opt_with_para_t(df$log_x[-1])

#分散の計算
my_fun_res2_t <- garch_f_mu(df$log_x[-1], mu = my_fun_res_t$par[4],params = my_fun_res_t$par[1:3], error = "ToS_stand",dis_para = tos_para)
exp(my_fun_res$par[1:3])



t_g_para <- exp(my_fun_res$par[1:3])
t_g_mu <- my_fun_res$par[4]
t_g_df <- exp(my_fun_res$par[5])+2

plot_d  <- data.frame(dt = df$dt[-1],
                      my_fun_res2 = my_fun_res2_mu %>% sqrt(),
                      my_fun_tos2_with_para2_t = my_fun_res2_t %>% sqrt(),
                      my_fun_tos2_with_para2 = my_fun_res2_with_para2_mu %>% sqrt()
                      )


ggplot(plot_d %>% gather(key=fun,value,-dt)) + 
  geom_line(aes(x = dt, y= value,color=fun)) +
  scale_color_hue(name = "関数", labels = c(expression(Norm),expression(ToS),expression(t)))+
  theme_bw() 

