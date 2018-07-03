IS.t.outs <- c()
SMC.t.outs <- c()
parameters <- c()

SMC_t_pre <- function(theta){
  y<-rt(n=10000, df=theta[1], 
                        ncp=theta[2])
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
IS_t_pre <- function(){
  
  rfa1 <- rnorm(10000,mean = VaR.true.FA[1],sd=1)
  rfa25 <- rnorm(10000,mean = VaR.true.FA[2],sd=1)
  rfa5 <- rnorm(10000,mean = VaR.true.FA[3],sd=1)
  
  #weightを計算する
  w1 <- dt(rfa1,df=fit$par2[1],ncp=fit$par2[2])/dnorm(rfa1,mean = VaR.true.FA[1],sd=1)
  w25 <- dt(rfa25,df=fit$par2[1],ncp=fit$par2[2])/dnorm(rfa25,mean = VaR.true.FA[2],sd=1)
  w5 <- dt(rfa5,df=fit$par2[1],ncp=fit$par2[2])/dnorm(rfa5,mean = VaR.true.FA[3],sd=1)
  
  #99%点での計算 100~10000までサンプル数を増やして行う
  out1<-cbind( rfa1,  w1/10000)
  # サンプルの小さい順にならべかえる
  A <- out1[sort.list(out1[,1]),]
  # weightの累積和を並べる
  A <- cbind(A, cumsum(A[,2]))
  # 累積和が0.01に一番近いサンプルが99%VaR
  # v1までのサンプルからES0.01の推定値を求める
  v1 <- A[which.min(abs(A[,3]-0.01)),1]
  es1<- sum(apply(A[1:which.min(abs(A[,3]-0.01)),1:2], 1, prod))/0.01
  out1 <- c(v1, es1)
  
  
  out25<-cbind(rfa25,  w25/10000)
  A <- out25[sort.list(out25[,1]),]
  A <- cbind(A, cumsum(A[,2]))
  v25 <- A[which.min(abs(A[,3]-0.025)),1]
  es25<- sum(apply( A[1:which.min(abs(A[,3]-0.025)),1:2], 1, prod))/0.025
  out25 <- c(v25, es25)
  
  out5<-cbind( rfa5,  w5/10000)
  A <- out5[sort.list(out5[,1]),]
  A <- cbind(A, cumsum(A[,2]))
  v5 <- A[which.min(abs(A[,3]-0.05)),1]
  es5<- sum(apply( A[1:which.min(abs(A[,3]-0.05)),1:2], 1, prod))/0.05
  out5 <- c(v5, es5)
  
  return(out = cbind(t(out1),t(out25),t(out5)))
  
}

IS.t.outs <- c()
SMC.t.outs <- c()
parameters <- c()

for(i in 2:(length(df$log_x)-249)){
  rt <- df$log_x[i:(i+249)]
  rt <- rt[rt!=0]
  
  first_parameter <- c(0, 0)
  if(i != 2){
    first_parameter <- c(old_parameter[1],  old_parameter[2])
  }
  try(fit <- mle.dt(rt,ini=c(first_parameter[1],first_parameter[2])))
  if(first_parameter[1] == old_parameter[1] &first_parameter[2] == old_parameter[2]){
    fit <- mle.dt(rt,ini=c(0,0))
  }
  SMC.t.out <- SMC.t_pre(fit$par2)
  
  
  #真値計算?
  #99%,97.5%,95%の各点に対して，先ほどの関数を用いて求める
  VaR1.fa <- qt(0.01, df=fit$par2[1], 
                ncp=fit$par2[2])
  VaR25.fa <- qt(0.025, df=fit$par2[1], 
                   ncp=fit$par2[2])
  VaR5.fa <- qt(0.05, df=fit$par2[1], 
                  ncp=fit$par2[2])
  VaR.true.FA <- c(VaR1.fa ,VaR25.fa ,VaR5.fa )
  
  
  IS.t.out <- c()
  while(is.null(IS.t.out)){
    try(IS.t.out <- IS_t_pre())
  }
  #-------
  
  
  IS.t.outs <- rbind(IS.t.outs,IS.t.out)
  SMC.t.outs <- rbind(SMC.t.outs,SMC.t.out)
  parameters <- rbind(parameters, fit$par2)
  
  print(i)
  #return(IS.fa.out,IS.norm.outSMC.fa.out,SMC.norm.out)
  old_parameter <- fit$par
}

result_t <- cbind(dt=df$dt[251:length(df$dt)],SMC.t.outs,IS.t.outs,parameters)
#save(list=c("result_t"),file="data/20180701_rolling_result_t.Rdata")

t_roll_var <- data.frame(d_logx=df$log_x[251:length(df$dt)],result_t[,c(2,3,8,10)])
var_matrix <- apply(t_roll_var,1,function(x) x[1] < x[2:5]) %>% t
colSums(var_matrix)/dim(var_matrix)[1]

Ec_99 <- length(df$log_x[251:length(df$dt)]) * 0.01
Ec_97.5 <- length(df$log_x[251:length(df$dt)]) * 0.025

ee <- var_matrix * result_t[,c(5,6,9,11)]

sum((df$log_x[251:length(df$dt)][which(var_matrix[,1] == TRUE)] / ee[,1][which(var_matrix[,1] == TRUE)]  -1)/Ec_99)
sum((df$log_x[251:length(df$dt)][which(var_matrix[,2] == TRUE)] / ee[,2][which(var_matrix[,2] == TRUE)]  -1)/Ec_97.5)
sum((df$log_x[251:length(df$dt)][which(var_matrix[,3] == TRUE)] / ee[,3][which(var_matrix[,3] == TRUE)]  -1)/Ec_99)
sum((df$log_x[251:length(df$dt)][which(var_matrix[,4] == TRUE)] / ee[,4][which(var_matrix[,4] == TRUE)]  -1)/Ec_97.5)




load("data/20180530_rolling_result_useoldpara.Rdata")
colnames(result) <- c("dt","IS_VaR_fa_0.01","IS_ES_fa_0.01",
                      "IS_VaR_fa_0.025","IS_ES_fa_0.025",
                      "IS_VaR_fa_0.05","IS_ES_fa_0.05",
                      "IS_VaR_norm_0.01","IS_ES_norm_0.01",
                      "IS_VaR_norm_0.025","IS_ES_norm_0.025",
                      "IS_VaR_norm_0.05","IS_ES_norm_0.05",
                      "SMC_VaR_fa_0.01","SMC_VaR_fa_0.025","SMC_VaR_fa_0.05",
                      "SMC_ES_fa_0.01","SMC_ES_fa_0.025","SMC_ES_fa_0.05",
                      "SMC_VaR_norm_0.01","SMC_VaR_norm_0.025","SMC_VaR_norm_0.05",
                      "SMC_ES_norm_0.01","SMC_ES_norm_0.025","SMC_ES_norm_0.05",
                      "mu","sigma","lambda","delta")


result_para <- result[,c(1,26,27,28,29)]  %>% data.frame()
load("data/location_time0701")

load("data/20180701_2_rolling_result_useoldpara.Rdata")


colnames(result) <- c("IS_VaR_fa_0.01","IS_ES_fa_0.01",
                      "IS_VaR_fa_0.025","IS_ES_fa_0.025",
                      "IS_VaR_fa_0.05","IS_ES_fa_0.05",
                      "IS_VaR_norm_0.01","IS_ES_norm_0.01",
                      "IS_VaR_norm_0.025","IS_ES_norm_0.025",
                      "IS_VaR_norm_0.05","IS_ES_norm_0.05",
                      "SMC_VaR_fa_0.01","SMC_VaR_fa_0.025","SMC_VaR_fa_0.05",
                      "SMC_ES_fa_0.01","SMC_ES_fa_0.025","SMC_ES_fa_0.05",
                      "SMC_VaR_norm_0.01","SMC_VaR_norm_0.025","SMC_VaR_norm_0.05",
                      "SMC_ES_norm_0.01","SMC_ES_norm_0.025","SMC_ES_norm_0.05")



result_VaR <- result[,c(1,3,5,7,9,11,13,14,15,19,20,21)]
result_ES <- result[,c(2,4,6,8,10,12,16,17,18,22,23,24)]

result_VaR_new <- cbind(result_VaR,result_location_time_var[,c(1,3,5)],
                        result_location_time_es[,c(1,3,5)])
result_ES_new <- cbind(result_ES,
                       result_location_time_var[,c(2,4,5)],
                       result_location_time_es[,c(2,4,6)])
result_VaR_plot99 <- result_VaR_new[,c(1,4,7,10,13)]
plot_d <- cbind(dt = df$dt[251:length(df$dt)],result_VaR_plot99)
colnames(plot_d) <- c("dt","IS_VaR_fa","IS_VaR_norm",
                                 "SMC_VaR_fa","SMC_VaR_norm",
                                 "IS_VaR_fa(shift)")



p1 <- ggplot(plot_d %>% gather(variable,value,-dt)) + geom_path(aes(dt,value,color=variable))+
  theme_bw() + theme(legend.position = "none")+
  theme(axis.title.x = element_text(size=25),axis.title.y = element_text(size=25))+
  theme(axis.text.x = element_text(size=25),axis.text.y = element_text(size=25)) +
  theme(legend.title = element_text(size=15),legend.text = element_text(size=15))

result_ES_plot99 <- result_ES_new[,c(1,4,7,10,13)]


plot_d <- cbind(dt = df$dt[251:length(df$dt)],result_ES_plot99)

colnames(plot_d) <- c("dt","IS_ToS","IS_norm",
                      "SMC_ToS","SMC_norm",
                      "IS_ToS(shift)")

p2 <- ggplot(plot_d %>% gather(variable,value,-dt)) + geom_path(aes(dt,value,color=variable))+
  theme_bw() + theme(legend.position = "bottom")+
  theme(axis.title.x = element_text(size=25),axis.title.y = element_text(size=25))+
  theme(axis.text.x = element_text(size=25),axis.text.y = element_text(size=25)) +
  theme(legend.title = element_text(size=15),legend.text = element_text(size=15))
  
  
  grid.arrange(p1, p2,
               ncol = 1) 

  