library(tseries)
#パッケージのインストールと読み込み
#持ってないパッケージはインストールする
targetPackages <- c('zoo', 'xts','Quandl',
                    'quantmod','ggplot2','grid',"reshape2",'scales',
                    'dplyr','moments','xtable','gridExtra','snow',
                    'tidyr',
                    'parallel',"doParallel") 
newPackages <- targetPackages[!(targetPackages %in% installed.packages()[,"Package"])]
if(length(newPackages)) install.packages(newPackages, repos = "http://cran.us.r-project.org")
for(package in targetPackages) library(package, character.only = T)

#データの読み込み
n225 <- read.csv("data/nky.csv",header=TRUE,skip=4)
y <- NULL
#終値(1日の最後の値段)を使う
y$Close <- n225$PX_LAST
#日付データをDate型に変換
y$ymd <- as.POSIXct(n225$Date)
#データフレームにする(行列の列に名前がついているもの)
#ggplotはdata.frameのデータにしか使えないので注意
df <-data.frame(dt=y$ymd, x=y$Close)

#日経平均の対数収益率をplot
df$log_x <- c(NA,diff(log(df$x))*100)
ggplot(df[-1,],aes(dt,log_x))+geom_line()+
  scale_x_datetime(breaks = date_breaks("6 months"))+
  labs(x="Date",y="log return")+
  theme_bw()+
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black"))



# GARCHモデル推定関数
garch_f_opt <- function(data, order = c(1, 1), error = "norm"){
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
  optim(params, garch_op)
}

# garchモデルのsigmaを計算する関数
garch_f <- function(data, params, order = c(1, 1), error = "norm"){
  t <- length(data)
  sigma_t <- rep(0, t)
  if(length(params) != order[1] + order[2] + 1 ){
    stop("The number of parameters is wrong")
  }
  p <- order[1]
  q <- order[2]
  #sigma_tの計算
  #計算できないところまでは、全体のボラティリティを与える
  for(i in 1:p){
    sigma_t[i] <- sum(data^2)/t
  }
  for(i in (p+1):t){
    sigma_t[i] <- exp(params[1]) +
      exp(params[2:(p+1)]) * sigma_t[(i-p):(i-1)] +
      exp(params[(p+2):(p+q+1)]) * data[(i-q):(i-1)]^2
  }
  return(sigma_t)
}


R_fun_res <- garch(df$log_x[-1])
my_fun_res <- garch_f_opt(df$log_x[-1])
my_fun_res2 <- garch_f(df$log_x[-1], my_fun_res$par)
R_fun_res$coef
exp(my_fun_res$par)

plot(R_fun_res$fitted.values[,1] , type="l", main="MLE GARCH volatility ", ylab="", minor.ticks=F)
plot(garch_f(df$log_x[-1], my_fun_res$par) , type="l", main="MLE GARCH volatility ", ylab="", minor.ticks=F)

plot_d  <- data.frame(dt = df$dt[-1],
                R_fun = R_fun_res$fitted.values[,1], my_fun = my_fun_res2 %>% sqrt())

ggplot(plot_d %>% gather(key=fun,value,-dt)) + 
  geom_line(aes(x = dt, y= value,color=fun)) +
  theme_bw()

