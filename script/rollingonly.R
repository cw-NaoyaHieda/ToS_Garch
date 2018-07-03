# 重点サンプリングの結果がなんか違うのでやり直してみる

IS.fa_pre <- function(par){
  
  f <- function(x, theta, par){
    exp(theta*x)*dfas2(x, mu=par[1], sigma=par[2],
                       lambda=par[3], delta=par[4])
  }
  #weightを計算するためにMを計算する．99%,97.5%,95%をまとめて行う
  M1 <- integrate(f, -30, 30, theta=theta.val1, par=par)$value
  M25 <- integrate(f, -30, 30, theta=theta.val25, par=par)$value
  M5 <- integrate(f, -30, 30, theta=theta.val5, par=par)$value
  
  #weightを計算する
  w1 <- exp(-theta.val1*rfa1)*M1
  w25 <- exp(-theta.val25*rfa25)*M25
  w5 <- exp(-theta.val5*rfa5)*M5
  
  #99%点での計算 100~10000までサンプル数を増やして行う
  out1<-cbind( rfa1,  w1/length(w1))
  # サンプルの小さい順にならべかえる
  A <- out1[sort.list(out1[,1]),]
  # weightの累積和を並べる
  A <- cbind(A, cumsum(A[,2]))
  # 累積和が0.01に一番近いサンプルが99%VaR
  # v1までのサンプルからES0.01の推定値を求める
  v1 <- A[which.min(abs(A[,3]-0.01)),1]
  es1<- sum(apply(A[1:which.min(abs(A[,3]-0.01)),1:2],1,prod))/0.01
  out1 <- c(v1, es1)
  
  
  out25<-cbind(rfa25,  w25/length(w25))
  A <- out25[sort.list(out25[,1]),]
  A <- cbind(A, cumsum(A[,2]))
  v25 <- A[which.min(abs(A[,3]-0.025)),1]
  es25<- sum(apply(A[1:which.min(abs(A[,3]-0.025)),1:2],1,prod))/0.025
  out25 <- c(v25, es25)
  
  out5<-cbind( rfa5,  w5/length(w5))
  A <- out5[sort.list(out5[,1]),]
  A <- cbind(A, cumsum(A[,2]))
  v5 <- A[which.min(abs(A[,3]-0.05)),1]
  es5<- sum(apply(A[1:which.min(abs(A[,3]-0.05)),1:2],1,prod))/0.05
  out5 <- c(v5, es5)
  
  return(out = cbind(t(out1),t(out25),t(out5)))
  
}


rIS_SIR <- function(n, par, par2, theta){
  ## 正規分布を提案分布に
  q <- rnorm(n,mean=par2,sd=15)
  # 重点分布の密度関数の分子
  f <- function(x, theta, par){
    exp(theta*x)*dfas2(x, mu=par[1], sigma=par[2],
                       lambda=par[3], delta=par[4])
  }
  # 分母
  M <- integrate(f, -30, 30, theta=theta, par=par %>% as.numeric())$value
  ### 指数変換した重点分布の密度関数  
  d.IS <- function(x, theta,par){
    # 重点分布の密度関数の分子
    f <- function(x, theta, par){
      exp(theta*x)*dfas2(x, mu=par[1], sigma=par[2],
                         lambda=par[3], delta=par[4])
    }
    # 分母
    return( f(x, theta, par)/M )
  }
  ## 重み
  w <- sapply(q, 
              d.IS, theta =theta, par=par) /dnorm(q, mean=par2, sd=15)
  w <- w/sum(w)
  ## resample
  q.resample <- Resample1(q, weight=w, NofSample = n)
  list( q=q.resample, w=w)
}

#Fsaからの乱数 SIR並列処理
rfa_SIR<- function(n, mu, sigma, lambda, delta)
{
  ## 正規分布を提案分布に
  q <- rnorm(n,mean=mu %>% as.numeric(),sd=5*sigma %>% as.numeric())
  ## 重み
  w <- sapply(q, dfas2, mu=mu, sigma=sigma, lambda=lambda, delta=delta) %>% as.numeric()/
    dnorm(q, mean=mu %>% as.numeric(), sd=5*sigma %>% as.numeric()) %>% as.numeric()
  ## 合計が1になるように重みを基準化
  w <- w/sum(w)
  ## 重みに従ってresample
  q.resample <- Resample1(q, weight=w, NofSample = n)
  list(q,q=q.resample, w=w)
}
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


cl <- makeCluster(detectCores()-1)  # クラスタの作成
registerDoParallel(cl)
cl_l <- detectCores()-1

result <- pforeach::pforeach(i = 1:dim(result_para)[1], .combine = rbind)({
  result_para_now <- result_para[i,]%>% as.numeric()
  #真値計算?
  #99%,97.5%,95%の各点に対して，先ほどの関数を用いて求める
  VaR1.fa <- qfas(0.01, mu=result_para_now[2], sigma=result_para_now[3],
                  lambda=result_para_now[4], delta = result_para_now[5])
  VaR25.fa <- qfas(0.025, mu=result_para_now[2], sigma=result_para_now[3],
                   lambda=result_para_now[4], delta = result_para_now[5])
  VaR5.fa <- qfas(0.05, mu=result_para_now[2], sigma=result_para_now[3],
                  lambda=result_para_now[4], delta = result_para_now[5])
  VaR.true.FA <- c(VaR1.fa ,VaR25.fa ,VaR5.fa )
  #単純モンテカルロ
  SMC.fa.out <- SMC.fa_pre(result_para_now[-1])
  
  #--------------
 
  # 99%,97.5%,95%それぞれのVaRと平均が一致するthetaを取得
  theta.val1<- find.theta(0.01, result_para_now[-1])
  theta.val25<- find.theta(0.025, result_para_now[-1])
  theta.val5<- find.theta(0.05, result_para_now[-1])
  
  out.fa<-c()
  rfa.IS.1<-rIS_SIR(n=20000, par=result_para_now[-1], par2=VaR.true.FA[1], theta=theta.val1)
  rfa.IS.25<-rIS_SIR(n=20000, par=result_para_now[-1], par2=VaR.true.FA[2], theta=theta.val25)
  rfa.IS.5<-rIS_SIR(n=20000, par=result_para_now[-1], par2=VaR.true.FA[3], theta=theta.val5)
  # サンプリングしたものを入力としてFA分布の重点サンプリングを行う
  rfa1 <- sample(rfa.IS.1$q, 10000)
  rfa25 <- sample(rfa.IS.25$q, 10000)
  rfa5 <- sample(rfa.IS.5$q, 10000)
  
  #clusterExport(cl,list("rfa1","rfa25","rfa5"))
  #IS.fa.out <- NULL
  #重点サンプリング
  #while(is.null(IS.fa.out)){
  IS.fa.out <- IS.fa_pre(result_para_now[-1])
  #}
  
  #-------
  rt <- df$log_x[c((i+1):(i+250))]
  theta <- c(mean(rt), sd(rt))
  SMC.norm.out <- SMC.norm_pre(theta)
  IS.norm.out <- IS.norm_pre(theta)
  
  c(IS.fa.out,IS.norm.out,SMC.fa.out,SMC.norm.out)
  
  
  
})

stopCluster(cl)
#result <- cbind(dt=df$dt[251:length(df$dt)],IS.fa.outs,IS.norm.outs,SMC.fa.outs,SMC.norm.outs,parameters)
#save(list=c("result"),file="data/20180629_rolling_result_useoldpara.Rdata")
#save(list=c("result"),file="data/20180701_rolling_result_useoldpara.Rdata")


result_tmp <- c()
for(i in 1:dim(result_para)[1]){
  rt <- df$log_x[(i+1):(i+250)]
  theta <- c(mean(rt), sd(rt))

  SMC.norm.out <- SMC.norm_pre(theta)
  IS.norm.out <- IS.norm_pre(theta)
  result_tmp <- rbind(result_tmp,data.frame(SMC.norm.out,IS.norm.out))
  print(i)
}
load("data/20180629_rolling_result_useoldpara.Rdata")


result <- cbind(result[,c(1:6)],result_tmp[,c(7:12)],result[,c(13:18)],result_tmp[,c(1:6)])

#save(list=c("result"),file="data/20180630_rolling_result_useoldpara.Rdata")









result_only_IS_fa <- pforeach::pforeach(i = 1:dim(result_para)[1], .combine = rbind, .cores = 45)({
#result_only_IS_fa <- c()
  
#for(i in 1:dim(result_para)[1]){
  par <- result_para[i,-1]%>% as.numeric()
  result_para_now <- result_para[i,]%>% as.numeric()
  #99%,97.5%,95%の各点に対して，先ほどの関数を用いて求める
  VaR1.fa <- qfas(0.01, mu=result_para_now[2], sigma=result_para_now[3],
                  lambda=result_para_now[4], delta = result_para_now[5])
  VaR25.fa <- qfas(0.025, mu=result_para_now[2], sigma=result_para_now[3],
                   lambda=result_para_now[4], delta = result_para_now[5])
  VaR5.fa <- qfas(0.05, mu=result_para_now[2], sigma=result_para_now[3],
                  lambda=result_para_now[4], delta = result_para_now[5])
  VaR.true.FA <- c(VaR1.fa ,VaR25.fa ,VaR5.fa )
  
  theta.val1<- find.theta(0.01, par)
  theta.val25<- find.theta(0.025, par)
  theta.val5<- find.theta(0.05, par)
  
  rfa.IS.1<-rIS_SIR(n=20000, par=result_para_now[-1], par2=VaR.true.FA[1], theta=theta.val1)
  rfa.IS.25<-rIS_SIR(n=20000, par=result_para_now[-1], par2=VaR.true.FA[2], theta=theta.val25)
  rfa.IS.5<-rIS_SIR(n=20000, par=result_para_now[-1], par2=VaR.true.FA[3], theta=theta.val5)
  
  f <- function(x, theta, par){
    exp(theta*x)*dfas2(x, mu=par[1], sigma=par[2],
                       lambda=par[3], delta=par[4])
  }
  #weightを計算するためにMを計算する．99%,97.5%,95%をまとめて行う
  M1 <- integrate(f, -30, 30, theta=theta.val1, par=par)$value
  M25 <- integrate(f, -30, 30, theta=theta.val25, par=par)$value
  M5 <- integrate(f, -30, 30, theta=theta.val5, par=par)$value
  rfa1 <- sample(rfa.IS.1$q, 10000)
  rfa25 <- sample(rfa.IS.25$q, 10000)
  rfa5 <- sample(rfa.IS.5$q, 10000)
  
  #weightを計算する
  w1 <- exp(-theta.val1*rfa1)*M1
  w25 <- exp(-theta.val25*rfa25)*M25
  w5 <- exp(-theta.val5*rfa5)*M5
  
  #99%点での計算 100~10000までサンプル数を増やして行う
  out1<-cbind(rfa1,  w1/N)
  # サンプルの小さい順にならべかえる
  A <- out1[sort.list(out1[,1]),]
  # weightの累積和を並べる
  A <- cbind(A, cumsum(A[,2]))
  # 累積和が0.01に一番近いサンプルが99%VaR
  v1 <- A[which.min(abs(A[,3]-0.01)),1]
  # v1までのサンプルからES0.01の推定値を求める
  if(is.null(dim(A[1:which.min(abs(A[,3]-0.01)),]))){
    es1 <- sum(prod(A[1:which.min(abs(A[,3]-0.01)),c(1:2)]))/0.01
  }else{
    es1 <- sum(apply(A[1:which.min(abs(A[,3]-0.01)),1:2],1,prod))/0.01
  }
  
  out25<-cbind(rfa25,  w25/N)
  # サンプルの小さい順にならべかえる
  A <- out25[sort.list(out25[,1]),]
  # weightの累積和を並べる
  A <- cbind(A, cumsum(A[,2]))
  
  v25<-A[which.min(abs(A[,3]-0.025)),1]
  if(is.null(dim(A[1:which.min(abs(A[,3]-0.025)),]))){
    es25 <- sum(prod(A[1:which.min(abs(A[,3]-0.025)),c(1:2)]))/0.025
  }else{
    es25<- sum(apply(A[1:which.min(abs(A[,3]-0.025)),1:2],1,prod))/0.025
  }
  
  
  out5<-cbind(rfa5,  w5/N)
  # サンプルの小さい順にならべかえる
  A <- out5[sort.list(out5[,1]),]
  # weightの累積和を並べる
  A <- cbind(A, cumsum(A[,2]))
  
  v5<-A[which.min(abs(A[,3]-0.05)),1]
  if(is.null(dim(A[1:which.min(abs(A[,3]-0.05)),]))){
    es5 <- sum(prod(A[1:which.min(abs(A[,3]-0.05)),c(1:2)]))/0.05
  }else{
    es5<- sum(apply(A[1:which.min(abs(A[,3]-0.05)),1:2],1,prod))/0.05
  }
  
  c(v1,es1,v25,es25,v5,es5)
  
  #-------
  #result_only_IS_fa <- rbind(result_only_IS_fa,IS.fa.out)
  #print(i)
  
})

#save(list=c("result_only_IS_fa"),file="data/result_only_IS_fa")
load("data/20180701_rolling_result_useoldpara.Rdata")

result <- cbind(result_only_IS_fa,result[,c(-1:-6)])
#save(list=c("result"),file="data/20180701_2_rolling_result_useoldpara.Rdata")

