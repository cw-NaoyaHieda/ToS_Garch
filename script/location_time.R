# ローリング推定したパラメータを使ってlocationsfihtした重点サンプリングでのVaRの推定
# 時間かかるので回すときは注意

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

cl <- makeCluster(detectCores()-1)  # クラスタの作成
#clusterExport(cl,c('s_inverse','dfas2','dfas','Cf','Sf'))


result_location_time_var <- pforeach::pforeach(i = 1:dim(result_para)[1], .combine = rbind)({
  result_para_now <- result_para[i,] %>% as.numeric()
  VaR1.fa <- qfas(0.01, mu=result_para_now[2], sigma=result_para_now[3],
                  lambda=result_para_now[4], delta = result_para_now[5])
  VaR25.fa <- qfas(0.025, mu=result_para_now[2], sigma=result_para_now[3],
                  lambda=result_para_now[4], delta = result_para_now[5])
  VaR5.fa <- qfas(0.05, mu=result_para_now[2], sigma=result_para_now[3],
                  lambda=result_para_now[4], delta = result_para_now[5])
  
  VaR.true.FA <- c(VaR1.fa ,VaR25.fa ,VaR5.fa )
  
  mu.val1 <- VaR.true.FA[1]
  mu.val25 <- VaR.true.FA[2]
  mu.val5 <- VaR.true.FA[3]

### 次の関数を新しく用意した。
rIS_SIR2_ <- function(n, par, par2){
  ## 正規分布を提案分布に par2 = c( mu, sigma) 
  q <- rnorm(n, mean=par2[1], sd=par2[2])
  ## 重み
  #hosts <- rep('localhost',12)
  #hosts <- 7
  #cl <- makeCluster(hosts, "SOCK")
  w <- sapply(q, 
                 dfas2, mu = par[1], sigma= par[2] , lambda=par[3], delta = par[4]
  )/dnorm(q, mean=par2[1], sd=par2[2])
  w <- w/sum(w)
  ## resample
  q.resample <- Resample1(q, weight=w, NofSample = n)
  list( q=q.resample, w=w)
}

rfa.IS.1<-rIS_SIR2_(n=20000, par=c( mu.val1, result_para_now[3:5]),
                    par2 = c(mu.val1, sd(rt))) 
rfa.IS.25<-rIS_SIR2_(n=20000, par=c( mu.val25, result_para_now[3:5]),
                    par2 = c(mu.val1, sd(rt))) 
rfa.IS.5<-rIS_SIR2_(n=20000, par=c( mu.val5, result_para_now[3:5]),
                    par2 = c(mu.val1, sd(rt))) 

rfa1 <- sample(rfa.IS.1$q, 10000)
rfa25 <- sample(rfa.IS.25$q, 10000)
rfa5 <- sample(rfa.IS.5$q, 10000)



# 重点サンプリングの関数の関数も新しいものにした。
f <- function(x, par, par2){
  #par2 f # par g
  dfas2(x, mu=par[1], sigma=par[2],
        lambda=par[3], delta=par[4])/
    dfas2(x, mu=par2[1], sigma=par[2],
          lambda=par[3], delta=par[4])}



#weightを計算する
w1 <- sapply(rfa1,f, par = result_para_now[2:5], par2=c(mu.val1, result_para_now[3:5]))
w25 <- sapply(rfa25,f, par = result_para_now[2:5], par2=c(mu.val25, result_para_now[3:5]))
w5 <- sapply(rfa5,f, par = result_para_now[2:5], par2=c(mu.val5, result_para_now[3:5]))

#99%点での計算 100~10000までサンプル数を増やして行う
N<-10000
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

data.frame(v1,es1,v25,es25,v5,es5)
})

result_location_time_es <- pforeach::pforeach(i = 1:dim(result_para)[1], .combine = rbind)({
  result_para_now <- result_para[i,] %>% as.numeric()
  ES1.fa <- find.ES(p=0.01, par=c(mu=result_para_now[2], sigma=result_para_now[3],
                    lambda=result_para_now[4], delta = result_para_now[5]))
  ES25.fa <- find.ES(p=0.025, par=c(mu=result_para_now[2], sigma=result_para_now[3],
                                    lambda=result_para_now[4], delta = result_para_now[5]))
  ES5.fa <- find.ES(p=0.05, par=c(mu=result_para_now[2], sigma=result_para_now[3],
                                  lambda=result_para_now[4], delta = result_para_now[5]))
  #まとめる
  ES.true.FA <- c(ES1.fa ,ES25.fa ,ES5.fa )
  
  mu.val1 <- ES.true.FA[1]
  mu.val25 <- ES.true.FA[2]
  mu.val5 <- ES.true.FA[3]
  
  ### 次の関数を新しく用意した。
  rIS_SIR2_ <- function(n, par, par2){
    ## 正規分布を提案分布に par2 = c( mu, sigma) 
    q <- rnorm(n, mean=par2[1], sd=par2[2])
    ## 重み
    #hosts <- rep('localhost',12)
    #hosts <- 7
    #cl <- makeCluster(hosts, "SOCK")
    w <- sapply(q, 
                dfas2, mu = par[1], sigma= par[2] , lambda=par[3], delta = par[4]
    )/dnorm(q, mean=par2[1], sd=par2[2])
    w <- w/sum(w)
    ## resample
    q.resample <- Resample1(q, weight=w, NofSample = n)
    list( q=q.resample, w=w)
  }
  
  rfa.IS.1<-rIS_SIR2_(n=20000, par=c( mu.val1, result_para_now[3:5]),
                      par2 = c(mu.val1, sd(rt))) 
  rfa.IS.25<-rIS_SIR2_(n=20000, par=c( mu.val25, result_para_now[3:5]),
                       par2 = c(mu.val1, sd(rt))) 
  rfa.IS.5<-rIS_SIR2_(n=20000, par=c( mu.val5, result_para_now[3:5]),
                      par2 = c(mu.val1, sd(rt))) 
  
  rfa1 <- sample(rfa.IS.1$q, 10000)
  rfa25 <- sample(rfa.IS.25$q, 10000)
  rfa5 <- sample(rfa.IS.5$q, 10000)
  
  
  
  # 重点サンプリングの関数の関数も新しいものにした。
  f <- function(x, par, par2){
    #par2 f # par g
    dfas2(x, mu=par[1], sigma=par[2],
          lambda=par[3], delta=par[4])/
      dfas2(x, mu=par2[1], sigma=par[2],
            lambda=par[3], delta=par[4])}
  
  
  N<-10000
  #weightを計算する
  w1 <- sapply(rfa1,f, par = result_para_now[2:5], par2=c(mu.val1, result_para_now[3:5]))
  w25 <- sapply(rfa25,f, par = result_para_now[2:5], par2=c(mu.val25, result_para_now[3:5]))
  w5 <- sapply(rfa5,f, par = result_para_now[2:5], par2=c(mu.val5, result_para_now[3:5]))
  
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
  
  data.frame(v1,es1,v25,es25,v5,es5)
})


stopCluster(cl)

#save(list = c("result_location_time_var",
#              "result_location_time_es"),
#     file="/home/naoya/Desktop/ToS_Garch/data/location_time")

#save(list = c("result_location_time_var",
#                            "result_location_time_es"),
#                   file="/home/naoya/Desktop/ToS_Garch/data/location_time0701")

#save(list = c("result_location_time_var",
#                            "result_location_time_es"),
#                   file="/home/naoya/Desktop/ToS_Garch/data/location_time0701_2")
