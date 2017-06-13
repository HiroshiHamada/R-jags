
#関数edu2が計算用，rr2は進学者数の変化に伴う結果をplotする関数

edu2<-function(x,n1,n0){
  #x:進学者数    n1:上層人数　　　n0:下層人数
  df1 <- data.frame(score = rnorm(n1, mean = 52, sd = 1),id = rep(1,n1))
  df0 <- data.frame(score = rnorm(n0, mean = 50, sd = 1),id = rep(0,n0))
  #df1,上層n1人の成績データ      df0:下層n0人の成績データ    
  #score,idは変数名ラベル．結合のため同じ名前を使う
  df <- rbind(df1, df0) 
  #データフレームdfへと，df1とdf0を行結合して格納する
  #同じ列名でないとマッチしないので注意する
  df<-df[order(df$score, decreasing = TRUE), ]  #成績降順で並び替え
  x1=head(df$id, n=x) #上位x人のid抽出
  rr <- (sum(x1)/x)/(1-sum(x1)/x)#相対リスク比.上層進学率/下層進学率  
  return(rr)}

rr2<-function(n1,n0){
  #n1:上層人数　　　n0:下層人数
  end<-n1+n0
  result<-Inf#結果の格納．最初は下層進学者0人のため，infを入れておく
  for(i in 1:end){
  result<- c(result,edu2(i,n1,n0))}#結果をfor文でresultに追加
  #ここで先に定義した関数edu2を使う．進学者数をiで変化させる
    xop<-0:end #plot用にx軸を定義
  plot(xop,result)
}
 
rr2(1000,1000)


rr2(5000,5000)#時間がかかるので注意


