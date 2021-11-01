dfxx<-NULL
for (k in c(seq(0, 0.9, by=0.1), seq(1, 10, by=1))){
  for (i in c(1:100)){
    
    x<-runif(100) * 100
    y<-x + k * runif(100) * 100
    
    r2<-caret::postResample(x, y)
    item<-data.frame(k=k, cor=cor(x, y), rmse=r2[1], mae=r2[3], r2=r2[2])
    if (is.null(dfxx)){
      dfxx<-item
    }else{
      dfxx<-rbind(dfxx, item)
    }
  }
}
library(ggplot2)
ggplot(dfxx)+geom_point(aes(x=rmse, y=r2, color=cor))

ggplot(dfxx)+geom_point(aes(x=rmse, y=r2, color=cor))
ggplot(dfxx)+geom_point(aes(x=cor, y=r2, color=cor))

ggplot(evaluated_metrics_all)+geom_point(aes(x=RMSE, y=Rsquared, color=factor(formulas)))
ggplot(evaluated_metrics_all)+geom_point(aes(x=cor, y=Rsquared, color=factor(formulas)))

plot(dfxx$cor, dfxx$r2)
plot(dfxx$rmse, dfxx$mae)
plot(dfxx$rmse, dfxx$r2)
plot(dfxx$mae, dfxx$r2)

