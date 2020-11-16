bind<-function(df1, df2){
  if (is.null(df1)){
    df1<-df2
  }else{
    df1<-rbindlist(list(df1, df2))
  }
  return(df1)
}

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

min_dist<-function(x, y, points){
  min(sqrt((x-points$x)^2+(y-points$y)^2), na.rm = T)
}

