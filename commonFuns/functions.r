bind<-function(df1, df2){
  if (is.null(df1)){
    df1<-df2
  }else{
    df1<-rbindlist(list(df1, df2), use.names=TRUE)
  }
  return(df1)
}

bind_dplyr<-function(df1, df2){
  if (is.null(df1)){
    df1<-df2
  }else{
    df1<-dplyr::bind_rows(df1, df2)
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


NDquntil <- function(nD, level) {
  n <- floor(nD * level)
  if (n > nD) 
    n <- nD
  return(n)
}
in_Ellipsoid <- stats::qchisq(0.95, 2)
