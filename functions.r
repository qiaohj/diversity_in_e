bind<-function(df1, df2){
  if (is.null(df1)){
    df1<-df2
  }else{
    df1<-dplyr::bind_rows(df1, df2)
  }
  return(df1)
}