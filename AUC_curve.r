Roc_curve <- function(x) {
  
  
  res <- data.frame()
  for(i in seq(min((x$Sum)),(max((x$Sum))),length.out = 100)){
    #z <- log10(temp$V5+1)>=i
    z <- x$Sum >=i
    a=length(which(x$V4=="Experimental"&z==T))
    b=length(which(x$V4=="Experimental"&z==F))
    sensi = a/(a+b)
    
    c=length(which(x$V4=="random"&z==T))
    d=length(which(x$V4=="random"&z==F))
    speci=d/(c+d)
    res <- rbind(res,c(i,speci,sensi))
  }
  
  colnames(res) <- c("threshold","Specificity","Sensitivity")
  
  a <- ggplot(data=x,aes((Sum),fill = V4))+ 
    geom_density(alpha=0.6) +
    theme_bw() + 
    xlab(paste("Average",feature,"signal in",slop)) + 
    scale_fill_discrete(name = "")
  
  b <- ggplot(res,aes(1-Specificity,Sensitivity,fill = threshold)) + 
    geom_line(col="grey")+ 
    geom_point(col = "grey",shape = 21) +
    theme_bw() +
    scale_fill_continuous(low = "white",high = "red") + 
    scale_x_continuous(name = "False Positive Rate", labels = scales::percent) + 
    scale_y_continuous(name = "True Positive Rate", labels = scales::percent) +
    geom_abline(slope = 1,intercept = 0,lty=2,col="grey")

  
  c <- ggplot(res %>% melt.data.frame(id.vars = c("threshold")),aes(threshold,value,col = variable)) + 
    geom_point() + 
    scale_y_continuous(name = "Value",labels = scales::percent) +
    scale_color_discrete(name="")+
    theme_bw()
  
  grid.arrange(a,c,b)

}


