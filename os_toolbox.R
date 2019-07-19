## Toolbox for distribution modeling using Biomod2
## @author Yujing Yan

# mask raster
intersect_mask <- function(x){ # @ x:raster stack
  values_x <- getValues(x)
  inter_x <- values_x %*% rep(1,nlayers(x))
  mask <- setValues(subset(x,1),values = (inter_x>0))
  return(mask)
}

# generate evaluation matrix from model output
extract_test <- function (x) {  # @models.evaluation@val
  test <- matrix(nrow = 72, ncol = 5)
  test <- data.frame(test)
  evalname <- rownames(x)
  mod_eval_names <- c("name")
  for (i in 1: (nrow(x)/4)){                        
    mod_eval_names[i] <- evalname[4*i-3]
    test[i,2] <- x[(4*i-3),1]
    test[i,3] <- x[(4*i-3),2]
  }
  evalname1 <- matrix(ncol = 5, nrow = nrow(test))
  for (i in 1:nrow(test)){
    evalname1[i,] <- unlist(strsplit(mod_eval_names[i],"\\."))
  }
  test[,1] <- evalname1[,3]
  test[,4:5] <- evalname1[,4:5]
  colnames(test) <- c('Model', 'TSS', 'ROC','RUN','PA')
  return(test)
}

# convert evaluation matrix from wide to long
GetEvaluate <- function(eva){
  eval_test <- extract_test(eva)
  eval <- melt(eval_test[,1:3], id.vars = "Model",    
               measure.vars = c("TSS", "ROC"))
  colnames(eval) <- c("Model", "Eval", "Value")
  return (eval)
}

#e <- extent(ra)
#i<-"H:/distribution/future climate/5026/ip/ip26bi507.asc"
cropAndwrite <- function(i){
  cat("\n",i)
  x <- raster(i)
  y <- crop(x,e)
  writeRaster(y,i, format='ascii',overwrite = T)
}

getmean <- function(c,files,path){
  gn <- files[grep(c,files)]
  rg <- lapply(gn,importraster)
  rg <- stack(rg)
  rasters10 <- sum(rg)
  writeRaster(rasterm,paste(path,"bio",c,sep=""),format='ascii')
}

plotelec <- function(x,year){
  x <- x[,-1]
  x[4,] <- x[2,]
  x[5,] <- x[1,]
  x[6,] <- x[3,]
  x <- x[-1,]
  x <- x[-1,]
  x <- x[-1,]
  rownames(x) <- c("3", "2","1")
  colnames(x) <- c("<3000","3000-3500","3500-4000","4000-4500","4500-5000",">5000")
  x <- as.matrix(x)
  
  png(paste("os1_plot/ele_os1_",year,"1.jpg",sep=""),width=1200,height=900,bg="white")
  par(mfrow = c(1,1),cex.axis = 3.5, cex.lab=4, mar = c(20,11,5,3))#,mar = c(20,11,5,3))#,mar=c(3,3,3,3)
  
  barplot(x,col=c("grey22","grey48","grey69"),border=NA,ylim=c(0,45),las=2,xaxt="n",width=0.9)#,mgp=c(17.5,2,0),las=2,tcl=-0.8)
  
  text(1:6, par("usr")[3]-3, srt = 45, adj = 1,
       labels = colnames(plotele), xpd = TRUE,cex=3.5) 
  title(xlab="Elevation(m)",mgp=c(15.5,2,0))
  title(ylab="Percentage (%)",mgp=c(7,3.5,0)) 
  legend(x=4.2,y=45 ,legend=c("Predicted gain","Predicted loss","Predicted stable"),fill=c("grey69","grey48","grey22"),border=NA,bty="n",x.intersp =0.7, y.intersp = 0.9,cex=3.5)
  dev.off()
}