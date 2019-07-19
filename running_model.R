library(biomod2)
library(maptools)
library(raster)
library(rgdal)
library(sp)
library(spatial)
library(reshape2)

wd = ""
setwd(wd)
source("biomod_toolbox.R")
#### 0.Read in data ####
proj <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 ")
newpro <- CRS("+proj=aea +lat_1=25 +lat_2=47 +lat_0=30 +lon_0=110 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs ")
envList <- list.files("project/",full.names=T) # environment data list
envNames <- list.files("project/",full.names=F) # environment data name list
environ <- stack(lapply(envList,raster,crs=proj)) # read in environment data
imExpl <- intersect_mask(environ)
## keep only all cells that are defined for all layers
maskEx <- mask(environ, imExpl)
myExpl <- stack(maskEx)
myRespXY <- read.csv("Occ.csv",stringsAsFactors = F)  # read in occurrence data
myRespXY <- myRespXY[,2:3]
colnames(myRespXY) <- c("Long","Lat")
myResp <- as.numeric(rep(1,nrow(myRespXY)))
myRespName <- "Os1"

#### 1. Running model ####
### Initialisation
n_absences <- 10*sum(myResp==1,na.rm=TRUE) # number of randomly generated absence points
BiomodData1 <- BIOMOD_FormatingData(resp.var = myResp,
                                     expl.var = myExpl,
                                     resp.xy = myRespXY,
                                     resp.name = myRespName,
                                     PA.nb.rep = 2,
                                     PA.nb.absences = n_absences,
                                     PA.strategy = 'random')
### Options definition
#myBiomodOption <- BIOMOD_ModelingOptions()
### Modelling
BiomodModelOut1 <- BIOMOD_Modeling(BiomodData1,
                                    models = c('GAM','MARS','MAXENT','RF','CTA','GLM'),
                                    #models.options = myBiomodOption,
                                    NbRunEval=5,
                                    DataSplit=75,
                                    Yweights=NULL,
                                    VarImport=3,
                                    models.eval.meth = c('TSS','ROC'),
                                    SaveObj = TRUE,
                                    rescal.all.models = FALSE,
                                    modeling.id="allchina")
### save models evaluation scores and variables importance 
var_imp <- as.matrix(data.frame(BiomodModelOut1@variables.importances@val))
var_mean <- data.frame(apply(var_imp,1,mean))

# select environment factor that contribute more to the model and rebuild the model
subclim <- as.numeric(which(var_mean[,1]>0.05))
print(length(subclim))
write.csv(var_mean,file="var_mean1.txt")

myExpl <- unstack(myExpl)
myExpl_sub1 <- stack(myExpl[subclim])

BiomodData_sub1 <- BIOMOD_FormatingData(resp.var = myResp,
                                     expl.var = myExpl_sub,
                                     resp.xy = myRespXY,
                                     resp.name = myRespName,
                                     PA.nb.rep = 2,
                                     PA.nb.absences = n_absences,
                                     PA.strategy = 'random')
### Options definition
#myBiomodOption <- BIOMOD_ModelingOptions()
### Modelling
BiomodModelOut_sub1 <- BIOMOD_Modeling(BiomodData_sub1,
                                    models = c('GAM','MARS','MAXENT','RF','CTA','GLM'),
                                    #models.options = myBiomodOption,
                                    NbRunEval=5,
                                    DataSplit=75,
                                    Yweights=NULL,
                                    VarImport=3,
                                    models.eval.meth = c('TSS','ROC'),
                                    SaveObj = TRUE,
                                    rescal.all.models = F,
                                    modeling.id="allchina_sub")

### Building ensemble-models
BiomodEM1 <- BIOMOD_EnsembleModeling(modeling.output = BiomodModelOut_sub1,
                                        chosen.models = 'all',
                                        em.by = "all",
                                        eval.metric = c('ROC'),
                                        eval.metric.quality.threshold = c(0.95),
                                        prob.mean = F,
                                        prob.cv = F,
                                        prob.ci = F,
                                        #prob.ci.alpha = 0.05,
                                        prob.median = F,
                                        committee.averaging = T,
                                        prob.mean.weight = T,
                                        prob.mean.weight.decay = 'proportional' )

#### 2.Project to current environment ####
BiomodEF1 <- BIOMOD_EnsembleForecasting(
  BiomodEM,
  #projection.output = BiomodProj5,
  new.env = myExpl_sub,
  proj.name = "current",
  binary.meth = 'ROC',
  total.consensus = TRUE)

envList <- list.files("project/",full.names=T)
environc <- lapply(envList,raster,crs=proj)  #current environment

for(n in 1:2){
  sp.n <- spnames[n]
  cat("\n",sp.n,"model...")   
  listOut <- list.files(paste(sp.n,"/",sep=""),full.names=T)
  emOut <- listOut[grepl("ensemble",listOut)==T]
  #cat("\n",emOut)
  emModelOut <- load(emOut)
  myBiomodEM <- get(emModelOut)
  rm(list = c(emModelOut, 'emModelOut'))
  ### Do ensemble-models projections on current varaiable
  subclim <- get(paste("subclim",n,sep=""))
  #cat("\n",subclim)
  myExplc <- stack(environc[subclim])
  BiomodEF <- BIOMOD_EnsembleForecasting(
    EM.output = myBiomodEM,
    new.env = myExplc,
    proj.name = "current",
    binary.meth = 'ROC',
    total.consensus = TRUE)
}

#### 3. Evaluate model and plot result ####
eva <- t(as.matrix(data.frame(BiomodModelOut_sub1@models.evaluation@val)))
eval <- GetEvaluate(eva)
boxplot(Value~Model+Eval, data = eval)

eva_os1 <- eva   #using sub out data
eval_os1 <- eval
write.csv(eval_os1, file="eval_os1.txt")

colb <- c(rep(c("white","black"),6))
pchs <- c(rep(1,6),rep(19,6))
borb <- c(rep(c("black","white"),6))
wcol <- c(rep("black",12))

xlabb <- c("CTA","GAM","GLM","MARS","MaxEnt","RF")
eval_os1t <- eval_os1[1:72,]

png(filename="evalos1.jpg",width=1200,height=900,bg="white")
par(mfrow = c(1,1),cex.axis = 3.5, cex.lab=4,mar = c(15,12,5,3),mgp = c(7,1,0))#,mar=c(3,3,3,3)
boxplot(Value~Eval+Model, data = eval_os1,ylab="Model Accuracy",ylim=c(0.6,1),border=borb,col=colb,whiskcol=wcol,staplecol=wcol,outpch=pchs,outcol=wcol,xaxt="n",las=2)#,col=colb)#,las=2,tcl=-0.8)#names=xlabb,mgp=c(10,9,0))
axis(1,at=c(1.5,3.5,5.5,7.5,9.5,11.5),labels=F) #
text(c(1.5,3.5,5.5,7.5,9.5,11.5), par("usr")[3]-0.035, srt = 0, 
     labels = xlabb, xpd = TRUE,cex=3.5)
#legend(x=10.3,y=0.65 ,legend=c("AUC", "TSS"),fill=c("black","white"),x.intersp =0.7, y.intersp = 0.9,cex=3.5,bty="n")
dev.off()
wilcox.test(x = eval_9adn$Value, y = eval_18ad$Value)

#### 4. Project to future ####
#_ 4.1 Calculate future climate change(temperature and rainfall) ####
cnames <- c("02","03","08","15","16","17","18") # environment variables
modelname <- c("bc","he","ip","mg","no","average")
futuretime <- c("5026","5085","7026","7085")
# envList: current env list
#read in data and summarize the relative change,write to raster and table
for (n in 1:7){
  cat("\n",envList[[n]])
  rasterc <- raster(envList[[n]])
  for (time.n in futuretime[1:4]){
    for (model.n in modelname[1:6]){
      pathf <- paste("future climate/",time.n,"/",model.n,"/",sep="")
      Listf <- list.files(pathf,full.names=T)
      rasterf <- raster(Listf[[n]])
      rasterdiffer <- (rasterf-rasterc)/rasterc
      ra.name <- paste(time.n,model.n,cnames[n],sep="")
      cat("\n",Listf[[n]])
      cat("\n",ra.name)
      writeRaster(rasterdiffer, paste("future climate/climatedivide/", ra.name, sep=""), format="ascii",overwrite=T)
      differ_freq <- freq(rasterdiffer)
      capture.output(differ_freq,
                     file=file.path(paste("future climate/differdivide/",ra.name,".txt", sep="")))
    }
  }
}

#summarize the difference and write to txt files
for (factor.n in cnames){
  temp_matrix <- matrix(nrow=24,ncol=5)
  for (time.n in futuretime){
    ntime <- which(futuretime[]==time.n)
    for (model.n in modelname){
      nmodel =  which(modelname[]==model.n)
      n = nmodel+6*(ntime-1)
      temp_name <- paste(time.n,model.n,factor.n,sep="")
      cat("\n",temp_name)
      temp_txt <- read.table(paste("future climate/differ/",temp_name,".txt",sep=""),fill=T)
      len <- nrow(temp_txt)
      temp_txt <- temp_txt[-len,]
      sumcount <- sum(temp_txt$count)
      temp_matrix[n,1] <- max(temp_txt$value)  #get the max value
      temp_matrix[n,2] <- min(temp_txt$value)  #get the min value
      temp_matrix[n,3] <- sum(apply(temp_txt,1,cumprod)[2,])/sumcount  #mean
      temp_matrix[n,4] <- temp_txt[which(temp_txt$count == max(temp_txt$count)),1] #mode
      temp_matrix[n,5] <- max(temp_txt$count)/sumcount  #percentage of mode
    }
  }
  rownames(temp_matrix) <- rep(modelname,times=4)
  colnames(temp_matrix) <- c("max","min","mean","mode","percent")
  capture.output(temp_matrix,
                 file=file.path(paste("future climate/",factor.n,".txt", sep="")))
}

#_ 4.2 Project distribution to future ####
sp.n <- "Os1"
listOut <- list.files(paste(sp.n,"/",sep=""),full.names=T)
emOut <- listOut[grepl("ensemble",listOut)==T]
emModelOut <- load(emOut)
myBiomodEM <- get(emModelOut)
rm(list = c(emModelOut, 'emModelOut'))
for (time.n in futuretime[1:4]){
  cat("\n",time.n)
  for (model.n in modelname){
    cat("\n",model.n)
    pathf <- paste("future climate/",time.n,"/",model.n,"/",sep="")
    projname <- paste(model.n,time.n,sep="")
    Listf <- list.files(pathf,full.names=T)
    Listf <- Listf[-12]
    Listf <- Listf[-11]
    environ <- lapply(Listf,raster)
    #subclim <- get(paste("subclim",n,sep=""))
    #cat("\n",subclim)
    myExplf <- stack(environ[subclim])
    BiomodEF <- BIOMOD_EnsembleForecasting(
      EM.output = myBiomodEM,
      new.env = myExplf,
      proj.name = projname,
      binary.meth = 'ROC',
      total.consensus = TRUE)
  }
}

#### 5. Plot, calculate distribution change ####
#_ 5.1 Plot current distribution#####
#spnames <- c("Os","Os2","Os3","Os4")
modelname <- c("bc","he","ip","mg","no","average")
futuretime <- c("5026","5085","7026","7085")

# Binary habitat suitability distribution
png(filename="baseline.jpg",width=1200,height=900,bg="white")
rbin <- raster(paste(sp.n,"/proj_current/proj_current_",
                     sp.n,"_ensemble_ROCbin.grd",
                     sep=""),band=1)
plot(rbin,box=F,axes=F,legend.width=1,main=n)
myRespXY <- get(paste("myRespXY",n,sep=""))
points(myRespXY$X, myRespXY$Y,col="red")
dev.off()

# Habitat suitability distribution
png(filename="baselinep.jpg",width=1200,height=900,bg="white")
rbin <- raster(paste(sp.n,"/proj_current/proj_current_",
                       sp.n,"_ensemble.grd",
                       sep=""),band=1)
plot(rbin,box=F,axes=F,legend.width=1,main=n)
myRespXY <- get(paste("myRespXY",n,sep=""))
points(myRespXY$X, myRespXY$Y,col="red")
dev.off()

# Manual transform the habitat suitability
rbin <- raster(paste(sp.n,"/proj_current/proj_current_",
                     sp.n,"_ensemble.grd",
                     sep=""),band=1)
png(filename="baselinep2.jpg",width=1200,height=900,bg="white")
par(mfrow=c(1,1))
rbin1 <- BinaryTransformation(rbin, 750)
plot(rbin1)
bou_county <- shapefile("county/county_wgs/county_boundary.shp")
plot(bou_county,add=T)
points(myRespXY$X, myRespXY$Y,col="red")
dev.off()

#_ 5.2 Calculate area change#######
#__calculate using binary map##
modelname <- c("bc","he","ip","mg","no","average")
futuretime <- c("5026","5085","7026","7085")
#n=1
changearea <- matrix(nrow=24,ncol=10)
rownamechange <- rep("NA",times=24)
spdisc <- raster(paste(sp.n,"/proj_current/proj_current_", # current distribution
                          sp.n,"_ensemble_ROCbin.grd",
                          sep=""),band=1)
for (time.n in futuretime){
  for (model.n in modelname){
    fscenario <- paste("/proj_",model.n,time.n,sep="")
    cat('/n',fscenario)
    spdisf <- raster(paste(sp.n,fscenario,fscenario,"_",
                           sp.n,"_ensemble_ROCbin.grd",
                           #sp.n,"_ensemble.grd",
                           sep=""),band=1)
    #spdisf <- BinaryTransformation(spdisf, 750)
    changecf <- BIOMOD_RangeSize(spdisc,spdisf,SpChange.Save=paste("change",sp.n,model.n,time.n,sep=""))
    changearea[n,] <- changecf$Compt.By.Models
    rownamechange[n] <- paste(sp.n,time.n,model.n,sep="")
    writeRaster(changecf$Diff.By.Pixel[[1]],paste(sp.n,"/binary/change",model.n,time.n,".asc",sep=""),format="ascii")
    n=n+1
  }
}

rownames(changearea) <- rownamechange
colnames(changearea) <- colnames(changecf$Compt.By.Models)
write.csv(changearea,file=paste(sp.n,"changearea.txt",sep=""))

#Synthesize results, convert to gain, stable, loss
synmodelchange <- list()
for (i in 1:4){
  fgain <- list()
  floss <- list()
  fstable1 <- list()
  fstable0 <- list()
  for (n in 1:6){
    fraster <- raster(paste(sp.n,"/binary/change",modelname[n],futuretime[i],".asc",sep=""))
    cat('\n',futuretime[i],modelname[n])
    cgain <- c(-3,0.5,0, 0.5,1.5,1 )
    closs <- c(-3,-1.5,1, -1.5,1.5,0)
    cstable1 <- c(-3,-1.5,0, -1.5,-0.5,1, -0.5,1.5,0 )
    cstable0 <- c(-3,-0.5,0, -0.5,0.5,1, 0.5,1.5,0 )
    cgain <- matrix(cgain, ncol=3, byrow=TRUE)
    closs <- matrix(closs, ncol=3, byrow=TRUE)
    cstable1 <- matrix(cstable1, ncol=3, byrow=TRUE)
    cstable0 <- matrix(cstable0, ncol=3, byrow=TRUE)
    fgain[[n]] <- reclassify(fraster, cgain)
    floss[[n]] <- reclassify(fraster, closs)
    fstable1[[n]] <- reclassify(fraster, cstable1)
    fstable0[[n]] <- reclassify(fraster, cstable0)
  }
  fgain <- sum(stack(fgain))
  floss <- sum(stack(floss))
  cgains <- c(-1,3.5,0, 3.5,7,1)
  closss <- c(-1,3.5,0, 3.5,7,-2)
  cgains <- matrix(cgains, ncol=3, byrow=TRUE)
  closss <- matrix(closss, ncol=3, byrow=TRUE)
  fstable1 <- sum(stack(fstable1))
  fstable0 <- sum(stack(fstable0))
  cstable1 <- c(-1,3.5,0, 3.5,7,-1)
  cstable0 <- c(-1,3.5,0, 3.5,7,0)
  cstable1 <- matrix(cstable1, ncol=3, byrow=TRUE)
  cstable0 <- matrix(cstable0, ncol=3, byrow=TRUE)
  rgain <- reclassify(fgain, cgains)
  rloss <- reclassify(floss, closss)
  rstable1 <- reclassify(fstable1, cstable1)
  rstable0 <- reclassify(fstable0, cstable0)
  sumraster <- rgain+rloss+rstable1+rstable0
  synmodelchange[[i]] <- freq(sumraster)
  writeRaster(sumraster,file=paste(sp.n,"/plot/change",futuretime[i],".asc",sep=""),format="ascii")
}

#_ 5.3 Calculate elevation change####
# Current elevation distribution 
ele <- crop(ele,extent(binp))
rbin1 <- raster(paste(sp.n,"/proj_current/proj_current_",
                      sp.n,"_ensemble_ROCbin.grd",
                      sep=""),band=1)
v <- c(-Inf,0.5,0,0.5,1.5,10000)#0=not exist, 1=distribution area margin,2=core distribution area
v <- matrix(v,ncol=3,byrow=TRUE)
binpc1 <- reclassify(rbin1,v)
eled <- binpc1+ele
elef <- freq(eled)
write.csv(elef,file=paste(sp.n,"/plot/current",".txt",sep="")) 

# Future elevation change
v <- c(-Inf,-1.5,10000, -1.5,-0.5,100000,-0.5,0.5,0, 0.5,1.5,1000000)#0=not exist, 1=distribution area margin,2=core distribution area
v <- matrix(v,ncol=3,byrow=TRUE)
sp.n <- "Os1"
for (i in 1:4){
  binf <- raster(paste(sp.n,"/plot/change",futuretime[i],".asc",sep=""))
  binpc <- reclassify(binf,v)
  eled <- binpc+ele
  elef <- freq(eled)
  write.csv(elef,file=paste(sp.n,"/plot/ele",futuretime[i],".txt",sep="")) 
}

# plot
for (i in 1:4){
  elef <- read.csv(paste(sp.n,"/plot/ele",futuretime[i],".txt",sep=""),stringsAsFactors=FALSE)
  elef <- elef[-length(elef),]
  elel <- list()
  elel[[1]]<- elef[which((elef$value>10000)&(elef$value<100000)),2:3]
  elel[[1]][,1] <- elel[[1]][,1]-10000
  elel[[2]]<- elef[which((elef$value>100000)&(elef$value<1000000)),2:3]
  elel[[2]][,1] <- elel[[2]][,1]-100000
  elel[[3]]<- elef[which((elef$value>1000000)),2:3]
  elel[[3]][,1] <- elel[[3]][,1]-1000000
  cutn <- c(0,3000,3500,4000,4500,5000,9000)
  eler <- matrix(nrow=6,ncol=6)
  for (n in 1:6){
    eler[n,1] <- sum(elel[[1]][which((elel[[1]][,1]>cutn[n])&(elel[[1]][,1]<cutn[n+1])),2])
    eler[n,2] <- round(eler[n,1]/489404,2)*100
    eler[n,3] <- sum(elel[[2]][which((elel[[2]][,1]>cutn[n])&(elel[[2]][,1]<cutn[n+1])),2])
    eler[n,4] <- round(eler[n,3]/489404,2)*100
    eler[n,5] <- sum(elel[[3]][which((elel[[3]][,1]>cutn[n])&(elel[[3]][,1]<cutn[n+1])),2])
    eler[n,6] <- round(eler[n,5]/489404,2)*100
  }
  plotele <- t(cbind(eler[,2],eler[,4],eler[,6]))
  rownames(plotele) <- c("Predicted loss","Predicted stable", "Predicted gain")
  colnames(plotele) <- c("<3000","3000-3500","3500-4000","4000-4500","4500-5000",">5000")
  #plotele2 <- melt(plotele1,id.var=row.names)
  #colnames(plotele2) <- c("Elevation","State","Value")  
  png(filename=paste(sp.n,"/plot/ele_",sp.n,"_1",futuretime[i],".jpg",sep=""),width=1200,height=900,bg="white")
  par(mfrow = c(1,1), family = 'san serif',cex.axis = 3.5, cex.lab=4,mar = c(20,11,5,3))#,mar=c(3,3,3,3)
  
  #par(mfrow=c(1,1))
  #barplot(plotele,col=c("#FF851B","gold","yellowgreen"),border=NA,ylim=c(0,45),xlab="Elevation(m)",width=3,mgp=c(17.5,2,0),las=2,tcl=-0.8)
  barplot(plotele,col=c("#FF851B","gold","yellowgreen"),border=NA,ylim=c(0,45),xaxt="n",width=0.9)#,mgp=c(17.5,2,0),las=2,tcl=-0.8)

  text(1:6, par("usr")[3]-3, srt = 45, adj = 1,
       labels = colnames(plotele), xpd = TRUE,cex=3.5) 
  title(xlab="Elevation(m)",mgp=c(15.5,2,0))
  title(ylab="Percentage (%)",mgp=c(7,3.5,0)) 
  
  #ggplot(plotele2,aes(x=Elevation, Y=Value,fill=State))+geom_bar(position="stack")
  #title(ylab="Percentage(%)",mgp=c(7,3.5,0))
  #legend(x=14.1,y=45 ,legend=rownames(plotele),fill=c("#FF851B","gold","yellowgreen"),bty="n",x.intersp =0.7, y.intersp = 0.9,cex=3.5)
  dev.off()
  write.csv(plotele,file=paste(sp.n,"/plot/sumele",futuretime[i],".txt",sep=""))
}

plotelec(sumele5085,"5085")
plotelec(sumele7026,"7026")
plotelec(sumele7085,"7085")

# plot 5026 separately
family = 'serif'
sumele5026 <- sumele5026[,-1]
sumele5026[4,] <- sumele5026[2,]
sumele5026[5,] <- sumele5026[1,]
sumele5026[6,] <- sumele5026[3,]
sumele5026 <- sumele5026[-1,]
sumele5026 <- sumele5026[-1,]
sumele5026 <- sumele5026[-1,]
rownames(sumele5026) <- c("3", "2","1")
colnames(sumele5026) <- c("<3000","3000-3500","3500-4000","4000-4500","4500-5000",">5000")
sumele5026 <- as.matrix(sumele5026)

png(filename=paste("os1_plot/ele_os1_50261.jpg",sep=""),width=1200,height=900,bg="white")
par(mfrow = c(1,1), cex.axis = 3.5, cex.lab=4,mar = c(20,11,5,3))#,mar=c(3,3,3,3)

barplot(sumele5026,col=c("grey22","grey48","grey69"),border=NA,ylim=c(0,45),las=2,xaxt="n",width=0.9)#,mgp=c(17.5,2,0),las=2,tcl=-0.8)

text(1:6, par("usr")[3]-3, srt = 45, adj = 1,
     labels = colnames(plotele), xpd = TRUE,cex=3.5) 
title(xlab="Elevation(m)",mgp=c(15.5,2,0))
title(ylab="Percentage (%)",mgp=c(7,3.5,0)) 
legend(x=4.2,y=45 ,legend=c("Predicted gain","Predicted loss","Predicted stable"),fill=c("grey69","grey48","grey22"),border=NA,bty="n",x.intersp =0.7, y.intersp = 0.9,cex=3.5)
dev.off()

for (n in 6:20){
  listname <- paste("List",n,sep="")
  lapply(get(listname),cropAndwrite,extent(ra)) #string to variable
}

##_ 5.4 Calculate the change in temperature and precipitation under different elevation##########
ele <- raster("crop1/alt.asc")
v <- c(-Inf,1000,NA,1000,3000,0,3000,3500,10000,3500,4000,20000,4000,4500,30000,4500,5000,40000,5000,Inf,50000)#0=not exist, 1=distribution area margin,2=core distribution area
v <- matrix(v,ncol=3,byrow=TRUE)
elec <- reclassify(ele,v)
sp.n <- "Os1"
cnames <- c("02","03","08","15","16","17","18")
modelname <- c("bc","he","ip","mg","no","average")
futuretime <- c("5026","5085","7026","7085")

# factor.n = "08"
# i=1
p=6

for (factor.n in cnames){
  temp_matrix <- matrix(nrow=144,ncol=6)
  for (i in 1:4){
    #for (p in 1:6){
    envl <- list()
    envd <- raster(paste("future climate/climatediffer/",futuretime[i],modelname[p],factor.n,".asc",sep=""))
    envele <- envd+elec
    envf <- freq(envele)
    envf <- data.frame(envf[-nrow(envf),])
    envm <- max(envf$value)
    cat("\n",envm)
    sumcounta <- sum(envf$count)
    envl[[1]]<- envf[which((envf$value>-5000)&(envf$value<5000)),1:2]
    envl[[2]]<- envf[which((envf$value>5000)&(envf$value<15000)),1:2]
    envl[[2]][,1] <- envl[[2]][,1]-10000
    
    envl[[3]]<- envf[which((envf$value>15000)&(envf$value<25000)),1:2]
    envl[[3]][,1] <- envl[[3]][,1]-20000
    envl[[4]]<- envf[which((envf$value>25000)&(envf$value<35000)),1:2]
    envl[[4]][,1] <- envl[[4]][,1]-30000
    envl[[5]]<- envf[which((envf$value>35000)&(envf$value<45000)),1:2]
    envl[[5]][,1] <- envl[[5]][,1]-40000
    envl[[6]]<-envf[which((envf$value>45000)),1:2]
    envl[[6]][,1] <- envl[[6]][,1]-50000
  }
  rownames(temp_matrix) <- rep(c("e3000","e3500","e4000","e4500","e5000","e6000"),times=24)
  colnames(temp_matrix) <- c("max","min","mean","mode","percent","percentall")
  capture.output(temp_matrix,
                 file=file.path(paste("future climate/",factor.n,"ele.txt", sep="")))
}

cabin1 <- raster("Os1/proj_current/proj_current_Os1_ensemble.grd",band=1)
cabin2 <- raster("Os1/proj_current/proj_current_Os1_ensemble_ROCbin.grd",band=1)
cabin3 <- raster("Os1/proj_current/proj_current_Os1_ensemble_ROCbin.grd",band=2)

