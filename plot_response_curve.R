# Plot response curve of SDM

library(biomod2)
library(maptools)
library(raster)
library(rgdal)
library(sp)
library(spatial)
library(ggplot2)

######## 0.Preparation ######
envList <- list.files("project/",full.names=T)
envNames <- list.files("project/",full.names=F)
environ <- lapply(envList,raster,crs=proj)
alt <- raster("crop1/alt.asc")
bio18 <- raster("project/bio18.asc")
bio08 <- raster("project/bio08.asc")
bio16 <- raster("project/bio16.asc")
herb <- raster("project/consensus_full_class_6.asc")
extent(alt)

myRespXY <- read.csv("Occ.csv",stringsAsFactors = F, header = T)
colnames(myRespXY) <- c("Long","Lat")
myResp <- as.numeric(rep(1,nrow(myRespXY)))

myRespName <- "Os1"
n_absences <- 10*sum(myResp==1,na.rm=TRUE)

########1. Generate data for response curve #######
generate_model <- function(ExplName){
  myRespName <- paste("Os_",ExplName,sep="")
  ModelId <- paste("resp_",ExplName,sep="")
  Expl <- get(ExplName)
  BioData <- BIOMOD_FormatingData(resp.var = myResp,
                                      expl.var = Expl,
                                      resp.xy = myRespXY,
                                      resp.name = myRespName,
                                      PA.nb.rep = 2,
                                      PA.nb.absences = n_absences,      # random get the pseudo-absence data
                                      PA.strategy = 'random',
                                      na.rm = T)
  ModelOut <- BIOMOD_Modeling(BioData,
                                     models = c('GAM','MARS','MAXENT','RF','CTA','GLM'),
                                     #models.options = myBiomodOption,
                                     NbRunEval=5,
                                     DataSplit=75,
                                     Yweights=NULL,
                                     VarImport=3,
                                     models.eval.meth = c('TSS','ROC'),
                                     SaveObj = TRUE,
                                     rescal.all.models = FALSE,
                                     modeling.id=ModelId)
  return (ModelOut)
}

Out_alt <- generate_model("Os_alt",alt,"resp_alt")
Out_bio08 <- generate_model("bio08")
Out_bio16 <- generate_model("bio16")
Out_bio18 <- generate_model("bio18")
Out_herb <- generate_model("herb")

BioData_alt <- BIOMOD_FormatingData(resp.var = myResp,
                                   expl.var = alt,
                                   resp.xy = myRespXY,
                                   resp.name = Os,
                                   PA.nb.rep = 2,
                                   PA.nb.absences = 1880,      # random get the pseudo-absence data
                                   PA.strategy = 'random',
                                   na.rm = T)

i = 19
BioData_c19 <- BIOMOD_FormatingData(resp.var = myRespn,
                                   expl.var = cropr[[i]],
                                   resp.xy = myRespXYn,
                                   resp.name = paste("Os_c",i,sep=""),
                                   PA.nb.rep = 2,
                                   PA.nb.absences = 1880,      # random get the pseudo-absence data
                                   PA.strategy = 'random',
                                   na.rm = T)


i = 9
BioData_c9 <- BIOMOD_FormatingData(resp.var = myRespn,
                                   expl.var = cropr[[i]],
                                   resp.xy = myRespXYn,
                                   resp.name = paste("Os_c",i,sep=""),
                                   PA.nb.rep = 2,
                                   PA.nb.absences = 1880,      # random get the pseudo-absence data
                                   PA.strategy = 'random',
                                   na.rm = T)

i <- 25
BioData_c25 <- BIOMOD_FormatingData(resp.var = myRespn,
                                    expl.var = cropr[[i]],
                                    resp.xy = myRespXYn,
                                    resp.name = paste("Os_c",i,sep=""),
                                    PA.nb.rep = 2,
                                    PA.nb.absences = 1880,      # random get the pseudo-absence data
                                    PA.strategy = 'random',
                                    na.rm = T)

i <- 17
BioData_c17 <- BIOMOD_FormatingData(resp.var = myRespn,
                                    expl.var = cropr[[i]],
                                    resp.xy = myRespXYn,
                                    resp.name = paste("Os_c",i,sep=""),
                                    PA.nb.rep = 2,
                                    PA.nb.absences = 1880,      # random get the pseudo-absence data
                                    PA.strategy = 'random',
                                    na.rm = T)

########## 2. Running Model ####################
Out_c1 <- BIOMOD_Modeling(BioData_c1, 
                           models = c('GAM','MARS','MAXENT','RF','GLM'),
                           NbRunEval = 5,
                           DataSplit = 75,      # 75% training, 25% testing
                           Yweights = NULL,
                           VarImport = 3,       # 
                           models.eval.meth = c('TSS','ROC'),
                           SaveObj = TRUE,
                           rescal.all.models = F,
                           do.full.models = F)

Out_c9 <- BIOMOD_Modeling(BioData_c9, 
                          models = c('GAM','MARS','MAXENT','RF','GLM'),
                          NbRunEval = 5,
                          DataSplit = 75,      # 75% training, 25% testing
                          Yweights = NULL,
                          VarImport = 3,       # 
                          models.eval.meth = c('TSS','ROC'),
                          SaveObj = TRUE,
                          rescal.all.models = F,
                          do.full.models = F)

Out_c19 <- BIOMOD_Modeling(BioData_c19, 
                           models = c('GAM','MARS','MAXENT','RF','GLM'),
                           NbRunEval = 5,
                           DataSplit = 75,      # 75% training, 25% testing
                           Yweights = NULL,
                           VarImport = 3,       # 
                           models.eval.meth = c('TSS','ROC'),
                           SaveObj = TRUE,
                           rescal.all.models = F,
                           do.full.models = F)

Out_c25 <- BIOMOD_Modeling(BioData_c25, 
                           models = c('GAM','MARS','MAXENT','RF','GLM'),
                           NbRunEval = 5,
                           DataSplit = 75,      # 75% training, 25% testing
                           Yweights = NULL,
                           VarImport = 3,       # 
                           models.eval.meth = c('TSS','ROC'),
                           SaveObj = TRUE,
                           rescal.all.models = F,
                           do.full.models = F)

Out_c17 <- BIOMOD_Modeling(BioData_c17, 
                           models = c('GAM','MARS','MAXENT','RF','GLM'),
                           NbRunEval = 5,
                           DataSplit = 75,      # 75% training, 25% testing
                           Yweights = NULL,
                           VarImport = 3,       # 
                           models.eval.meth = c('TSS','ROC'),
                           SaveObj = TRUE,
                           rescal.all.models = F,
                           do.full.models = F)


Model_c1 <- BIOMOD_LoadModels(Out_c1)

Model_c9 <- BIOMOD_LoadModels(Out_c9)

Model_c19 <- BIOMOD_LoadModels(Out_c19)

Model_c25 <- BIOMOD_LoadModels(Out_c25)

Model_c17 <- BIOMOD_LoadModels(Out_c17)
rm(list=Model_c17)
rm(list=Model_c9)
rm(list=Model_c1)
rm(list=Model_c19)
rm(list=Model_c25)

####### 3. Plot ############
#######C1 ########


ModelOut = "Out_bio18"
xlabn = "Precipitation of Warmest Quarter"
xrange = c(0,1800)


get_respdata <- function(ModelOut,xlabn,xrange){
  ModelOut1 <- get(ModelOut)
  loadmodels <- BIOMOD_LoadModels(ModelOut1)
  resplot <- response.plot2(models = loadmodels,
                            Data = get_formal_data(ModelOut1,'expl.var'),
                            show.variables = ModelOut1@expl.var.names,
                            do.bivariate = F,
                            plot = F,
                            fixed.var.metric = 'mean',
                            legend = F)
  #rm(list=loadmodels)
  resplot <- data.frame(resplot[[1]])
  nameres <- colnames(resplot)
  res_mean <- matrix(nrow=100,ncol=5) 
  modelname <- c("GAM","MAXENT","MARS","GLM")
  res_mean[,1] <- resplot[,1]
  
  png(filename=paste(ModelOut,".jpg",sep=""),width=1200,height=900,bg="white")
  par(mfrow = c(1,1),mar = c(16,15,5,3),cex.axis = 4.5, cex.lab=5)
  plot(NULL,NULL,type="n",ylim = c(0,1),xlim = xrange, mgp=c(9.5,3.5,0),tcl=-0.8, las =1,xlab=xlabn, ylab = "")
  title(ylab="Probability of Occurrence",mgp=c(11,3.5,0))
  for (n in 1:4){
    num_n <- nameres[grep(modelname[n],nameres)]
    res_n <- resplot[,num_n]
    for (i in 1:12){
      lines(res_mean[,1],res_n[,i],col = "#E0EEEE")
    }
    res_mean[,n+1] <- apply(res_n,1,mean)
  }
  lines(res_mean[,1],res_mean[,2],type="l",
        col = "#FF4136", lwd = 2.5)
  lines(res_mean[,1],res_mean[,3],type="l",
        col = "#3D9970", lwd = 2.5)
  lines(res_mean[,1],res_mean[,4],type="l",
        col = "#FF851B", lwd = 2.5)
  lines(res_mean[,1],res_mean[,5],type="l",
        col = "#7AC5CD", lwd = 2.5)
  
  legend(x=17,y=1.05,c("GAM","GLM", "MARS", "MaxEnt"),cex=4,  #bio18
        col = c("#FF4136","#7AC5CD", "#FF851B", "#3D9970"),
        x.intersp =0.7, y.intersp = 0.9,bty = "n",lty = 1,lwd=4,seg.len=1.1)
  dev.off()
  colnames(res_mean) <- c("prob",modelname)
  return (res_mean)
}

Model_alt <- BIOMOD_LoadModels(Out_alt)
Model_bio08 <- BIOMOD_LoadModels(Out_bio08)
Model_bio16 <- BIOMOD_LoadModels(Out_bio16)
Model_bio18 <- BIOMOD_LoadModels(Out_bio18)
Model_herb <- BIOMOD_LoadModels(Out_herb)

respdata_alt <- get_respdata("Out_alt")
respdata_bio08 <- get_respdata("Out_bio08",xlabn = "MTWQ (°C)", c(-3,26))
respdata_bio16 <- get_respdata("Out_bio16",xlabn = "PWEQ (mm)", c(0,1000))
respdata_bio18 <- get_respdata("Out_bio18",xlabn = "PWAQ (mm)", c(0,1000))
respdata_herb <- get_respdata("Out_herb",xlabn ="HV", c(0,100))

rm(list=Model_alt)
rm(list=Model_bio08)
rm(list=Model_bio16)
rm(list=Model_bio18)
rm(list=Model_herb)


oldpar <- par()
par(mfrow = c(1,2),mar = c(7,7,5,2), family = 'serif',cex.axis = 2, cex.lab=2)

plot(respdata_alt[,7],respdata_alt[,1],type="l",ylim = c(0,1),
     col = "#E0EEEE", lwd = 2,
     xlab="Elevation", ylab = "Probability of occurrence")

#xlim = c(0,30),
for (i in 1:10){
  lines(resplot_c1[,1],res_GAM[,i],col = "#E0EEEE")
}
for (i in 1:10){
  lines(resplot_c1[,1],res_MAXENT[,i],col = "#E0EEEE")
}
#for (i in 1:10){
#  lines(resplot_c1[,1],res_CTA[,i],col = "#E0EEEE")
#}
for (i in 1:10){
  lines(resplot_c1[,1],res_MARS[,i],col = "#E0EEEE")
}
for (i in 1:10){
  lines(resplot_c1[,1],res_GLM[,i],col = "#E0EEEE")
}
lines(resplot_c1[,1],res_mean_c1[,1],type="l",
      col = "#FF4136", lwd = 2)
lines(resplot_c1[,1],res_mean_c1[,2],type="l",
      col = "#3D9970", lwd = 2)
lines(resplot_c1[,1],res_mean_c1[,3],type="l",
      col = "#FF851B", lwd = 2)
lines(resplot_c1[,1],res_mean_c1[,4],type="l",
      col = "#7AC5CD", lwd = 2)
#lines(resplot_c1[,1],res_mean_c1[,5],type="l",
#      col = "#7AC5CD", lwd = 2)




resplot_c1 <- response.plot2(models = Model_c1,
                             Data = get_formal_data(Out_c1,'expl.var'),
                             show.variables = Out_c1@expl.var.names,
                             do.bivariate = F,
                             plot = F,
                             fixed.var.metric = 'mean',
                             legend = F)

resplot_c1 <- data.frame(resplot_c1[[1]])
nameres <- colnames(resplot_c1)
res_mean_c1 <- data.frame(rep(1,times=100)) 
n_GAM <- nameres[grep("GAM",nameres)]
res_GAM <- resplot_c1[,n_GAM]
#MAXENT
n_MAXENT <- nameres[grep("MAXENT",nameres)]
res_MAXENT <- resplot_c1[,n_MAXENT]
#CTA
#n_CTA <- nameres[grep("CTA",nameres)]
#res_CTA <- resplot_c1[,n_CTA]
#MARS
n_MARS <- nameres[grep("MARS",nameres)]
res_MARS <- resplot_c1[,n_MARS]
#GLM
n_GLM <- nameres[grep("GLM",nameres)]
res_GLM <- resplot_c1[,n_GLM]

res_mean_c1[,1] <- apply(res_GAM,1,mean)
res_mean_c1[,2] <- apply(res_MAXENT,1,mean)
#res_mean_c1[,3] <- apply(res_CTA,1,mean)
res_mean_c1[,3] <- apply(res_MARS,1,mean)
res_mean_c1[,4] <- apply(res_GLM,1,mean)

colnames(res_mean_c1) <- c("GAM", "MAXENT", "MARS","GLM")

oldpar <- par()
par(mfrow = c(1,2),mar = c(7,7,5,2), family = 'serif',cex.axis = 2, cex.lab=2)

plot(resplot_c1[,1],res_mean_c1[,1],type="l",ylim = c(0,1),
     col = "#E0EEEE", lwd = 2,
     xlab="Elevation", ylab = "Probability of occurrence")

#xlim = c(0,30),
for (i in 1:10){
  lines(resplot_c1[,1],res_GAM[,i],col = "#E0EEEE")
}
for (i in 1:10){
  lines(resplot_c1[,1],res_MAXENT[,i],col = "#E0EEEE")
}
#for (i in 1:10){
#  lines(resplot_c1[,1],res_CTA[,i],col = "#E0EEEE")
#}
for (i in 1:10){
  lines(resplot_c1[,1],res_MARS[,i],col = "#E0EEEE")
}
for (i in 1:10){
  lines(resplot_c1[,1],res_GLM[,i],col = "#E0EEEE")
}
lines(resplot_c1[,1],res_mean_c1[,1],type="l",
      col = "#FF4136", lwd = 2)
lines(resplot_c1[,1],res_mean_c1[,2],type="l",
      col = "#3D9970", lwd = 2)
lines(resplot_c1[,1],res_mean_c1[,3],type="l",
      col = "#FF851B", lwd = 2)
lines(resplot_c1[,1],res_mean_c1[,4],type="l",
      col = "#7AC5CD", lwd = 2)
#lines(resplot_c1[,1],res_mean_c1[,5],type="l",
#      col = "#7AC5CD", lwd = 2)



###########c17##########
resplot_c17 <- response.plot2(models = Model_c17,
                             Data = get_formal_data(Out_c17,'expl.var'),
                             show.variables = Out_c17@expl.var.names,
                             do.bivariate = F,
                             plot = F,
                             fixed.var.metric = 'mean',
                             legend = F)

resplot_c17 <- data.frame(resplot_c17[[1]])
nameres <- colnames(resplot_c17)
res_mean_c17 <- data.frame(rep(1,times=100)) 
n_GAM <- nameres[grep("GAM",nameres)]
res_GAM <- resplot_c17[,n_GAM]
#MAXENT
n_MAXENT <- nameres[grep("MAXENT",nameres)]
res_MAXENT <- resplot_c17[,n_MAXENT]
#CTA
#n_CTA <- nameres[grep("CTA",nameres)]
#res_CTA <- resplot_c17[,n_CTA]
#MARS
n_MARS <- nameres[grep("MARS",nameres)]
res_MARS <- resplot_c17[,n_MARS]
#GLM
n_GLM <- nameres[grep("GLM",nameres)]
res_GLM <- resplot_c17[,n_GLM]

res_mean_c17[,1] <- apply(res_GAM,1,mean)
res_mean_c17[,2] <- apply(res_MAXENT,1,mean)
#res_mean_c17[,3] <- apply(res_CTA,1,mean)
res_mean_c17[,3] <- apply(res_MARS,1,mean)
res_mean_c17[,4] <- apply(res_GLM,1,mean)

colnames(res_mean_c17) <- c("GAM", "MAXENT", "MARS","GLM")

plot(resplot_c17[,1],res_mean_c17[,1],type="l",ylim = c(0,1),
     xlim = c(0,1700),col = "#E0EEEE", lwd = 2,
     xlab="Precipitation of Wettest Quarter", ylab = "Probability of Occurrence")

#
for (i in 1:10){
  lines(resplot_c17[,1],res_GAM[,i],col = "#E0EEEE")
}
for (i in 1:10){
  lines(resplot_c17[,1],res_MAXENT[,i],col = "#E0EEEE")
}
#for (i in 1:10){
#  lines(resplot_c17[,1],res_CTA[,i],col = "#E0EEEE")
#}
for (i in 1:10){
  lines(resplot_c17[,1],res_MARS[,i],col = "#E0EEEE")
}
for (i in 1:10){
  lines(resplot_c17[,1],res_GLM[,i],col = "#E0EEEE")
}
lines(resplot_c17[,1],res_mean_c17[,1],type="l",
      col = "#FF4136", lwd = 2)
lines(resplot_c17[,1],res_mean_c17[,2],type="l",
      col = "#3D9970", lwd = 2)
lines(resplot_c17[,1],res_mean_c17[,3],type="l",
      col = "#FF851B", lwd = 2)
lines(resplot_c17[,1],res_mean_c17[,4],type="l",
      col = "#7AC5CD", lwd = 2)
#lines(resplot_c17[,1],res_mean_c17[,5],type="l",
#      col = "#7AC5CD", lwd = 2)


legend(x=1200,y=1.05,c("GAM", "MAXENT", "MARS","GLM"),
       col = c("#FF4136", "#3D9970", "#FF851B","#7AC5CD"),
       x.intersp =0.3, y.intersp = 0.3,bty = "n",lty = 1,lwd=2,seg.len=0.8)

###########c9##########
resplot_c9 <- response.plot2(models = Model_c9,
                              Data = get_formal_data(Out_c9,'expl.var'),
                              show.variables = Out_c9@expl.var.names,
                              do.bivariate = F,
                              plot = F,
                              fixed.var.metric = 'mean',
                              legend = F)
resplot_c9 <- data.frame(resplot_c9[[1]])
nameres <- colnames(resplot_c9)
res_mean_c9 <- data.frame(rep(1,times=100)) 
n_GAM <- nameres[grep("GAM",nameres)]
res_GAM <- resplot_c9[,n_GAM]
#MAXENT
n_MAXENT <- nameres[grep("MAXENT",nameres)]
res_MAXENT <- resplot_c9[,n_MAXENT]
#CTA
#n_CTA <- nameres[grep("CTA",nameres)]
#res_CTA <- resplot_c9[,n_CTA]
#MARS
n_MARS <- nameres[grep("MARS",nameres)]
res_MARS <- resplot_c9[,n_MARS]
#GLM
n_GLM <- nameres[grep("GLM",nameres)]
res_GLM <- resplot_c9[,n_GLM]

res_mean_c9[,1] <- apply(res_GAM,1,mean)
res_mean_c9[,2] <- apply(res_MAXENT,1,mean)
#res_mean_c9[,3] <- apply(res_CTA,1,mean)
res_mean_c9[,3] <- apply(res_MARS,1,mean)
res_mean_c9[,4] <- apply(res_GLM,1,mean)

colnames(res_mean_c9) <- c("GAM", "MAXENT", "MARS","GLM")

plot(resplot_c9[,1],res_mean_c9[,1],type="l",ylim = c(0,1),
     #xlim = c(0,2000),
     col = "#E0EEEE", lwd = 2,
     xlab="Mean Temperature of Wettest Quarter", ylab = "Probability of Occurrence")

#
for (i in 1:10){
  lines(resplot_c9[,1],res_GAM[,i],col = "#E0EEEE")
}
for (i in 1:10){
  lines(resplot_c9[,1],res_MAXENT[,i],col = "#E0EEEE")
}
#for (i in 1:10){
#  lines(resplot_c9[,1],res_CTA[,i],col = "#E0EEEE")
#}
for (i in 1:10){
  lines(resplot_c9[,1],res_MARS[,i],col = "#E0EEEE")
}
for (i in 1:10){
  lines(resplot_c9[,1],res_GLM[,i],col = "#E0EEEE")
}
lines(resplot_c9[,1],res_mean_c9[,1],type="l",
      col = "#FF4136", lwd = 2)
lines(resplot_c9[,1],res_mean_c9[,2],type="l",
      col = "#3D9970", lwd = 2)
lines(resplot_c9[,1],res_mean_c9[,3],type="l",
      col = "#FF851B", lwd = 2)
lines(resplot_c9[,1],res_mean_c9[,4],type="l",
      col = "#7AC5CD", lwd = 2)
#lines(resplot_c9[,1],res_mean_c9[,5],type="l",
#      col = "#7AC5CD", lwd = 2)



###########c19##########
resplot_c19 <- response.plot2(models = Model_c19,
                              Data = get_formal_data(Out_c19,'expl.var'),
                              show.variables = Out_c19@expl.var.names,
                              do.bivariate = F,
                              plot = F,
                              fixed.var.metric = 'mean',
                              legend = F)

resplot_c19 <- data.frame(resplot_c19[[1]])
nameres <- colnames(resplot_c19)
res_mean_c19 <- data.frame(rep(1,times=100)) 
n_GAM <- nameres[grep("GAM",nameres)]
res_GAM <- resplot_c19[,n_GAM]
#MAXENT
n_MAXENT <- nameres[grep("MAXENT",nameres)]
res_MAXENT <- resplot_c19[,n_MAXENT]
#CTA
#n_CTA <- nameres[grep("CTA",nameres)]
#res_CTA <- resplot_c19[,n_CTA]
#MARS
n_MARS <- nameres[grep("MARS",nameres)]
res_MARS <- resplot_c19[,n_MARS]
#GLM
n_GLM <- nameres[grep("GLM",nameres)]
res_GLM <- resplot_c19[,n_GLM]

res_mean_c19[,1] <- apply(res_GAM,1,mean)
res_mean_c19[,2] <- apply(res_MAXENT,1,mean)
#res_mean_c19[,3] <- apply(res_CTA,1,mean)
res_mean_c19[,3] <- apply(res_MARS,1,mean)
res_mean_c19[,4] <- apply(res_GLM,1,mean)

colnames(res_mean_c19) <- c("GAM", "MAXENT", "MARS","GLM")

plot(resplot_c19[,1],res_mean_c19[,1],type="l",ylim = c(0,1),
     xlim = c(0,1700),
     col = "#E0EEEE", lwd = 2,
     xlab="Precipitation of Warmest Quarter", ylab = "Probability of Occurrence")

#
for (i in 1:10){
  lines(resplot_c19[,1],res_GAM[,i],col = "#E0EEEE")
}
for (i in 1:10){
  lines(resplot_c19[,1],res_MAXENT[,i],col = "#E0EEEE")
}
#for (i in 1:10){
#  lines(resplot_c19[,1],res_CTA[,i],col = "#E0EEEE")
#}
for (i in 1:10){
  lines(resplot_c19[,1],res_MARS[,i],col = "#E0EEEE")
}
for (i in 1:10){
  lines(resplot_c19[,1],res_GLM[,i],col = "#E0EEEE")
}
lines(resplot_c19[,1],res_mean_c19[,1],type="l",
      col = "#FF4136", lwd = 2)
lines(resplot_c19[,1],res_mean_c19[,2],type="l",
      col = "#3D9970", lwd = 2)
lines(resplot_c19[,1],res_mean_c19[,3],type="l",
      col = "#FF851B", lwd = 2)
lines(resplot_c19[,1],res_mean_c19[,4],type="l",
      col = "#7AC5CD", lwd = 2)
#lines(resplot_c19[,1],res_mean_c19[,5],type="l",
#      col = "#7AC5CD", lwd = 2)

legend(x=1200,y=1.05,c("GAM", "MAXENT", "CTA", "MARS","GLM"),
       col = c("#FF4136", "#3D9970", "#FF851B","#6B229B","#7AC5CD"),
       x.intersp =0.3, y.intersp = 0.3,bty = "n",lty = 1,lwd=2,seg.len=0.8)
###########c25##########
resplot_c25 <- response.plot2(models = Model_c25,
                              Data = get_formal_data(Out_c25,'expl.var'),
                              show.variables = Out_c25@expl.var.names,
                              do.bivariate = F,
                              plot = F,
                              fixed.var.metric = 'mean',
                              legend = F)

resplot_c25 <- data.frame(resplot_c25[[1]])
nameres <- colnames(resplot_c25)
res_mean_c25 <- data.frame(rep(1,times=100)) 
n_GAM <- nameres[grep("GAM",nameres)]
res_GAM <- resplot_c25[,n_GAM]
#MAXENT
n_MAXENT <- nameres[grep("MAXENT",nameres)]
res_MAXENT <- resplot_c25[,n_MAXENT]
#CTA
#n_CTA <- nameres[grep("CTA",nameres)]
#res_CTA <- resplot_c25[,n_CTA]
#MARS
n_MARS <- nameres[grep("MARS",nameres)]
res_MARS <- resplot_c25[,n_MARS]
#GLM
n_GLM <- nameres[grep("GLM",nameres)]
res_GLM <- resplot_c25[,n_GLM]

res_mean_c25[,1] <- apply(res_GAM,1,mean)
res_mean_c25[,2] <- apply(res_MAXENT,1,mean)
#res_mean_c25[,3] <- apply(res_CTA,1,mean)
res_mean_c25[,3] <- apply(res_MARS,1,mean)
res_mean_c25[,4] <- apply(res_GLM,1,mean)

colnames(res_mean_c25) <- c("GAM", "MAXENT", "MARS","GLM")

plot(resplot_c25[,1],res_mean_c25[,1],type="l",ylim = c(0,1),
     #xlim = c(0,2000),
     col = "#E0EEEE", lwd = 2,
     xlab="Herbaceous Vegetation", ylab = "Probability of Occurrence")

#
for (i in 1:10){
  lines(resplot_c25[,1],res_GAM[,i],col = "#E0EEEE")
}
for (i in 1:10){
  lines(resplot_c25[,1],res_MAXENT[,i],col = "#E0EEEE")
}
#for (i in 1:10){
#  lines(resplot_c25[,1],res_CTA[,i],col = "#E0EEEE")
#}
for (i in 1:10){
  lines(resplot_c25[,1],res_MARS[,i],col = "#E0EEEE")
}
for (i in 1:10){
  lines(resplot_c25[,1],res_GLM[,i],col = "#E0EEEE")
}
lines(resplot_c25[,1],res_mean_c25[,1],type="l",
      col = "#FF4136", lwd = 2)
lines(resplot_c25[,1],res_mean_c25[,2],type="l",
      col = "#3D9970", lwd = 2)
lines(resplot_c25[,1],res_mean_c25[,3],type="l",
      col = "#FF851B", lwd = 2)
lines(resplot_c25[,1],res_mean_c25[,4],type="l",
      col = "#7AC5CD", lwd = 2)
#lines(resplot_c25[,1],res_mean_c25[,5],type="l",
#      col = "#7AC5CD", lwd = 2)
