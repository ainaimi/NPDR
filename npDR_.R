# load libraries/functions
args=(commandArgs(TRUE))
setwd(".")

packages <- c("foreach","doParallel","doRNG","boot","rmutil","mvtnorm","gam","sandwich",
              "devtools","glmnet","data.table","rpart","ranger","nnet","arm","earth","e1071")
userLib <-  "~/R/R_LIBS_USER"
.libPaths(userLib)

ifelse(!dir.exists(userLib), dir.create(userLib), FALSE)

for (package in packages) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package,lib=userLib, repos='http://lib.stat.cmu.edu/R/CRAN')
  }
}

for (package in packages) {
  library(package, character.only=T, lib.loc=.libPaths())
}

devtools::install_github("ecpolley/SuperLearner",lib=userLib, repos='http://lib.stat.cmu.edu/R/CRAN')
library(SuperLearner)

install.packages("tmle",lib=userLib,repos='http://lib.stat.cmu.edu/R/CRAN')
library(tmle)

time<-proc.time()
print(args)
set.seed(as.numeric(args[[1]]))
expit <- function(x){ exp(x)/(1+exp(x)) }
logit <- function(x){ log(x/(1-x)) }
nsim <- as.numeric(args[[2]])

cols <- c("ipwPMT","ipwNPT","ipwPMF","ipwNPF",
          "regPMT","regNPT","regPMF","regNPF",
          "aipwPMT","aipwNPT","aipwPMF","aipwNPF",
          "tmlePMT","tmleNPT","tmlePMF","tmleNPF")
res.est <- data.frame(matrix(nrow=nsim*5,ncol=length(cols)));colnames(res.est) <- cols; res.se <- res.est

##  true value
true<-6;true
npDR<-function(counter,pNum,bs=T,bootNum=100){
  source("create.SL.glmnet.R")
  source("create.SL.gam.Wrapper.R")
  source("rangerWrapper.R")
  source("rangerWrapper1.R")
  source("rangerWrapper2.R")
  source("rangerWrapper3.R")
  source("rangerWrapper4.R")
  source("rangerWrapper5.R")
  source("create.SL.xgboost.R")
  
  tune = list(ntrees = c(100, 500), max_depth = c(1, 2), minobspernode = 10,
              shrinkage = c(0.1, 0.01, 0.001))
  xgb_grid = create.SL.xgboost(tune = tune) #, env = sl_env
  create.SL.glmnet(alpha = c(0,.5)) #setting alpha = 0 gives LASSO, 1 gives elastic net.
  create.SL.gam(deg.gam = c(3,4,5)) 
  sl.lib<-c("SL.glmnet.0","SL.glmnet.0.5","SL.glmnet","SL.rpartPrune","SL.svm",
            "SL.gam.3","SL.gam.4","SL.gam.5","SL.glm.interaction","SL.earth",
            "SL.ranger","SL.ranger1","SL.ranger2","SL.ranger3","SL.ranger4","SL.ranger5",
            "SL.bayesglm","SL.xgboost","SL.mean",xgb_grid$names)
  # data management
  i<-counter
  samp<-rep(1:5,nsim*5)[counter]
  ss<-c(50,100,200,600,1200)
  n<-ss[samp]
  p=pNum
  cat("Now running iteration",i,"with a sample size of",n,'\n');flush.console()
  
  # confounders
  sigma<-matrix(0,nrow=4,ncol=4);diag(sigma)<-1
  x <- rmvnorm(n, mean=rep(0,4), sigma=sigma)
  
  z<-x
  z[,1]<-exp(x[,1]/2)
  z[,2]<-x[,2]/(1+exp(x[,1]))+10
  z[,3]<-(x[,1]*x[,3]/25+.6)^3
  z[,4]<-(x[,2]*x[,4]+20)^2

  # pdf("eFigure1.pdf",width=6,height=6)
  # GGally::ggpairs(data.frame(z))
  # dev.off()
  
  # design matrix for outcome model

  muMatT<-model.matrix(as.formula(paste("~(",paste("x[,",1:ncol(x),"]",collapse="+"),")^2")))
  
  parms3<-c(3.5,2.5,-1,5,2,2.5,1.5,1.5,1.5,1) #,4.25,-2
  parms4<-c(log(2),log(2.5),log(.5),log(1.5),log(1.75),log(1.5),log(1.25),log(1.25),log(1.25),log(1.25)) #,log(2.25),log(.25)
  
  beta<-parms3;beta<-c(120,beta)
  # design matrix for propensity score model
  piMatT<-model.matrix(as.formula(paste("~(",paste("x[,",1:ncol(x),"]",collapse="+"),")^2")))
  theta<-parms4;theta<-c(-.5,theta)
  mu <- muMatT%*%beta
  # propensity score model
  pi <- expit(piMatT%*%theta);r<-1-rbinom(n,1,pi)
  # outcome model: true expsoure effect = 6
  y<- r*6 + mu + rnorm(n,0,20)
  
  # induce misspecification
  muMatF<-model.matrix(as.formula(paste("~(",paste("z[,",1:ncol(x),"]",collapse="+"),")")))
  piMatF<-model.matrix(as.formula(paste("~",paste("z[,",1:ncol(x),"]",collapse="+"))))
  
  # correct specification
  dat <- data.frame(x,r,y); colnames(dat)[1:ncol(x)] <- paste("x",1:ncol(x),sep="")
  W <- muMatT[,-1];colnames(W) <- paste("W",1:ncol(muMatT[,-1]), sep="");A <- r;Y <- y
  tQForm<-as.formula(paste0("Y~A+", paste(paste0("W",1:ncol(muMatT[,-1])), collapse="+")))
  tgForm<-as.formula(paste0("A~", paste(paste0("W",1:ncol(piMatT[,-1])), collapse="+")))
  tmlePMT <- tmle(Y,A,W,family="gaussian",Qform=tQForm,gform=tgForm)
  
  folds<-c(2,2,3,5,5)[samp]
  cat("Number of cross-validation folds is",folds,'\n');flush.console()
  tmleNPT <- tmle(Y,A,W,family="gaussian",Q.SL.library=sl.lib,g.SL.library=sl.lib)
  
  SL_out.monT<-data.table(N=n,model="outcome",t(tmleNPT$Qinit$coef))
  SL_exp.monT<-data.table(N=n,model="exposure",t(tmleNPT$g$coef))
  
  cbind(t(SL_exp.monT)[-c(1,2),],t(SL_out.monT)[-c(1,2),])
  
  if(i==1&samp==1){
    write.table(SL_out.monT,"SL_outDat_True.txt",sep="\t",row.names=F)
    write.table(SL_exp.monT,"SL_expDat_True.txt",sep="\t",row.names=F)
  } else{
    write.table(SL_out.monT,"SL_outDat_True.txt",sep="\t",row.names=F,col.names=F,append=T)
    write.table(SL_exp.monT,"SL_expDat_True.txt",sep="\t",row.names=F,col.names=F,append=T)
  }
  
  pihatPT <- tmlePMT$g$g1W
  swPT<-dat$r*(mean(dat$r)/tmlePMT$g$g1W) + (1-dat$r)*((1-mean(dat$r))/(1-tmlePMT$g$g1W))
  pihatNPT <- tmleNPT$g$g1W
  swNPT<-dat$r*(mean(dat$r)/tmleNPT$g$g1W) + (1-dat$r)*((1-mean(dat$r))/(1-tmleNPT$g$g1W))
  muhatPT  <- tmlePMT$Qinit$Q[,2]*A+tmlePMT$Qinit$Q[,1]*(1-A)
  muhatPT1 <- tmlePMT$Qinit$Q[,2];muhatPT0 <- tmlePMT$Qinit$Q[,1]
  muhatNPT  <- tmleNPT$Qinit$Q[,2]*A+tmleNPT$Qinit$Q[,1]*(1-A)
  muhatNPT1 <- tmleNPT$Qinit$Q[,2];muhatNPT0 <- tmleNPT$Qinit$Q[,1]
  
  # misspecified
  dat <- data.frame(z,r,y); colnames(dat)[1:ncol(z)] <- paste("z",1:ncol(x),sep="")
  W <- muMatF[,-1];colnames(W) <- paste("W",1:ncol(muMatF[,-1]), sep="");A <- r;Y <- y
  tQForm<-as.formula(paste0("Y~A+", paste(paste0("W",1:ncol(muMatF[,-1])), collapse="+")))
  tgForm<-as.formula(paste0("A~", paste(paste0("W",1:ncol(piMatF[,-1])), collapse="+")))
  tmlePMF <- tmle(Y,A,W,family="gaussian",Qform=tQForm,gform=tgForm)
  
  folds<-c(2,2,3,5,5)[samp]
  tmleNPF <- tmle(Y,A,W,family="gaussian",Q.SL.library=sl.lib,g.SL.library=sl.lib)
  
  SL_out.monF<-data.table(N=n,model="outcome",t(tmleNPF$Qinit$coef))
  SL_exp.monF<-data.table(N=n,model="exposure",t(tmleNPF$g$coef))
  
  cbind(t(SL_exp.monF)[-c(1,2),],t(SL_out.monF)[-c(1,2),])
  
  if(i==1&samp==1){
    write.table(SL_out.monF,"SL_outDat_MisSpec.txt",sep="\t",row.names=F)
    write.table(SL_exp.monT,"SL_expDat_MisSpec.txt",sep="\t",row.names=F)
  } else{
    write.table(SL_out.monT,"SL_outDat_MisSpec.txt",sep="\t",row.names=F,col.names=F,append=T)
    write.table(SL_exp.monT,"SL_expDat_MisSpec.txt",sep="\t",row.names=F,col.names=F,append=T)
  }
  
  pihatPF <- tmlePMF$g$g1W
  swPF<-dat$r*(mean(dat$r)/tmlePMF$g$g1W) + (1-dat$r)*((1-mean(dat$r))/(1-tmlePMF$g$g1W))
  pihatNPF <- tmleNPF$g$g1W
  swNPF<-dat$r*(mean(dat$r)/tmleNPF$g$g1W) + (1-dat$r)*((1-mean(dat$r))/(1-tmleNPF$g$g1W))
  muhatPF  <- tmlePMF$Qinit$Q[,2]*A+tmlePMF$Qinit$Q[,1]*(1-A)
  muhatPF1 <- tmlePMF$Qinit$Q[,2];muhatPF0 <- tmlePMF$Qinit$Q[,1]
  muhatNPF  <- tmleNPF$Qinit$Q[,2]*A+tmleNPF$Qinit$Q[,1]*(1-A)
  muhatNPF1 <- tmleNPF$Qinit$Q[,2];muhatNPF0 <- tmleNPF$Qinit$Q[,1]
  
  # propensityPT<-data.frame(pihatPT,dat$r);propensityPT$model<-1;names(propensityPT)[1]<-"pi"
  # head(propensityPT)
  # propensityPF<-data.frame(pihatPF,dat$r);propensityPF$model<-2;names(propensityPF)[1]<-"pi"
  # head(propensityPF)
  # propensityNT<-data.frame(pihatNPT,dat$r);propensityNT$model<-3;names(propensityNT)[1]<-"pi"
  # head(propensityNT)
  # propensityNF<-data.frame(pihatNPF,dat$r);propensityNF$model<-4;names(propensityNF)[1]<-"pi"
  # head(propensityNF)
  # propensity<-rbind(propensityPT,propensityPF,propensityNT,propensityNF)
  # propensity$model<-factor(propensity$model,labels=c("Parametric Corr","Parametric Miss","Nonparametric Corr","Nonparametric Miss"))
  # propensity$x<-factor(propensity$dat.r,labels=c("Unexposed","Exposed"));propensity$dat.r<-NULL
  # head(propensity)
  # 
  # pdf(file="~/Dropbox/Documents/Research/Papers/SemiparametricInference/eFigure2.pdf",height=8,width=6)
  # ggplot(propensity,aes(pi,color=x)) + 
  #   theme_light() + 
  #   theme(panel.grid.major = element_blank(),
  #         panel.grid.minor = element_blank()) + 
  #   labs(x = "Propensity Score",y = "Density") + 
  #   scale_colour_brewer(palette="Set1") +
  #   facet_grid(model ~ .) + 
  #   geom_density()
  # #Accent, Dark2, Paired, Pastel1, Pastel2, Set1, Set2, Set3
  # dev.off()
  
  # muT<-data.frame(muhatPT,dat$y);muT$model<-1;names(muT)[1]<-"mu"
  # head(muT)
  # muF<-data.frame(muhatPF,dat$y);muF$model<-2;names(muF)[1]<-"mu"
  # head(muF)
  # muNPT<-data.frame(muhatNPT,dat$y);muNPT$model<-3;names(muNPT)[1]<-"mu"
  # head(muNPT)
  # muNPF<-data.frame(muhatNPF,dat$y);muNPF$model<-4;names(muNPF)[1]<-"mu"
  # head(muNPF)
  # muDat<-rbind(muT,muF,muNPT,muNPF)
  # head(muDat)
  # muDat$model<-factor(muDat$model,labels=c("Parametric Corr","Parametric Miss","Nonparametric Corr","Nonparametric Miss"))
  # tail(muDat)
  # 
  # pdf(file="~/Dropbox/Documents/Research/Papers/SemiparametricInference/eFigure3.pdf",height=10,width=6)
  # ggplot(muDat,aes(mu,dat.y)) +
  #   theme_light() +
  #   theme(panel.grid.major = element_blank(),
  #         panel.grid.minor = element_blank()) +
  #   labs(x = "Predicted Y",y = "Observed Y") +
  #   scale_colour_brewer(palette="Set1") +
  #   facet_grid(model ~ .) +
  #   geom_point() + ylim(0,200) + xlim(0,200) +
  #   geom_abline(slope=1, intercept=c(0,0),color="red",
  #               linetype="dashed",size=.25)
  # #Accent, Dark2, Paired, Pastel1, Pastel2, Set1, Set2, Set3
  # dev.off()
  
  # compute estimators
  res.est$ipwPMT[i] <- coef(lm(y~r,data=dat,weights=swPT))[2] 
  res.est$ipwNPT[i] <- coef(lm(y~r,data=dat,weights=swNPT))[2] 
  res.est$regPMT[i] <- mean(muhatPT1 - muhatPT0)
  res.est$regNPT[i] <- mean(muhatNPT1 - muhatNPT0)
  
  res.est$aipwPMT[i]  <- mean((((2*dat$r-1)*(dat$y - muhatPT))/((2*dat$r-1)*pihatPT + (1-dat$r)) + muhatPT1 - muhatPT0))
  res.est$aipwNPT[i]  <- mean((((2*dat$r-1)*(dat$y - muhatNPT))/((2*dat$r-1)*pihatNPT + (1-dat$r)) + muhatNPT1 - muhatNPT0))
  
  res.est$tmlePMT[i]<-tmlePMT$estimates$ATE$psi
  res.est$tmleNPT[i]<-tmleNPT$estimates$ATE$psi
  
  res.est$ipwPMF[i] <- coef(lm(y~r,data=dat,weights=swPF))[2] 
  res.est$ipwNPF[i] <- coef(lm(y~r,data=dat,weights=swNPF))[2] 
  res.est$regPMF[i] <- mean(muhatPF1 - muhatPF0)
  res.est$regNPF[i] <- mean(muhatNPF1 - muhatNPF0)
  
  res.est$aipwPMF[i]  <- mean((((2*dat$r-1)*(dat$y - muhatPF))/((2*dat$r-1)*pihatPF + (1-dat$r)) + muhatPF1 - muhatPF0))
  res.est$aipwNPF[i]  <- mean((((2*dat$r-1)*(dat$y - muhatNPF))/((2*dat$r-1)*pihatNPF + (1-dat$r)) + muhatNPF1 - muhatNPF0))
  
  res.est$tmlePMF[i]<-tmlePMF$estimates$ATE$psi
  res.est$tmleNPF[i]<-tmleNPF$estimates$ATE$psi
  
  # compute closed-form SEs
  mod1<-lm(y~r,data=dat,weights=swPT)
  mod2<-lm(y~r,data=dat,weights=swNPT)
  res.se$ipwPMT[i] <- vcovHC(mod1)[2,2] 
  res.se$ipwNPT[i] <- vcovHC(mod2)[2,2] 
  res.se$aipwPMT[i] <- sd((((2*dat$r-1)*(dat$y - muhatPT))/((2*dat$r-1)*pihatPT + (1-dat$r)) + muhatPT1 - muhatPT0))
  res.se$aipwNPT[i] <- sd((((2*dat$r-1)*(dat$y - muhatNPT))/((2*dat$r-1)*pihatNPT + (1-dat$r)) + muhatNPT1 - muhatNPT0))
  res.se$tmlePMT[i] <- sqrt(tmlePMT$estimates$ATE$var.psi)*sqrt(n)
  res.se$tmleNPT[i] <- sqrt(tmleNPT$estimates$ATE$var.psi)*sqrt(n)
  
  mod1<-lm(y~r,data=dat,weights=swPF)
  mod2<-lm(y~r,data=dat,weights=swNPF)
  res.se$ipwPMF[i] <- vcovHC(mod1)[2,2] 
  res.se$ipwNPF[i] <- vcovHC(mod2)[2,2] 
  res.se$aipwPMF[i] <- sd((((2*dat$r-1)*(dat$y - muhatPF))/((2*dat$r-1)*pihatPF + (1-dat$r)) + muhatPF1 - muhatPF0))
  res.se$aipwNPF[i] <- sd((((2*dat$r-1)*(dat$y - muhatNPF))/((2*dat$r-1)*pihatNPF + (1-dat$r)) + muhatNPF1 - muhatNPF0))
  res.se$tmlePMF[i] <- sqrt(tmlePMF$estimates$ATE$var.psi)*sqrt(n)
  res.se$tmleNPF[i] <- sqrt(tmleNPF$estimates$ATE$var.psi)*sqrt(n)
  
  # compute bootstrap CIs
  datB<- data.frame(z,x,r,y); colnames(datB)[1:(2*ncol(z))] <- c(paste("z",1:ncol(x),sep=""),paste("x",1:ncol(x),sep=""))
  cat("Bootstrapping Parametric",'\n');flush.console()
  plugin1 <- function(d,j,dim){
    x<-d[j,c(paste0("x",1:p))]
    mMat<-model.matrix(as.formula(paste("~(",paste("x[,",1:dim,"]",collapse="+"),")")))[,-1]
    dat<-data.table(y=d[j,]$y,r=d[j,]$r,mMat);mumod<-glm(y~.,data=dat)
    dat$r<-1;muhat1 <- predict(mumod,newdata=dat)
    dat$r<-0;muhat0 <- predict(mumod,newdata=dat)
    meanPT<-mean(muhat1-muhat0)
    x<-d[j,c(paste0("z",1:p))]
    mMat<-model.matrix(as.formula(paste("~(",paste("x[,",1:dim,"]",collapse="+"),")")))[,-1]
    dat<-data.table(y=d[j,]$y,r=d[j,]$r,mMat);mumod<-glm(y~.,data=dat)
    dat$r<-1;muhat1 <- predict(mumod,newdata=dat)
    dat$r<-0;muhat0 <- predict(mumod,newdata=dat)
    meanPF<-mean(muhat1-muhat0)
    res<-c(meanPT,meanPF)
  }
  bs1 <- boot(datB,plugin1,R=bootNum,dim=p)
  res.se$regPMT[i] <- sd(bs1$t[,1])*sqrt(n)
  res.se$regPMF[i] <- sd(bs1$t[,2])*sqrt(n)
    
  if (bs==T){
    cat("Bootstrapping Nonparametric",'\n');flush.console()
    plugin2<-function(d,dim,f){
      j<-sample(1:nrow(d),nrow(d),replace=T);dat<-d[j,]
      x<-dat[,c(paste0("x",1:p))]
      W<-model.matrix(as.formula(paste("~(",paste("x[,",1:dim,"]",collapse="+"),")")))[,-1]
      Y<-dat$y;A<-as.matrix(dat$r) 
      X=data.frame(cbind(A,W))
      names(X)<-c("A",paste0("x",1:ncol(W)))
      mumod<-SuperLearner(Y,X,family=gaussian,SL.library=sl.lib,cvControl=list(V=f))
      X$A<-1;muhat1 <- predict(mumod,newdata=X,onlySL=T)$pred
      X$A<-0;muhat0 <- predict(mumod,newdata=X,onlySL=T)$pred
      gT<-mean(muhat1-muhat0)
      x<-dat[,c(paste0("z",1:p))]
      W<-model.matrix(as.formula(paste("~(",paste("x[,",1:dim,"]",collapse="+"),")")))[,-1]
      Y<-dat$y;A<-as.matrix(dat$r) 
      X=data.frame(cbind(A,W))
      names(X)<-c("A",paste0("x",1:ncol(W)))
      mumod<-SuperLearner(Y,X,family=gaussian,SL.library=sl.lib,cvControl=list(V=f))
      X$A<-1;muhat1 <- predict(mumod,newdata=X,onlySL=T)$pred
      X$A<-0;muhat0 <- predict(mumod,newdata=X,onlySL=T)$pred
      gF<-mean(muhat1-muhat0)
      g<-c(gT,gF)
      return(g)
    }
    coreNum<-detectCores()
    cl <- makeCluster(coreNum)
    registerDoParallel(cl)
    clusterCall(cl, function(x) .libPaths(x), .libPaths())
    clusterExport(cl, sl.lib)
    pak<-c("tmle","SuperLearner","glmnet","rpart","ranger","nnet","arm","earth")
    results<-foreach(i=1:bootNum,.packages=pak) %dorng% {
      b<-plugin2(dat,dim=p,f=folds);b
    }
    stopCluster(cl)
    #results<-lapply(1:bootNum,function(x) plugin2(dat,dim=p,f=folds,misspec=mSpec))
    results<-do.call(rbind,results)
    res.se$regNPT[i] <- sd(results[,1])*sqrt(n)
    res.se$regNPF[i] <- sd(results[,2])*sqrt(n)
  }
  # print updating results
  res.cov <- res.est-1.96*res.se/sqrt(n) < true & true < res.est+1.96*res.se/sqrt(n)
  res.width <- (res.est+1.96*res.se/sqrt(n)) - (res.est-1.96*res.se/sqrt(n))
  tmp <- data.frame(rbind(c(n,apply(res.est-true,2,mean,na.rm=T)),
               c(n,apply((res.est-true)^2,2,mean,na.rm=T)),
               c(n,apply(res.cov,2,mean,na.rm=T)),
               c(n,apply(res.width,2,mean,na.rm=T))))
  rownames(tmp)<-c("bias","rmse","cov","width");colnames(tmp)[1]<-"N";print(round(tmp,3));cat('\n');flush.console()
  setDT(tmp, keep.rownames = TRUE)[];colnames(tmp)[1] <- "type"
  
  if(i==1&samp==1){
    write.table(tmp,"results.txt",sep="\t",row.names=F)
  } else{
    write.table(tmp,"results.txt",sep="\t",row.names=F,col.names=F,append=T)
  }
  return(tmp)
}

results<-lapply(1:(nsim*5), function(x) npDR(x,pNum=4,bs=F,bootNum=100))

cores<-detectCores()
print(cores)
results<-mclapply(1:(nsim*5), function(x) npDR(x,pNum=4,bs=F,bootNum=100),mc.cores=cores)
results<-do.call(rbind,results)
