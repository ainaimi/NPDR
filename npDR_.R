# load libraries/functions
#args=(commandArgs(TRUE))
#setwd(".")

packages <- c("foreach","doParallel","doRNG","boot","rmutil","mvtnorm","gam","sandwich","ggplot2",
              "devtools","glmnet","data.table","rpart","ranger","nnet","arm","earth","e1071")
# userLib <-  "~/R/R_LIBS_USER"
# .libPaths(userLib)

# ifelse(!dir.exists(userLib), dir.create(userLib), FALSE)

for (package in packages) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package, repos='http://lib.stat.cmu.edu/R/CRAN')
  }
}

for (package in packages) {
  library(package, character.only=T)
}

devtools::install_github("ecpolley/SuperLearner")
library(SuperLearner)

install.packages("tmle",repos='http://lib.stat.cmu.edu/R/CRAN')
library(tmle)

time<-proc.time()
expit <- function(x){ exp(x)/(1+exp(x)) }
logit <- function(x){ log(x/(1-x)) }
nsim <- 1

cols <- c("ipwPMT","ipwNPT","ipwPMF","ipwNPF",
          "regPMT","regNPT","regPMF","regNPF",
          "aipwPMT","aipwNPT","aipwPMF","aipwNPF",
          "tmlePMT","tmleNPT","tmlePMF","tmleNPF")
res.est <- data.frame(matrix(nrow=nsim*3,ncol=length(cols)));colnames(res.est) <- cols; res.se <- res.est

ranger_learner <- create.Learner("SL.ranger", list(min.node.size = c(30,60)))
xgboost_learner <- create.Learner("SL.xgboost", list(minobspernode = c(30,60)))

sl.lib <- c(ranger_learner$names,xgboost_learner$names)

##  true value
true<-6;true
npDR<-function(counter,pNum,bs=T,bootNum=100){
  set.seed(counter)
  # data management
  i<-counter
  samp<-rep(1:3,nsim*3)[counter]
  ss<-c(50,200,1200)
  n<-ss[samp]
  p=pNum
  cat("Now running iteration",i,"with a sample size of",n,'\n');flush.console()
  
  # confounders
  sigma<-matrix(0,nrow=p,ncol=p);diag(sigma)<-1
  x <- rmvnorm(n, mean=rep(0,p), sigma=sigma)
  
  z<-x
  z[,1]<-exp(x[,1]/2)
  z[,2]<-x[,2]/(1+exp(x[,1]))+10
  z[,3]<-(x[,1]*x[,3]/25+.6)^3
  z[,4]<-(x[,2]*x[,4]+20)^2

   # pdf("eFigure1.pdf",width=6,height=6)
   # GGally::ggpairs(data.frame(z))
   # dev.off()
  
  # design matrix for outcome model

  muMatT<-model.matrix(as.formula(paste("~(",paste("x[,",1:ncol(x),"]",collapse="+"),")")))
  
  parms3<-c(3.5,2.5,-1,5) #,4.25,-2
  parms4<-c(log(2),log(2.5),log(.5),log(1.5)) #,log(2.25),log(.25)
  
  beta<-parms3;beta<-c(120,beta)
  # design matrix for propensity score model
  piMatT<-model.matrix(as.formula(paste("~(",paste("x[,",1:ncol(x),"]",collapse="+"),")")))
  theta<-parms4;theta<-c(-.5,theta)
  mu <- muMatT%*%beta
  # propensity score model
  pi <- expit(piMatT%*%theta);
  r<-1-rbinom(n,1,pi)
  # outcome model: true expsoure effect = 6
  y <- r*6 + mu + rnorm(n,0,6)
  
  # induce misspecification
  muMatF<-model.matrix(as.formula(paste("~(",paste("z[,",1:ncol(x),"]",collapse="+"),")")))
  piMatF<-model.matrix(as.formula(paste("~",paste("z[,",1:ncol(x),"]",collapse="+"))))
  
  # correct specification
  dat <- data.frame(x,r,y); colnames(dat)[1:ncol(x)] <- paste("x",1:ncol(x),sep="")
  W <- muMatT[,-1];colnames(W) <- paste("W",1:ncol(muMatT[,-1]), sep="");A <- r;Y <- y
  tQForm<-as.formula(paste0("Y~A+", paste(paste0("W",1:ncol(muMatT[,-1])), collapse="+")))
  tgForm<-as.formula(paste0("A~", paste(paste0("W",1:ncol(piMatT[,-1])), collapse="+")))
  tmlePMT <- tmle(Y,A,W,family="gaussian",Qform=tQForm,gform=tgForm)
  
  folds<-c(20,10,5)[samp]
  cat("Number of cross-validation folds is",folds,'\n');flush.console()
  tmleNPT <- tmle(Y,A,W,family="gaussian",Q.SL.library=sl.lib,g.SL.library=sl.lib)
  
  tmleNPT$g$coef
  
  pihatPT <- tmlePMT$g$g1W
  pihatPT <- ifelse(pihatPT < quantile(pihatPT,c(0.025)),quantile(pihatPT,c(0.025)),pihatPT)
  pihatPT <- ifelse(pihatPT > quantile(pihatPT,c(1-0.025)),quantile(pihatPT,c(1-0.025)),pihatPT)  
  
  swPT<-dat$r*(mean(dat$r)/tmlePMT$g$g1W) + (1-dat$r)*((1-mean(dat$r))/(1-tmlePMT$g$g1W))
  pihatNPT <- tmleNPT$g$g1W
  pihatNPT <- ifelse(pihatNPT < quantile(pihatNPT,c(0.025)),quantile(pihatNPT,c(0.025)),pihatNPT)
  pihatNPT <- ifelse(pihatNPT > quantile(pihatNPT,c(1-0.025)),quantile(pihatNPT,c(1-0.025)),pihatNPT)  
  
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
  
  pihatPF <- tmlePMF$g$g1W
  pihatPF <- ifelse(pihatPF < quantile(pihatPF,c(0.025)),quantile(pihatPF,c(0.025)),pihatPF)
  pihatPF <- ifelse(pihatPF > quantile(pihatPF,c(1-0.025)),quantile(pihatPF,c(1-0.025)),pihatPF)    
    
  swPF<-dat$r*(mean(dat$r)/tmlePMF$g$g1W) + (1-dat$r)*((1-mean(dat$r))/(1-tmlePMF$g$g1W))
  pihatNPF <- tmleNPF$g$g1W
  pihatNPF <- ifelse(pihatNPF < quantile(pihatNPF,c(0.025)),quantile(pihatNPF,c(0.025)),pihatNPF)
  pihatNPF <- ifelse(pihatNPF > quantile(pihatNPF,c(1-0.025)),quantile(pihatNPF,c(1-0.025)),pihatNPF)    
  
  swNPF<-dat$r*(mean(dat$r)/tmleNPF$g$g1W) + (1-dat$r)*((1-mean(dat$r))/(1-tmleNPF$g$g1W))
  muhatPF  <- tmlePMF$Qinit$Q[,2]*A+tmlePMF$Qinit$Q[,1]*(1-A)
  muhatPF1 <- tmlePMF$Qinit$Q[,2];muhatPF0 <- tmlePMF$Qinit$Q[,1]
  muhatNPF  <- tmleNPF$Qinit$Q[,2]*A+tmleNPF$Qinit$Q[,1]*(1-A)
  muhatNPF1 <- tmleNPF$Qinit$Q[,2];muhatNPF0 <- tmleNPF$Qinit$Q[,1]
  
  # #propensity score plot
  # plotDat <- rbind(cbind(data.frame(type="Parametric True",PS=pihatPT,Exposure=as.factor(dat$r),mu=muhatPT,Outcome=dat$y)),
  #                  cbind(data.frame(type="Parametric Misspecified",PS=pihatPF,Exposure=as.factor(dat$r),mu=muhatPF,Outcome=dat$y)),
  #                  cbind(data.frame(type="Nonparametric True",PS=pihatNPT,Exposure=as.factor(dat$r)),mu=muhatNPT,Outcome=dat$y),
  #                  cbind(data.frame(type="Nonarametric Misspecified",PS=pihatNPF,Exposure=as.factor(dat$r),mu=muhatNPF,Outcome=dat$y)))
  # 
  # group.colors <- c(`1` = "red", `0` = "blue")
  # pdf("eFigure1.pdf",width=6,height=6)
  # ggplot(plotDat) + 
  #   geom_density(aes(x=PS,group=Exposure,color=Exposure)) + 
  #   facet_wrap(~type) + 
  #   scale_color_manual(values=group.colors) +
  #   xlab("Propensity Score") + ylab("Density")
  # dev.off()
  # 
  # pdf("eFigure2.pdf",width=6,height=6)
  # ggplot(plotDat) + 
  #   geom_point(aes(x=mu,y=Outcome),size=.5,alpha=.25) + 
  #   facet_wrap(~type) + 
  #   xlab("Predicted Outcome") + ylab("Observed Outcome")
  # dev.off()
  
  # compute estimators
  res.est$ipwPMT[i] <- coef(lm(y~r,data=dat,weights=swPT))[2] 
  res.est$ipwNPT[i] <- coef(lm(y~r,data=dat,weights=swNPT))[2] 
  res.est$regPMT[i] <- mean(muhatPT1 - muhatPT0)
  res.est$regNPT[i] <- mean(muhatNPT1 - muhatNPT0)
  
  lower_bound <- min(y) - max(y)
  upper_bound <- max(y) - min(y)
  
  aipwPMT <- mean((((2*dat$r-1)*(dat$y - muhatPT))/((2*dat$r-1)*pihatPT + (1-dat$r)) + muhatPT1 - muhatPT0))
  aipwPMT <- ifelse(aipwPMT>upper_bound,upper_bound,aipwPMT)
  aipwPMT <- ifelse(aipwPMT<lower_bound,lower_bound,aipwPMT)
  
  aipwNPT <- mean((((2*dat$r-1)*(dat$y - muhatNPT))/((2*dat$r-1)*pihatNPT + (1-dat$r)) + muhatNPT1 - muhatNPT0))
  aipwNPT <- ifelse(aipwNPT>upper_bound,upper_bound,aipwNPT)
  aipwNPT <- ifelse(aipwNPT<lower_bound,lower_bound,aipwNPT)
  
  res.est$aipwPMT[i]  <- aipwPMT
  res.est$aipwNPT[i]  <- aipwNPT
  
  res.est$tmlePMT[i]<-tmlePMT$estimates$ATE$psi
  res.est$tmleNPT[i]<-tmleNPT$estimates$ATE$psi
  
  res.est$ipwPMF[i] <- coef(lm(y~r,data=dat,weights=swPF))[2] 
  res.est$ipwNPF[i] <- coef(lm(y~r,data=dat,weights=swNPF))[2] 
  res.est$regPMF[i] <- mean(muhatPF1 - muhatPF0)
  res.est$regNPF[i] <- mean(muhatNPF1 - muhatNPF0)
  
  ## truncate PS for AIPW and IPW 
  ## bound AIPW estimator uusing \psi \in [y_min - y_max, y_max - y_min]
  
  aipwPMF <- mean((((2*dat$r-1)*(dat$y - muhatPF))/((2*dat$r-1)*pihatPF + (1-dat$r)) + muhatPF1 - muhatPF0))
  aipwPMF <- ifelse(aipwPMF>upper_bound,upper_bound,aipwPMF)
  aipwPMF <- ifelse(aipwPMF<lower_bound,lower_bound,aipwPMF)
  
  aipwNPF <- mean((((2*dat$r-1)*(dat$y - muhatNPF))/((2*dat$r-1)*pihatNPF + (1-dat$r)) + muhatNPF1 - muhatNPF0))
  aipwNPF <- ifelse(aipwNPF>upper_bound,upper_bound,aipwNPF)
  aipwNPF <- ifelse(aipwNPF<lower_bound,lower_bound,aipwNPF)
  
  
  res.est$aipwPMF[i]  <- aipwPMF
  res.est$aipwNPF[i]  <- aipwNPF
  
  res.est$tmlePMF[i]<-tmlePMF$estimates$ATE$psi
  res.est$tmleNPF[i]<-tmleNPF$estimates$ATE$psi
  
  # compute closed-form SEs
  mod1<-lm(y~r,data=dat,weights=swPT)
  mod2<-lm(y~r,data=dat,weights=swNPT)
  res.se$ipwPMT[i] <- sqrt(vcovHC(mod1,type = "HC")[2,2])
  res.se$ipwNPT[i] <- sqrt(vcovHC(mod2,type = "HC")[2,2])
  res.se$aipwPMT[i] <- sd((((2*dat$r-1)*(dat$y - muhatPT))/((2*dat$r-1)*pihatPT + (1-dat$r)) + muhatPT1 - muhatPT0))/sqrt(n)
  res.se$aipwNPT[i] <- sd((((2*dat$r-1)*(dat$y - muhatNPT))/((2*dat$r-1)*pihatNPT + (1-dat$r)) + muhatNPT1 - muhatNPT0))/sqrt(n)
  res.se$tmlePMT[i] <- sqrt(tmlePMT$estimates$ATE$var.psi)
  res.se$tmleNPT[i] <- sqrt(tmleNPT$estimates$ATE$var.psi)
  
  mod1<-lm(y~r,data=dat,weights=swPF)
  mod2<-lm(y~r,data=dat,weights=swNPF)
  res.se$ipwPMF[i] <- sqrt(vcovHC(mod1,type="HC")[2,2])
  res.se$ipwNPF[i] <- sqrt(vcovHC(mod2,type="HC")[2,2])
  res.se$aipwPMF[i] <- sd((((2*dat$r-1)*(dat$y - muhatPF))/((2*dat$r-1)*pihatPF + (1-dat$r)) + muhatPF1 - muhatPF0))/sqrt(n)
  res.se$aipwNPF[i] <- sd((((2*dat$r-1)*(dat$y - muhatNPF))/((2*dat$r-1)*pihatNPF + (1-dat$r)) + muhatNPF1 - muhatNPF0))/sqrt(n)
  res.se$tmlePMF[i] <- sqrt(tmlePMF$estimates$ATE$var.psi)
  res.se$tmleNPF[i] <- sqrt(tmleNPF$estimates$ATE$var.psi)
  
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
  res.se$regPMT[i] <- sd(bs1$t[,1])
  res.se$regPMF[i] <- sd(bs1$t[,2])
    
  if (bs==T){
    cat("Bootstrapping Nonparametric",'\n');flush.console()
    plugin2<-function(d,j,dim,f){
      x<-d[j,c(paste0("x",1:p))]
      W<-model.matrix(as.formula(paste("~(",paste("x[,",1:dim,"]",collapse="+"),")")))[,-1]
      Y<-d[j,]$y;A<-as.matrix(d[j,]$r) 
      X=data.frame(cbind(A,W))
      names(X)<-c("A",paste0("x",1:ncol(W)))
      mumod<-SuperLearner(Y,X,family=gaussian,SL.library=sl.lib,cvControl=list(V=f))
      X$A<-1;muhat1 <- predict(mumod,newdata=X,onlySL=T)$pred
      X$A<-0;muhat0 <- predict(mumod,newdata=X,onlySL=T)$pred
      gT<-mean(muhat1-muhat0)
      
      x<-d[j,c(paste0("z",1:p))]
      W<-model.matrix(as.formula(paste("~(",paste("x[,",1:dim,"]",collapse="+"),")")))[,-1]
      Y<-d[j,]$y;A<-as.matrix(d[j,]$r) 
      X=data.frame(cbind(A,W))
      names(X)<-c("A",paste0("x",1:ncol(W)))
      mumod<-SuperLearner(Y,X,family=gaussian,SL.library=sl.lib,cvControl=list(V=f))
      X$A<-1;muhat1 <- predict(mumod,newdata=X,onlySL=T)$pred
      X$A<-0;muhat0 <- predict(mumod,newdata=X,onlySL=T)$pred
      gF<-mean(muhat1-muhat0)
      g<-c(gT,gF)
      return(g)
    }
    
    bs2 <- boot(datB,plugin2,R=bootNum,dim=p,f=folds)
    res.se$regNPT[i] <- sd(bs2$t[,1])
    res.se$regNPF[i] <- sd(bs2$t[,2])
  }
  # print updating results
  res.cov <- res.est-1.96*res.se < true & true < res.est+1.96*res.se
  res.width <- (res.est+1.96*res.se) - (res.est-1.96*res.se)
  tmp <- data.frame(rbind(c(n,apply(res.est-true,2,mean,na.rm=T)),
                          c(n,apply((res.est-true)^2,2,mean,na.rm=T)),
                          c(n,apply(res.cov,2,mean,na.rm=T)),
                          c(n,apply(res.width,2,mean,na.rm=T))))
  tmp.se <- data.frame(rbind(c(n,apply(res.se,2,mean,na.rm=T))))
  rownames(tmp)<-c("bias","rmse","cov","width");colnames(tmp)[1]<-"N";print(round(tmp,3));cat('\n');flush.console()
  colnames(tmp.se)[1]<-"N"
  setDT(tmp, keep.rownames = TRUE)[];colnames(tmp)[1] <- "type"
  
  if(i==1&samp==1){
    write.table(tmp,"results.txt",sep="\t",row.names=F)
    write.table(tmp.se,"results_se.txt",sep="\t",row.names=F)
  } else{
    write.table(tmp,"results.txt",sep="\t",row.names=F,col.names=F,append=T)
    write.table(tmp.se,"results_se.txt",sep="\t",row.names=F,col.names=F,append=T)
  }
  return(tmp)
}

start_time <- Sys.time()
cores<-detectCores()
print(cores)
results<-mclapply(1:(nsim*3), function(x) npDR(x,pNum=4,bs=T,bootNum=100),mc.cores=cores,mc.set.seed=F)
## run time:
Sys.time() - start_time