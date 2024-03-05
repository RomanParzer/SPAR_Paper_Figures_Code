#####################################################
### Testing sparse projected averaged regression SPAR
#####################################################

pacman::p_load(foreach,parallel,tidyr,dplyr)
source("../functions/multi_assign.R")
source("../functions/methods.R")

require(Matrix)
rateyedata <- data.frame(t(read.csv("../data/RatEyeAll.txt",sep="\t", header=TRUE, dec=".",row.names = 1)))
myq <- quantile(as.matrix(rateyedata),c(0.25))
out_q <- apply(rateyedata,2,function(mycol) max(mycol)<myq)
out_v <- apply(rateyedata,2, function(x) sd(x) / mean(x) * 100)
out_ind <- which((out_v < 2) | out_q)
# remove probes not sufficiently expressed and controls
mydata <- data.frame(rateyedata[,-unique(c(31043:31099,out_ind))])

generate_data_rateye <- function(mydata,p=ncol(mydata),ntest=0,snr=2,a=NULL,ind=NULL,beta=NULL,mu=1) {
  stopifnot(ncol(mydata)>=p)
  stopifnot(ntest<nrow(mydata))
  # first select active vars and generate beta
  bSb <- 0
  n <- nrow(mydata) - ntest
  test_ind <- sample(1:(n+ntest),ntest)
  
  if (is.null(beta)) {
    if (is.null(a)) {
      if (is.null(ind)) {
        a <- round(4*log(p))
      } else {
        a <- length(ind)
      }
    }
    if (is.null(ind)) {
      ind <- sample(1:p,a,replace=FALSE)
    }
    beta <- Matrix(data=c(0),p,1,sparse = TRUE)
    beta[ind] <- sample(c(1,-1),a,replace = TRUE,prob = c(0.6,0.4)) * (4*log(n)/sqrt(n) + abs(rnorm(a)))
  } else {
    ind <- which(beta!=0)
    a <- length(ind)
  }
  
  # then generate x with given cov structure and calculate bSb = beta' Sigma beta
  # efficient implementation for high p
  x <- as.matrix(mydata[,1:p])
  bSb <- crossprod(beta[ind],cov(x[,ind]))%*%beta[ind]
  stopifnot(nrow(x)==(n+ntest),ncol(x)==p)
  
  # finally generate responses
  sigma2 <- bSb/snr
  y <- mu + x[,ind,drop=FALSE]%*%beta[ind] + rnorm(n+ntest,0,sqrt(sigma2))
  
  return(list(x=x[-test_ind,],y=y[-test_ind],xtest=x[test_ind,],ytest=y[test_ind],
              ind=ind,sigma2=sigma2,beta=beta))
}

################ adapt #########################################################################

nrep <- 100

methods <- list("HOLP"=myOLS,
                "PCR"=myPCR,
                "PLS"=myPLS,
                "Ridge"=function(x,y,xtest){myElNet(x,y,xtest,alpha=0)},
                "LASSO"=function(x,y,xtest){myElNet(x,y,xtest,alpha=1)},
                "AdLASSO"=myAdLASSO,
                "ElNet"=myElNet,
                "SIS"=function(x,y,xtest){mySIS(x,y,xtest,iter=FALSE)},
                "RP_Gaus"=function(x,y,xtest){myRP(x,y,xtest,RP="Gaussian",nummod = 1)},
                "RP_Sparse"=function(x,y,xtest){myRP(x,y,xtest,RP="Sparse",nummod = 1)},
                "RP_CW"=function(x,y,xtest){myRP(x,y,xtest,RP="CWSparse",nummod = 1)},
                "RP_CW_Ensemble"=function(x,y,xtest){myRP(x,y,xtest,RP="CWSparse",nummod = 100)},
                "TARP"=myTARP,
                "SPAR"=function(x,y,xtest){mySPAR(x,y,xtest,nummods=20,nlambda=1)},
                "SPAR CV"=function(x,y,xtest){mySPAR(x,y,xtest,nummods=c(10,20,30,50,100),opt_par = "best")},
                "SPAR CV 1se"=function(x,y,xtest){mySPAR(x,y,xtest,nummods=c(10,20,30,50,100),opt_par = "1se")})


measures <- c("rMSPE","rMSPE_tr","Time","Precision","Recall","activeEE","passiveEE","numAct","AUC","pAUC")

simulation_settings <- expand_grid(p=c(500,2000,ncol(mydata)),ntest=round(nrow(mydata)/4), snr=c(10),
                                   act_setting=c("sparse","medium","dense"))
simulation_settings <- rbind(simulation_settings,expand_grid(p=2000,ntest=round(nrow(mydata)/4), snr=c(1,5,20),
                                                             act_setting="medium"))
simulation_settings <- simulation_settings %>% mutate(a = round(ifelse(act_setting=="sparse",2*log(p),
                                                                       ifelse(act_setting=="medium",nrow(mydata)/2+2*log(p),p/4))))

################## simulations ##############################################################################
unlink("../saved_results/log.txt")
# select number of cores
my.cluster <- parallel::makeCluster(4, type = "PSOCK", outfile = "../saved_results/log.txt")
doParallel::registerDoParallel(cl = my.cluster)
# foreach::getDoParRegistered()

nmethods <- length(methods)
nset <- nrow(simulation_settings)
nmeas <- length(measures)
res <- array(c(0),dim=c(nrep,nmeas,nmethods,nset),dimnames = list(reps=NULL,measures=measures,
                                                              method=names(methods),settings=paste("Scenario",1:nset)))
attributes(res)$settings <- simulation_settings
fits <- vector(mode="list",length = nset)

clusterExport(my.cluster,c('simulation_settings','nmeas','nmethods','methods','mydata','generate_data_rateye'), envir = environment())
clusterEvalQ(my.cluster, {  
  source("../functions/multi_assign.R")
  source("../functions/methods.R")
})

parres <- foreach(j = 1:nset) %:%
foreach(i=1:nrep) %dopar% {
    parresi <- list(res=matrix(c(0),nmeas,nmethods),fits = vector(mode="list",length = nmethods))
    names(parresi$fits) <- names(methods)
  
    c(p,ntest,snr,act_setting,a) %<-% simulation_settings[j,]
    set.seed((1234+i)^2)
    data <- generate_data_rateye(mydata,p,ntest,snr=snr,a=a)
    x <- data$x
    y <- data$y
    xtest <- data$xtest
    ytest <- data$ytest
    rMSPEconst <- mean((ytest-mean(y))^2)
    rMSPEconst_tr <- var(y)

    for (k in 1:nmethods) {
      set.seed((1234+i)^2 + k)
      tstamp <- Sys.time()
      callres <- tryCatch( methods[[k]](x,y,xtest),
                                                error=function(error_message) {
                                                  message(paste0("Error in ",names(methods)[k],": ",error_message))
                                                  return(NULL)
                                                })
      tmp_time <- as.numeric(Sys.time() - tstamp,units="secs")
      
      # plot(performance(prediction(as.numeric(abs(callres$beta)),as.numeric(data$beta!=0)),measure="tpr", x.measure="fpr"))
      
      if (is.null(callres$beta)) {
        parresi$res[,k] <- c(mean((callres$yhat-ytest)^2)/rMSPEconst, # rMSPE
                             mean((callres$yhat_tr-y)^2)/rMSPEconst_tr, # rMSPE_tr
                             tmp_time, # time
                             NA, # precision
                             NA, # recall
                             NA, # activeEE
                             NA, # passiveEE
                             NA, # numAct
                             NA, # AUC
                             NA) # pAUC
        
        
        
      } else {
        if ("dgCMatrix" %in% class(callres$beta)) {
          tmp_ind <- callres$beta@i+1
        } else {
          tmp_ind <- which(callres$beta!=0)
        }
        parresi$res[,k] <- c(mean((callres$yhat-ytest)^2)/rMSPEconst, # rMSPE
                             mean((callres$yhat_tr-y)^2)/rMSPEconst_tr, # rMSPE_tr
                             tmp_time, # time
                             ifelse(length(tmp_ind)==0,0,mean(tmp_ind %in% data$ind)), # precision
                             mean(data$ind %in% tmp_ind), # recall
                             mean((callres$beta[data$ind]-data$beta[data$ind])^2), # activeEE
                             mean(callres$beta[-data$ind]^2), # passiveEE
                             length(tmp_ind), # numAct
                             performance(prediction(as.numeric(abs(callres$beta)),as.numeric(data$beta!=0)),measure="auc")@y.values[[1]], # AUC
                             performance(prediction(as.numeric(abs(callres$beta)),as.numeric(data$beta!=0)),measure="auc",fpr.stop=nrow(mydata)/2/(p-a))@y.values[[1]]) # pAUC
      }
      
      
      # message(sprintf("Finished method %d / %d!",k,nmethods))
    }
    res[i,,,j] <- parresi$res
    res[i,,,j]
    cat(sprintf('Finished rep %d / %d for setting %d / %d at %s.\n',i,nrep,j,nset,Sys.time()))
    warnings()
    parresi
}

for (j in 1:nset) {
  # fits[[j]] <- vector("list",nrep)
  for (i in 1:nrep) {
    res[i,,,j] <- parres[[j]][[i]]$res
    # fits[[j]][[i]] <- parres[[j]][[i]]$fits
  }
}

saveRDS(list(res=res,fits=fits),
        file = sprintf("../saved_results/synth_rateye_SPAR_nset%d_reps%d_nmeth%d.rds",nset,nrep,nmethods))
parallel::stopCluster(cl = my.cluster)

warnings()
