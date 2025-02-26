#####################################################
### Testing sparse projected averaged regression SPAR
#####################################################

pacman::p_load(foreach,parallel,tidyr,dplyr,ROCR)
source("../functions/multi_assign.R")
source("../functions/data_generation.R")
source("../functions/methods.R")

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
                "SPAR CV 1se"=function(x,y,xtest){mySPAR(x,y,xtest,nummods=c(10,20,30,50,100),opt_par = "1se")},
                "HOLPScr"=myHOLPScr,
                "RF"=myRF,
                "OWL"=myOWL,
                "BAMP" = function(x,y,xtest){myBAMP(x,y,xtest,gamma=0.8)})


measures <- c("rMSPE","rMSPE_tr","Time","Precision","Recall","activeEE","passiveEE","numAct","AUC","pAUC")

simulation_settings <- expand_grid(n=200, p=c(500,2000,10^4), act_setting=c("sparse","medium","dense"), 
                                   ntest=1000, snr=c(10),  cov_setting=c("group"))

simulation_settings <- rbind(simulation_settings,expand_grid(n=200, p=c(2000), act_setting=c("sparse","medium","dense"), 
                                   ntest=1000, snr=c(10),  cov_setting=c("ind","comsym","ar1","group","factor","extreme")))

simulation_settings <- rbind(simulation_settings,expand_grid(n=c(100,400,800), p=2000, act_setting=c("sparse","medium","dense"), 
                                                             ntest=1000, snr=10,  cov_setting="group"))

simulation_settings <- rbind(simulation_settings,expand_grid(n=c(200), p=2000, act_setting=c("sparse","medium","dense"), 
                                                             ntest=1000, snr=c(1,5,20),  cov_setting="group"))

simulation_settings <- simulation_settings %>% mutate(a = round(ifelse(act_setting=="sparse",2*log(p),
                                                                       ifelse(act_setting=="medium",n/2+2*log(p),p/4))))

################## simulations ##############################################################################
unlink("../saved_results/log.txt")
# select number of cores eg detectCores()-1
my.cluster <- parallel::makeCluster(7, type = "PSOCK", outfile = "../saved_results/log.txt")
doParallel::registerDoParallel(cl = my.cluster)
# foreach::getDoParRegistered()

nmethods <- length(methods)
nset <- nrow(simulation_settings)
nmeas <- length(measures)
res <- array(c(0),dim=c(nrep,nmeas,nmethods,nset),dimnames = list(reps=NULL,measures=measures,
                                                              method=names(methods),settings=paste("Scenario",1:nset)))
attributes(res)$settings <- simulation_settings
fits <- vector(mode="list",length = nset)

clusterExport(my.cluster,c('simulation_settings','nmeas','nmethods','methods'), envir = environment())
clusterEvalQ(my.cluster, {  
  pacman::p_load(ROCR)
  source("../functions/multi_assign.R")
  source("../functions/data_generation.R")
  source("../functions/methods.R")
})
# # # #.# # - - - - - - START of AMP experiments - - - - - - - # # # # # # 
# # manually play around to check performance of BAMP (k=20)
tmp_res <- readRDS("../saved_results/sims_SPAR_24_nset45_reps100_nmeth17.rds")
j <- 2
j <- 16
i <- 1
parresi <- list(res=matrix(c(0),nmeas,nmethods),fits = vector(mode="list",length = nmethods))
names(parresi$fits) <- names(methods)

c(n,p,act_setting,ntest,snr,cov_setting,a) %<-% simulation_settings[j,]
set.seed((1234+i)^2)
data <- generate_data_linreg(n,p,cov_setting,ntest,snr=snr,a=a)
x <- data$x
y <- data$y
xtest <- data$xtest
ytest <- data$ytest
rMSPEconst <- mean((ytest-mean(y))^2)
rMSPEconst_tr <- var(y)

# for just BAMP only execute the following with
# k <- 20
for (k in 18:nmethods) {
  set.seed((1234+i)^2 + k)
  tstamp <- Sys.time()
  callres <- tryCatch( methods[[k]](x,y,xtest),
                       error=function(error_message) {
                         message(paste0("Error in ",names(methods)[k],": ",error_message))
                         return(NULL)
                       })
  tmp_time <- as.numeric(Sys.time() - tstamp,units="secs")
  
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
                         performance(prediction(as.numeric(abs(callres$beta)),as.numeric(data$beta!=0)),measure="auc",fpr.stop=n/2/(p-a))@y.values[[1]]) # pAUC
  }
  
  
  message(sprintf("Finished method %d / %d!",k,nmethods))
}
res[i,,,j] <- parresi$res
res[i,,1:17,j] <- tmp_res$res[i,,,j]
res[i,,,j]
cat(sprintf('Finished rep %d / %d for setting %d / %d at %s.\n',i,nrep,j,nset,Sys.time()))

# # # #.# # - - - - - - END of AMP experiments - - - - - - - # # # # # # 

parres <- foreach(j = 1:nset) %:%
foreach(i=1:nrep) %dopar% {
    parresi <- list(res=matrix(c(0),nmeas,nmethods),fits = vector(mode="list",length = nmethods))
    names(parresi$fits) <- names(methods)
  
    c(n,p,act_setting,ntest,snr,cov_setting,a) %<-% simulation_settings[j,]
    set.seed((1234+i)^2)
    data <- generate_data_linreg(n,p,cov_setting,ntest,snr=snr,a=a)
    x <- data$x
    y <- data$y
    xtest <- data$xtest
    ytest <- data$ytest
    rMSPEconst <- mean((ytest-mean(y))^2)
    rMSPEconst_tr <- var(y)
    # k <- 17
    for (k in 1:nmethods) {
      set.seed((1234+i)^2 + k)
      tstamp <- Sys.time()
      callres <- tryCatch( methods[[k]](x,y,xtest),
                                                error=function(error_message) {
                                                  message(paste0("Error in ",names(methods)[k],": ",error_message))
                                                  return(NULL)
                                                })
      tmp_time <- as.numeric(Sys.time() - tstamp,units="secs")
      
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
                             performance(prediction(as.numeric(abs(callres$beta)),as.numeric(data$beta!=0)),measure="auc",fpr.stop=n/2/(p-a))@y.values[[1]]) # pAUC
      }
      
      
      message(sprintf("Finished method %d / %d!",k,nmethods))
    }
    # res[i,,,j] <- parresi$res
    # res[i,,1:17,j] <- tmp_res$res[i,,,j]
    # res[i,,,j]
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
        file = sprintf("../saved_results/sims_SPAR_24_nset%d_reps%d_nmeth%d.rds",nset,nrep,nmethods))
parallel::stopCluster(cl = my.cluster)

warnings()
