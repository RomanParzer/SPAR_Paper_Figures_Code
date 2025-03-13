#####################################################
### Testing sparse projected averaged regression SPAR
#####################################################

pacman::p_load(foreach,parallel,tidyr,dplyr,R.matlab)
source("../functions/multi_assign.R")
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
                "SLOPE"=mySLOPE,
                "BAMP" = function(x,y,xtest){myBAMP(x,y,xtest,gamma=0.8)})

measures <- c("rMSPE","rMSPE_tr","Time","numAct")

rateyedata <- data.frame(t(read.csv("../data/RatEyeAll.txt",sep="\t", header=TRUE, dec=".",row.names = 1)))
myq <- quantile(as.matrix(rateyedata),c(0.25))
out_q <- apply(rateyedata,2,function(mycol) max(mycol)<myq)
out_v <- apply(rateyedata,2, function(x) sd(x) / mean(x) * 100)
out_ind <- which((out_v < 2) | out_q)
# remove probes not sufficiently expressed and controls
mydata1 <- data.frame(rateyedata[,-unique(c(21231,31043:31099,out_ind))],y=rateyedata[,21231])

# require(flare)
# data(eyedata)
# tmp_data <- data.frame(X = x, y = y)
# write.csv(tmp_data,"../data/RatEye200.txt")
mydata2 <- read.csv("../data/RatEye200.txt",row.names = 1)

facedata <- readMat("../data/face_data.mat")
mydata3 <- data.frame(t(facedata$images),y=facedata$poses[1,]) 

datasets <- list(rateye=mydata1,
                 rateye200=mydata2,
                 face=mydata3)
dataset_sizes <- tibble(dataset=names(datasets),
                            n = sapply(datasets,nrow),
                            p = sapply(datasets,ncol)-1)

################## simulations ##############################################################################
unlink("../saved_results/log.txt")
my.cluster <- parallel::makeCluster(7, type = "PSOCK", outfile = "../saved_results/log.txt")
doParallel::registerDoParallel(cl = my.cluster)
# foreach::getDoParRegistered()

nmethods <- length(methods)
nset <- length(datasets)
nmeas <- length(measures)
res <- array(c(0),dim=c(nrep,nmeas,nmethods,nset),dimnames = list(reps=NULL,measures=measures,
                                                              method=names(methods),datasets=names(datasets)))
fits <- vector(mode="list",length = nset)

clusterExport(my.cluster,c('datasets','nmeas','nmethods','methods'), envir = environment())
clusterEvalQ(my.cluster, {  
  source("../functions/multi_assign.R")
  source("../functions/methods.R")
})

parres <- foreach(j = 1:(nset)) %:%
foreach(i=1:nrep) %dopar% {
  parresi <- list(res=matrix(c(0),nmeas,nmethods),fits = vector(mode="list",length = nmethods))
  names(parresi$fits) <- names(methods)
  
    mydata <- datasets[[j]]
    ntot <- nrow(mydata)
    ntest <- ntot%/%4
    n <- ntot - ntest
    p <- ncol(mydata) - 1
    
    set.seed((1234+i)^2)
    testind <- sample(1:ntot,ntest,replace=FALSE)
    x <- as.matrix(mydata[-testind,-(p+1)])
    # remove cols with sd 0
    const_col_ind <- which(apply(x,2,sd)==0)
    if (length(const_col_ind)>0) {
      x <- x[,-const_col_ind]
    }
    xtest <- as.matrix(mydata[testind,-c(const_col_ind,p+1)])
    # xtest <- as.matrix(mydata[testind,-c(p+1)])
    y <- mydata$y[-testind]
    ytest <- mydata$y[testind]
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
    
      if (is.null(callres$beta)) {
        parresi$res[,k] <- c(mean((callres$yhat-ytest)^2)/rMSPEconst, # rMSPE
                             mean((callres$yhat_tr-y)^2)/rMSPEconst_tr, # rMSPE_tr
                             tmp_time, # time
                             NA) # numAct
      } else {
        if ("dgCMatrix" %in% class(callres$beta)) {
          tmp_ind <- callres$beta@i+1
        } else {
          tmp_ind <- which(callres$beta!=0)
        }
        parresi$res[,k] <- c(mean((callres$yhat-ytest)^2)/rMSPEconst, # rMSPE
                             mean((callres$yhat_tr-y)^2)/rMSPEconst_tr, # rMSPE_tr
                             tmp_time, # time
                             length(tmp_ind)) # numAct
      }
      parresi$res[,k]
      # message(sprintf("Finished method %d / %d!",k,nmethods))
    }
    # res[i,,,j] <- parresi$res
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

saveRDS(list(res=res,fits=fits,dataset_sizes=dataset_sizes), 
        file = sprintf("../saved_results/data_SPAR_25_nset%d_reps%d_nmeth%d.rds",nset,nrep,nmethods))

parallel::stopCluster(cl = my.cluster)

warnings()


