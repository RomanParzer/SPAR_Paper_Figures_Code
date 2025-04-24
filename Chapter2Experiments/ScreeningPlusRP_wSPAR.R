
# # install.packages("remotes")
# remotes::install_github("RomanParzer/SPAR@v1.1.1")
pacman::p_load(foreach,parallel,tidyr,dplyr,robustHD,SPAR)
source("../functions/multi_assign.R")
source("../functions/data_generation.R")
source("../functions/RPM_generation.R")

# compare prediction performance of different random projections

mybeta <- numeric(2000)
set.seed(1234)
mybeta[1:100] <- sample(c(-3:3)[-4],100,replace = TRUE)
simulation_settings <- list(n=200, p=c(2000), act_setting="medium", ntest=100, snr=10,  cov_setting="comsym",m=100,beta=mybeta)

# nrep <- 10^2
# maxnum <- 50
# rMSPE <- array(NA,dim=c(nrep,5,maxnum),dimnames = list(reps=1:nrep,methods=c("Screening","RP_CW","ScrRP_CW","ScrRP","SPAR"),nummod=1:maxnum))
# for (i in 1:nrep) {
#   c(n,p,act_setting,ntest,snr,cov_setting,m,beta) %<-% simulation_settings
#   set.seed((1234+i)^2)
#   data <- generate_data_linreg(n,p,cov_setting,ntest,snr=snr,m=m,beta=beta)
# 
#   x <- data$x
#   y <- data$y
#   z <- robustHD::standardize(x)
#   yz <- robustHD::standardize(y)
# 
#   ztest <- scale(data$xtest,center = attributes(z)$center, scale = attributes(z)$scale)
#   
#   eig <- eigen(tcrossprod(z),symmetric = TRUE)
#   if (sum(eig$values>1e-8) >= (n-1)) {
#     myinv <- tcrossprod(eig$vectors[,eig$values>1e-8]%*%diag(1/sqrt(eig$values[eig$values>1e-8])))
#     solve_res <- myinv%*%yz
#   } else {
#     solve_res <- solve(tcrossprod(z)+(sqrt(p)+sqrt(n))*diag(n),yz)
#   }
#   HOLP <- as.numeric(crossprod(z,solve_res))
#   preds <- matrix(c(0),ntest,3)
#   ens_preds <- array(c(0),dim=c(ntest,maxnum,4))
# 
#   ## Screening
#   ind_use <- which(abs(HOLP) > quantile(abs(HOLP),probs = (p-n/2)/p,type=1))
#   zz <- z[,ind_use]
#   zztest <- ztest[,ind_use]
#   betahat <- tryCatch( solve(crossprod(zz),crossprod(zz,yz)),
#                        error=function(error_message) {
#                          return(solve(crossprod(zz)+0.01*diag(ncol(zz)),crossprod(zz,yz)))
#                        })
#   preds[,1] <- attributes(yz)$center + attributes(yz)$scale*zztest%*%betahat
# 
#   ## Sc+RP_CW
#   ind_use <- which(abs(HOLP) > quantile(abs(HOLP),probs = (p-2*n)/p,type=1))
#   RP_CW <- generate_RPM(m,2*n)
#   RPM <-  RP_CW
#   zz <- tcrossprod(z[,ind_use],RPM)
#   zztest <- tcrossprod(ztest[,ind_use],RPM)
#   betahat <- tryCatch( solve(crossprod(zz),crossprod(zz,yz)),
#                        error=function(error_message) {
#                          return(solve(crossprod(zz)+0.01*diag(ncol(zz)),crossprod(zz,yz)))
#                        })
#   preds[,2] <- as.numeric(attributes(yz)$center + attributes(yz)$scale*zztest%*%betahat)
# 
#   ## Sc+RP
#   RPM@x <- HOLP[ind_use]
#   zz <- tcrossprod(z[,ind_use],RPM)
#   zztest <- tcrossprod(ztest[,ind_use],RPM)
#   betahat <- tryCatch( solve(crossprod(zz),crossprod(zz,yz)),
#                        error=function(error_message) {
#                          return(solve(crossprod(zz)+0.01*diag(ncol(zz)),crossprod(zz,yz)))
#                        })
#   preds[,3] <- as.numeric(attributes(yz)$center + attributes(yz)$scale*zztest%*%betahat)
# 
#   rMSPE[i,c(1,3,4),1] <- apply(preds,2,function(yhat) mean((yhat-data$ytest)^2)  )
# 
#   ## Ensembles
#   ms <- sample(seq(ceiling(log(p)),ceiling(n/2)),maxnum,replace=TRUE)
#   for (k in 1:maxnum) {
#     tmp_m <- ms[k]
#     ## Screening
#     ind_use <- sample(1:p,n/2,prob=abs(HOLP))
#     zz <- z[,ind_use]
#     zztest <- ztest[,ind_use]
#     betahat <- tryCatch( solve(crossprod(zz),crossprod(zz,yz)),
#                          error=function(error_message) {
#                            return(solve(crossprod(zz)+0.01*diag(ncol(zz)),crossprod(zz,yz)))
#                          })
#     ens_preds[,k,1] <- attributes(yz)$center + attributes(yz)$scale*zztest%*%betahat
# 
#     ## RP
#     RP_CW <- generate_RPM(tmp_m,p)
#     zz <- tcrossprod(z,RP_CW)
#     zztest <- tcrossprod(ztest,RP_CW)
#     betahat <- tryCatch( solve(crossprod(zz),crossprod(zz,yz)),
#                          error=function(error_message) {
#                            return(solve(crossprod(zz)+0.01*diag(ncol(zz)),crossprod(zz,yz)))
#                          })
#     ens_preds[,k,2] <- as.numeric(attributes(yz)$center + attributes(yz)$scale*zztest%*%betahat)
# 
# 
#     ## Sc+RP_CW
#     ind_use <- sample(1:p,2*n,prob=abs(HOLP))
#     RP_CW <- generate_RPM(tmp_m,2*n)
#     RPM <- RP_CW
#     zz <- tcrossprod(z[,ind_use],RPM)
#     zztest <- tcrossprod(ztest[,ind_use],RPM)
#     betahat <- tryCatch( solve(crossprod(zz),crossprod(zz,yz)),
#                          error=function(error_message) {
#                            return(solve(crossprod(zz)+0.01*diag(ncol(zz)),crossprod(zz,yz)))
#                          })
#     ens_preds[,k,3] <- as.numeric(attributes(yz)$center + attributes(yz)$scale*zztest%*%betahat)
# 
#     # Sc+RP
#     RPM@x <- HOLP[ind_use]
#     zz <- tcrossprod(z[,ind_use],RPM)
#     zztest <- tcrossprod(ztest[,ind_use],RPM)
#     betahat <- tryCatch( solve(crossprod(zz),crossprod(zz,yz)),
#                          error=function(error_message) {
#                            return(solve(crossprod(zz)+0.01*diag(ncol(zz)),crossprod(zz,yz)))
#                          })
#     ens_preds[,k,4] <- as.numeric(attributes(yz)$center + attributes(yz)$scale*zztest%*%betahat)
#   }
# 
#   rMSPE[i,1,-1] <- sapply(2:maxnum,function(nummod) mean((rowMeans(ens_preds[,1:nummod,1,drop=FALSE])-data$ytest)^2)  )
#   rMSPE[i,2,] <- sapply(1:maxnum,function(nummod) mean((rowMeans(ens_preds[,1:nummod,2,drop=FALSE])-data$ytest)^2)  )
#   rMSPE[i,3,-1] <- sapply(2:maxnum,function(nummod) mean((rowMeans(ens_preds[,1:nummod,3,drop=FALSE])-data$ytest)^2)  )
#   rMSPE[i,4,-1] <- sapply(2:maxnum,function(nummod) mean((rowMeans(ens_preds[,1:nummod,4,drop=FALSE])-data$ytest)^2)  )
#   
#   # SPAR
#   spar_res <- spar(x,y,nummods=1:maxnum)
#   rMSPE[i,5,] <- sapply(1:maxnum,function(tmp_num) mean((predict(spar_res,data$xtest,nummod = tmp_num)-data$ytest)^2)  )
#   
#   if (i%%10==0) {
#     message(sprintf("Finished rep %d / %d at %s.",i,nrep,Sys.time()))
#   }
# }
# saveRDS(rMSPE,"./Result_ScrRP_ens_SPAR.rds")
rMSPE <- readRDS("Result_ScrRP_ens_SPAR.rds")

nrep <- dim(rMSPE)[1]
df_rMSPE <- reshape2::melt(rMSPE,varnames=names(dimnames(rMSPE))) %>% rename(MSPE=value,Method=methods,`Number of Models`=nummod)
sum_df_rMSPE <- group_by(df_rMSPE,Method,`Number of Models`) %>% summarize(sd=sd(MSPE)/sqrt(nrep),MSPE=mean(MSPE))
levels(sum_df_rMSPE$Method ) <- c("Scr_HOLP","RP_CW","ScrRP_CW","ScrRP","SPAR")
sum_df_rMSPE %>% 
  filter(Method!="SPAR") %>%
  ggplot(aes(x=`Number of Models`,y=MSPE,col=Method,linetype=Method)) +
  # geom_ribbon(aes(ymin=MSPE-sd,ymax=MSPE+sd),alpha=0.1) +
  # coord_cartesian(ylim=c(315,360)) + 
  geom_line(size=1)
# ggsave("../plots/ScRP_comsym_medium.pdf", height = 4, width = 8*0.8)
# ggsave("../plots/ScRP_comsym_medium_wSPAR.pdf", height = 5, width = 8)
# ggsave("../plots/ScRP_comsym_medium_wSPAR_woCI.pdf", height = 5, width = 8)

