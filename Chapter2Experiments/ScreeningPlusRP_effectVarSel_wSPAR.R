
# # install.packages("remotes")
# remotes::install_github("RomanParzer/SPAR@v1.1.1")
pacman::p_load(foreach,parallel,tidyr,dplyr,robustHD,SPAR)
source("../functions/multi_assign.R")
source("../functions/data_generation.R")
source("../functions/RPM_generation.R")


mybeta <- numeric(2000)
set.seed(1234)
mybeta[1:100] <- sample(c(-3:3)[-4],100,replace = TRUE)
simulation_settings <- list(n=200, p=c(2000), act_setting="medium", ntest=100, snr=10,  cov_setting="comsym",m=100,beta=mybeta)

# nrep <- 10^2
# nummod <- 20
# nscreens <- 2^(c(-2,-1,-1/2,-0.1,0.1,1/2,1,2,3))
# rMSPE <- array(NA,dim=c(nrep,3,length(nscreens)),dimnames = list(reps=1:nrep,methods=c("Screening","ScrRP","SPAR"),nscreens=nscreens))
# 
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
# 
# 
#   ## Ensembles
#   ms <- sample(seq(ceiling(log(p)),ceiling(n/2)),nummod,replace=TRUE)
# 
#   for (j in 1:length(nscreens)) {
#     nscreen <- round(n*nscreens[j])
#     ens_preds <- array(c(0),dim=c(ntest,nummod,2))
#     for (k in 1:nummod) {
#       tmp_m <- ms[k]
#       ## Screening
#       ind_use <- sample(1:p,nscreen,prob=abs(HOLP))
#       zz <- z[,ind_use]
#       zztest <- ztest[,ind_use]
#       if (nscreen < n) {
#         betahat <- tryCatch( solve(crossprod(zz),crossprod(zz,yz)),
#                              error=function(error_message) {
#                                return(solve(crossprod(zz)+0.01*diag(ncol(zz)),crossprod(zz,yz)))
#                              })
#       } else {
#         eig <- eigen(tcrossprod(zz),symmetric = TRUE)
#         if (sum(eig$values>1e-8) >= (n-1)) {
#           myinv <- tcrossprod(eig$vectors[,eig$values>1e-8]%*%diag(1/sqrt(eig$values[eig$values>1e-8])))
#           solve_res <- myinv%*%yz
#         } else {
#           solve_res <- solve(tcrossprod(zz)+0.01*diag(n),yz)
#         }
#         betahat <- as.numeric(crossprod(zz,solve_res))
#       }
# 
#       ens_preds[,k,1] <- attributes(yz)$center + attributes(yz)$scale*zztest%*%betahat
# 
#       ## Sc+RP
#       ind_use <- sample(1:p,nscreen,prob=abs(HOLP))
#       if (tmp_m < nscreen) {
#         # RP_CW <- SPAR:::generate_RPM(tmp_m,nscreen)
#         RP_CW <- generate_RPM(tmp_m,nscreen)
#         RPM <- RP_CW
#         RPM@x <- HOLP[ind_use]
#       } else {
#         RPM <- Matrix(diag(1,nscreen),sparse=TRUE)
#       }
# 
#       zz <- tcrossprod(z[,ind_use],RPM)
#       zztest <- tcrossprod(ztest[,ind_use],RPM)
#       betahat <- tryCatch( solve(crossprod(zz),crossprod(zz,yz)),
#                            error=function(error_message) {
#                              return(solve(crossprod(zz)+0.01*diag(ncol(zz)),crossprod(zz,yz)))
#                            })
#       ens_preds[,k,2] <- as.numeric(attributes(yz)$center + attributes(yz)$scale*zztest%*%betahat)
#     }
# 
#     rMSPE[i,1,j] <- mean((rowMeans(ens_preds[,,1,drop=FALSE])-data$ytest)^2)
#     rMSPE[i,2,j] <- mean((rowMeans(ens_preds[,,2,drop=FALSE])-data$ytest)^2)
#     
#     # SPAR
#     spar_res <- spar(x,y,nummods=c(nummod),nscreen = nscreen,msup = min(nscreen,n/2))
#     rMSPE[i,3,j] <- mean((predict(spar_res,data$xtest)-data$ytest)^2) 
#   }
# 
# 
#   if (i%%10==0) {
#     message(sprintf("Finished rep %d / %d at %s.",i,nrep,Sys.time()))
#   }
# }
# 
# saveRDS(rMSPE,"./Result_ScrRP_nscreen_SPAR.rds")

rMSPE <- readRDS("Result_ScrRP_nscreen_SPAR.rds")
nrep <- dim(rMSPE)[1]
df_rMSPE <- reshape2::melt(rMSPE) %>% rename(MSPE=value,Method=methods,`Ratio of n of screened variables`=nscreens)
sum_df_rMSPE <- group_by(df_rMSPE,Method,`Ratio of n of screened variables`) %>% summarize(sd=sd(MSPE)/sqrt(nrep),MSPE=mean(MSPE))
levels(sum_df_rMSPE$Method) <- c("Scr_HOLP","ScrRP","SPAR")

minScr <- min(sum_df_rMSPE$MSPE[sum_df_rMSPE$Method=="Scr_HOLP"])
sum_df_rMSPE %>%
  # filter(Method!="SPAR") %>%
  ggplot(aes(x=`Ratio of n of screened variables`,y=MSPE,col=Method)) +
  geom_line(size=1) +
  geom_ribbon(aes(ymin=MSPE-sd,ymax=MSPE+sd),alpha=0.3) +
  geom_hline(data=data.frame(yintercept = min(sum_df_rMSPE$MSPE[sum_df_rMSPE$Method=="Scr_HOLP"])),aes(yintercept=yintercept),linetype=2) +
  geom_text(data=data.frame(x= 0.8,
                            y = minScr,
                            Method="Scr_HOLP",
                            label=round(minScr,1)),
            aes(x=x,y=y,label=label),col=1,nudge_y = -9) +
  geom_text(data=data.frame(x= 0.8,
                            y = minScr,
                            Method="ScrRP",
                            label=round(minScr,1)),
            aes(x=x,y=y,label=label),col=1,nudge_y = -0.7) +
  geom_text(data=data.frame(x= 0.8,
                            y = minScr,
                            Method="SPAR",
                            label=round(minScr,1)),
            aes(x=x,y=y,label=label),col=1,nudge_y = -0.7) +
  scale_x_log10() +
  theme(legend.position = "none") +
  facet_grid(Method~., scales = "free_y") 
# ggsave("../plots/ScRP_nsc_comsym_medium_no1.pdf", height = 5, width = 8)
# ggsave("../plots/ScRP_nsc_comsym_medium_no1_wSPAR.pdf", height = 5, width = 8)
# ggsave("../plots/ScRP_nsc_comsym_medium_no1_wSPAR_woCI.pdf", height = 5, width = 8)


