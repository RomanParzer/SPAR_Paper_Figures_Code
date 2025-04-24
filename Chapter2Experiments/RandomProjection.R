
pacman::p_load(foreach,parallel,tidyr,dplyr,robustHD)
source("../functions/multi_assign.R")
source("../functions/data_generation.R")
source("../functions/RPM_generation.R")

# compare prediction performance of different random projections

mybeta <- numeric(2000)
set.seed(1234)
mybeta[1:100] <- sample(c(-3:3)[-4],100,replace = TRUE)
simulation_settings <- list(n=200, p=c(2000), act_setting="medium", ntest=100, snr=10,  cov_setting="comsym",a=100,beta=mybeta)
c(n,p,act_setting,ntest,snr,cov_setting,a,beta) %<-% simulation_settings
m <- a

# nrep <- 10^2
# rMSPE <- matrix(NA,nrep,8)
# colnames(rMSPE) <- c("HOLP","Gaussian","Sparse","SparseCW","SparseCWSignH","SparseCWHolp","SparseCWSignB","SparseCWBeta")
# 
# for (i in 1:nrep) {
#   set.seed((1234+i)^2)
#   data <- generate_data_linreg(n,p,cov_setting,ntest,snr=snr,a=a,beta=beta)
# 
#   m <- a
#   x <- data$x
#   y <- data$y
#   z <- robustHD::standardize(x)
#   yz <- robustHD::standardize(y)
#   ztest <- scale(data$xtest,center = attributes(z)$center, scale = attributes(z)$scale)
# 
#   RP_CW <- generate_RPmatrixCW(m,p)
#   eig <- eigen(tcrossprod(z),symmetric = TRUE)
#   if (sum(eig$values>1e-8) >= (n-1)) {
#     myinv <- tcrossprod(eig$vectors[,eig$values>1e-8]%*%diag(1/sqrt(eig$values[eig$values>1e-8])))
#     solve_res <- myinv%*%yz
#   } else {
#     solve_res <- solve(tcrossprod(z)+(sqrt(p)+sqrt(n))*diag(n),yz)
#   }
#   HOLP <- as.numeric(crossprod(z,solve_res))
# 
#   RP_CW_SignH <- RP_CW_HOLP  <- RP_CW_SignB <- RP_CW_Beta <- RP_CW
#   RP_CW_SignH@x <- sign(HOLP)
#   RP_CW_HOLP@x <- HOLP
# 
#   RP_CW_SignB@x <- sign(data$beta)
#   # RP_CW_SignB@x <- ifelse(RP_CW_SignB@x==0,sample(c(-1,1),p,replace=TRUE),RP_CW_SignB@x)
#   goal_dims <- RP_CW_SignB@i +1
#   for (j in (1:m)[-unique(goal_dims[data$beta!=0])]) {
#     RP_CW_SignB@x[min(which(goal_dims==j))] <- sample(c(-1,1),1)
#   }
# 
#   RP_CW_Beta@x <- data$beta
# 
#   goal_dims <- RP_CW_Beta@i +1
#   for (j in (1:m)[-unique(goal_dims[data$beta!=0])]) {
#     RP_CW_Beta@x[min(which(goal_dims==j))] <- sample(c(-1,1),1)*min(abs(data$beta[data$beta!=0]))
#   }
# 
# 
#   reductions <- list("Gaussian"=matrix(rnorm(p*m),m,p),
#                      "Sparse"=generate_RPmatrix(m,p,psi = 1/3),
#                      "SparseCW"=RP_CW,
#                      "SparseCWSignH"=RP_CW_SignH,
#                      "SparseCWHolp"=RP_CW_HOLP,
#                      "SparseCWSignB"=RP_CW_SignB,
#                      "SparseCWBeta"=RP_CW_Beta)
# 
#   mspes <- lapply(reductions, function(Phi){
#     zz <- tcrossprod(z,Phi)
#     zztest <- tcrossprod(ztest,Phi)
#     if (m<n) {
#       betahat <- solve(crossprod(zz),crossprod(zz,yz))
#     } else {
#       betahat <- crossprod(zz,solve(tcrossprod(zz),yz,tol=1e-30))
#     }
#     pred <- attributes(yz)$center + attributes(yz)$scale*zztest%*%betahat
#     mean((pred-data$ytest)^2) #/ mean((mean(y)-data$ytest)^2)
#     })
# 
#   mspeHOLP <- mean((attributes(yz)$center + attributes(yz)$scale*ztest%*%HOLP-data$ytest)^2)# / mean((mean(y)-data$ytest)^2)
# 
#   rMSPE[i,] <- c(mspeHOLP,unlist(mspes))
# 
#   if (i%%10==0) {
#     message(sprintf("Finished rep %d / %d at %s.",i,nrep,Sys.time()))
#   }
# }
# saveRDS(rMSPE,"./Result_RP_final_sparse_SignB.rds")
rMSPE <- readRDS("Result_RP_final_sparse_SignB.rds")

df_rMSPE <- pivot_longer(as.data.frame(rMSPE),1:8,names_to = "Projection",values_to = "MSPE")
df_rMSPE$Projection <- factor(df_rMSPE$Projection,levels = levels(factor(df_rMSPE$Projection))[c(1,3,4,8,6,2,7,5)])
# lam1 = 1-rho+p rho, lamj = 1-rho
lamp <- 0.5
tau <- 1
a <- 100
th_bound <- sum(beta^2)*(lamp*(1-2*m/p)) +a/(p-1) *m*lamp*tau^2*(1-(m+1)/(p-1))
summ_df <- df_rMSPE %>% group_by(Projection) %>% summarise(medMSPE=median(MSPE))
medPBeta <- summ_df$medMSPE[8]
tmp_lab <- expression('C'['Th1'])
ggplot(df_rMSPE %>% filter(Projection!="HOLP"),aes(x=Projection,y=MSPE,fill=Projection)) +
  geom_boxplot() +
  theme(legend.position = "none") +
  geom_hline(yintercept=medPBeta,linetype=2) +
  geom_hline(yintercept=medPBeta+th_bound,linetype=2) +
  geom_segment(mapping=aes(x="SparseCW", y=medPBeta, xend="SparseCW", yend=medPBeta+th_bound), arrow=arrow(ends='both')) +
  annotate("text",x=3-0.25,y=300,label=as.character(tmp_lab),parse=TRUE)
# ggsave("../plots/RP_comsym_medium_sparseSignB.pdf", height = 4, width = 8)
