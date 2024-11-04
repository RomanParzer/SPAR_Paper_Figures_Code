
pacman::p_load(foreach,parallel,tidyr,dplyr,robustHD)
source("../functions/multi_assign.R")
source("../functions/data_generation.R")


mybeta <- numeric(2000)
set.seed(1234)
mybeta[1:100] <- sample(c(-3:3)[-4],100,replace = TRUE)
simulation_settings <- list(n=200, p=c(2000), act_setting="medium", ntest=100, snr=10,  cov_setting="comsym",m=100,beta=mybeta)
c(n,p,act_setting,ntest,snr,cov_setting,m,beta) %<-% simulation_settings

# nrep <- 10^2
# PrecRec <- data.frame(Method=NULL,Number=NULL,Recall=NULL,Precision=NULL,Sign=NULL, Cor=NULL)
# hist_df <- data.frame(HOLP=NULL, Ridge=NULL, Ridge_cv=NULL, Cor=NULL,predictors=NULL)
# 
# for (i in 1:nrep) {
#   set.seed((1234+i)^2)
#   data <- generate_data_linreg(n,p,cov_setting,ntest,snr=snr,m=m,beta=beta)
# 
#   x <- data$x
#   y <- data$y
# 
#   z <- robustHD::standardize(x)
#   yz <- robustHD::standardize(y)
# 
#   eig <- eigen(tcrossprod(z),symmetric = TRUE)
#   if (sum(eig$values>1e-8) >= (n-1)) {
#     myinv <- tcrossprod(eig$vectors[,eig$values>1e-8]%*%diag(1/sqrt(eig$values[eig$values>1e-8])))
#     solve_res <- myinv%*%yz
#   } else {
#     solve_res <- solve(tcrossprod(z)+(sqrt(p)+sqrt(n))*diag(n),yz)
#   }
#   HOLP <- crossprod(z,solve_res)
#   rid2 <- crossprod(z,solve(tcrossprod(z)+(sqrt(p)+sqrt(n))*diag(n),yz))
#   cvridgeres <- glmnet::cv.glmnet(z,yz,standardize=FALSE,alpha=0,folds=10)
#   myl <- cvridgeres$lambda.1se
#   rid_cv <- glmnet::glmnet(z,yz,standardize=FALSE,alpha=0, lambda = myl)$beta@x
#   cor_y <- cov(z,yz)
# 
#   coefs <- cbind(HOLP,rid2,rid_cv,cor_y)
#   colnames(coefs) <- c("HOLP","Ridge","Ridge_cv","Cor")
# 
#   abs_coef <- abs(coefs)
#   hist_df <- rbind(hist_df,data.frame(abs_coef,predictors=ifelse(1:p %in% data$ind,"active","non-active")))
# 
#   for (k in 1:4) {
#     coef <- coefs[,k]
#     abscoef <- abs(coef)
#     sortcoef <- sort(abscoef,decreasing = TRUE)
#     method <- c("HOLP","Ridge","Ridge_cv","Cor")[k]
#     tabres <- sapply(unique(round(exp(seq(0,log(p),length.out=100)))), function(l){
#       ind2 <- which(abscoef >= sortcoef[l])
#       coefa <- coef[data$ind[data$ind %in% ind2]] - mean(coef[data$ind[data$ind %in% ind2]])
#       coefb <- data$beta[data$ind[data$ind %in% ind2]] - mean(data$beta[data$ind[data$ind %in% ind2]])
#       sdb <- sd(coefb)
#       if (!is.na(sdb)) {
#         if (sdb==0) {
#           sdb <- sd(data$beta)
#         }
#       }
# 
#       c(l,
#         mean(data$ind %in% ind2),
#         mean(ind2 %in% data$ind),
#         mean(sign(coef[data$ind[data$ind %in% ind2]])==sign(data$beta[data$ind[data$ind %in% ind2]])),
#         # cor(coef[data$ind[data$ind %in% ind2]],data$beta[data$ind[data$ind %in% ind2]],method="kendall")
#         sum(coefa*coefb) / ((sum(data$ind %in% ind2)-1)*sd(coefa)*sdb)
#         )
#     })
#     PrecRec <- rbind(PrecRec,data.frame(Method=method,t(tabres)))
#   }
# 
#   if (i%%10==0) {
#     message(sprintf("Finished rep %d / %d at %s.",i,nrep,Sys.time()))
#   }
# }
# 
# colnames(PrecRec) <- c("Method","Number","Recall","Precision","Sign","Cor")
# 
# saveRDS(list(PrecRec=PrecRec,hist_df=hist_df),"./Result_Screening.rds")
save_obj <- readRDS("Result_Screening.rds")
PrecRec <- save_obj$PrecRec
hist_df <- save_obj$hist_df

PrecRec_long <- PrecRec %>% pivot_longer(3:6,names_to = "Measure",values_to = "Value")
PrecRec_long$Measure <- factor(PrecRec_long$Measure,levels = c("Precision","Recall","Sign","Cor"))
PrecRec_long$Method <- factor(PrecRec_long$Method,levels = c("Cor","Ridge_cv","Ridge","HOLP"))

PrecRec_long %>%
  ggplot(aes(x=Number,y=Value,col=Method,linetype=Method)) +
  geom_smooth(size=0.7,se = TRUE, alpha=0.4) +
  facet_grid(Measure~.) +
  scale_x_log10() +
  labs(x="Number of selected variables",y="") +
  geom_vline(xintercept=100,col="black",size=0.5,alpha=0.5)

# ggsave("../plots/Screening_comsym_medium.pdf", height = 6, width = 6)

hist_df_long <- hist_df %>%
  pivot_longer(1:4,names_to = "Method",values_to = "abs.coefs")
hist_df_long$Method <- factor(hist_df_long$Method,levels = c("Cor","Ridge_cv","Ridge","HOLP"))

hist_df_long %>%
  ggplot( aes(x=abs.coefs, color=predictors,linetype=predictors)) +
  geom_density(size=1) +
  scale_fill_brewer(palette = "Set1") +
  facet_wrap(Method~.,scales = "free",nrow=4, strip.position="right") +
  labs(x="Absolute coefficients",y="Density",color="Predictors",linetype="Predictors")

# ggsave("../plots/Hist_comsym_medium_all.pdf", height = 6, width = 6)


