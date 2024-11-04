###########################################
### benchmarks and competitors
##########################################

### overview: wrappers for functions providing predictions and set of active predictors (where applicable)
# 1 OLS (min Norm), 
# 2 PCR, 
# 3 PLS , 
# 4 Ridge, 
# 5 AdLASSO, 
# 6 Elastic Net with alpha=3/4, 
# 7 SIS and ISIS,
# 8 Stepwise forward BIC regression,
# 9 TARP,
# 10 SplitReg,
# 11 NCT / GCT,
# 12 RF?
# 13 SPAR wrapper
# 14 RP wrapper
# 15 HOLP Screening Wrapper

# # install.packages("remotes")
# remotes::install_github("RomanParzer/SPAR@v1.1.1")


pacman::p_load(pls, glmnet, SIS, MASS, SplitReg, robustHD, stringr,SPAR,Matrix)
source("../TARP-master/TARP.R")
source("../functions/RPM_generation.R")

# # 1 OLS (min Norm)
myOLS <- function(x,y,xtest) {
  p <- ncol(x)
  n <- nrow(x)
  z <- robustHD::standardize(x)
  yz <- robustHD::standardize(y)
  
  if (p < n) {
    betaz <- tryCatch( solve(crossprod(z),crossprod(z,yz)),
                               error=function(error_message) {
                                 return(solve(crossprod(z)+(sqrt(p)+sqrt(n))*diag(p),crossprod(z,yz)))
                               })
  } else {
    eig <- eigen(tcrossprod(z),symmetric = TRUE)
    if (sum(eig$values>1e-8) >= (n-1)) {
      myinv <- tcrossprod(eig$vectors[,eig$values>1e-8]%*%diag(1/sqrt(eig$values[eig$values>1e-8])))
      solve_res <- myinv%*%yz
    } else {
      solve_res <- solve(tcrossprod(z)+(sqrt(p)+sqrt(n))*diag(n),yz)
    }
    betaz <- crossprod(z,solve_res)
  }
  
  intercept <- as.numeric(attributes(yz)$center - attributes(yz)$scale*(attributes(z)$center/attributes(z)$scale)%*%betaz)
  beta <- attributes(yz)$scale* (betaz / attributes(z)$scale)
  
  yhat <- intercept + xtest%*%beta
  yhat_tr <- intercept + x%*%beta
  
  return(list(yhat=yhat, yhat_tr=yhat_tr, beta=beta,intercept=intercept))
}

# # 2 PCR
myPCR <- function(x,y,xtest) {
  pcrres <- pcr(y~.,data=data.frame(x,y=y),scale=TRUE,validation="CV",segments=10)
  nco <- which.min(pcrres$validation$PRESS)
  yhat <- predict(pcrres,ncomp = nco,newdata = data.frame(xtest))
  yhat_tr <- predict(pcrres,ncomp = nco,newdata = data.frame(x))
  coefs <- coef(pcrres,ncomp = nco,intercept = TRUE)
  
  return(list(yhat=yhat, yhat_tr=yhat_tr,nco=nco,beta=coefs[-1]/pcrres$scale,intercept=coefs[1]))
}

# # 3 PLS
myPLS <- function(x,y,xtest) {
  plsres <- plsr(y~.,data=data.frame(x,y=y),scale=TRUE,validation="CV",segments=10)
  nco <- which.min(plsres$validation$PRESS)
  yhat <- predict(plsres,ncomp = nco,newdata = data.frame(xtest))
  yhat_tr <- predict(plsres,ncomp = nco,newdata = data.frame(x))
  coefs <- coef(plsres,ncomp = nco,intercept = TRUE)
  
  return(list(yhat=yhat, yhat_tr=yhat_tr,nco=nco,beta=coefs[-1]/plsres$scale,intercept=coefs[1]))
}

# # 4 Ridge
# myRidge <- function(x,y,xtest) {
#   ridge_cv <- cv.glmnet(x,y,alpha=0,folds=10) # standardize = TRUE by default
#   # refit with optimal lambda
#   ridgeres <- glmnet(x,y,alpha=0, lambda = ridge_cv$lambda.1se)
#   
#   yhat <- predict(ridgeres,s = ridge_cv$lambda.1se,newx = xtest)
#   yhat_tr <- predict(ridgeres,s = ridge_cv$lambda.1se,newx = x)
# 
#   return(list(yhat=yhat, yhat_tr=yhat_tr,lambda=ridge_cv$lambda.1se,beta=ridgeres$beta,intercept=ridgeres$a0))
# }

# # 5 AdLASSO
myAdLASSO <- function(x,y,xtest) {
  ridge_cv <- cv.glmnet(x = x, y = y,nfold = 10,alpha = 0)
  best_ridge_coef <- as.numeric(coef(ridge_cv, s = ridge_cv$lambda.min))[-1]
  alasso_cv <- cv.glmnet(x=x,y=y,penalty.factor = 1 / abs(best_ridge_coef),alpha=1,nfold=10)
  # refit with optimal lambda
  alasso <- glmnet(x,y,penalty.factor = 1 / abs(best_ridge_coef),alpha=1, lambda = alasso_cv$lambda.1se)
  
  yhat <- predict(alasso,s = alasso_cv$lambda.1se,newx = xtest)
  yhat_tr <- predict(alasso,s = alasso_cv$lambda.1se,newx = x)
  
  return(list(yhat=yhat, yhat_tr=yhat_tr,lambda=alasso_cv$lambda.1se,beta=alasso$beta,intercept=alasso$a0))
}

# # 6 Elastic Net with alpha
myElNet <- function(x,y,xtest,alpha=3/4) {
  elnet_cv <- cv.glmnet(x,y,alpha=alpha,folds=10) # standardize = TRUE by default
  # refit with optimal lambda
  elnet <- glmnet(x,y,alpha=alpha, lambda = elnet_cv$lambda.1se)
  
  yhat <- predict(elnet,s = elnet_cv$lambda.1se,newx = xtest)
  yhat_tr <- predict(elnet,s = elnet_cv$lambda.1se,newx = x)
  
  return(list(yhat=yhat, yhat_tr=yhat_tr,lambda=elnet_cv$lambda.1se,beta=elnet$beta,intercept=elnet$a0))
}

# # 7 SIS and ISIS
mySIS <- function(x,y,xtest,iter=TRUE) {
  invisible(capture.output(
    SIS_res <- SIS(x,y,tune="cv",nfolds=10,iter=iter)
  ))
  
  yhat <- predict(SIS_res,newx=xtest,type="response")
  yhat_tr <- predict(SIS_res,newx=x,type="response")
  
  beta <- Matrix(c(0),nrow=ncol(x),ncol=1,sparse=TRUE)
  beta[SIS_res$ix,1] <- SIS_res$coef.est[-1]
  
  return(list(yhat=yhat, yhat_tr=yhat_tr,lambda=SIS_res$lambda,ind_screen=SIS_res$ix0,
              beta=beta,intercept=SIS_res$coef.est[1]))
}

# # 8 Stepwise forward BIC regression # takes quite long
myStepBIC <- function(x,y,xtest) {
  n <- nrow(x)
  p <- ncol(x)
  step_res <- stepAIC(lm(y~1,data=data.frame(x,y=y)),scope=paste0("~",paste0("X",1:(p-1),collapse = "",sep="+"),paste0("X",p)),
                      direction="forward",k=log(n),trace=0,steps=min(3*n/4,p))
  
  yhat <- predict(step_res,newdata = data.frame(xtest))
  yhat_tr <- step_res$fitted.values
  
  beta <- Matrix(c(0),nrow=p,ncol=1,sparse=TRUE)
  act_ind <- as.numeric(stringr::str_replace(names(step_res$coefficients[-1]),"X",""))
  beta[act_ind,1] <- step_res$coefficients[-1]

  return(list(yhat=yhat, yhat_tr=yhat_tr,intercept=step_res$coefficients[1],beta=beta))
}

# # 9 TARP (with different versions?) # some insight to beta, prec, recall?
myTARP <- function(x,y,xtest) {
  n <- nrow(x)
  z <- robustHD::standardize(x)
  ztest <- scale(xtest,center = attributes(z)$center, scale = attributes(z)$scale)
  yz <- robustHD::standardize(y)
  
  TARPres <- RIS_RP(z,yz,rbind(z,ztest),alpha=0.95)
  
  beta <- attributes(yz)$scale*TARPres$beta/attributes(z)$scale
  intercept <- attributes(yz)$center - sum(attributes(z)$center*beta)
  
  yhat <- xtest%*%beta + intercept
  yhat_tr <- x%*%beta + intercept
  return(list(yhat=yhat, yhat_tr=yhat_tr,beta=beta))
}

# # 10 SplitReg  (ensemble method), 
mySplitReg <- function(x,y,xtest,s=10) {
  splitreg_res <- SplitReg::cv.SplitReg(x,y,num_models = s, num_folds = 10)
  
  betas <- Matrix(splitreg_res$betas[,,splitreg_res$index_opt],sparse=TRUE)
  dimnames(betas) <- list(p=NULL,s=NULL)
  
  yhats <-  t(cbind(1,rbind(x,xtest))%*%rbind(splitreg_res$intercepts[1,,splitreg_res$index_opt],betas))
  coef <- predict(splitreg_res,type="coefficients")
  
  yhat <- predict(splitreg_res,newx = xtest)
  yhat_tr <- predict(splitreg_res,newx=x)
  
  return(list(yhat = yhat, yhat_tr = yhat_tr, 
              intercept=coef[1], beta=coef[-1],
              yhats = yhats, betas=betas))
}

# # 11 NCT / GCT,
myGCT <- function(x,y,xtest,phi=1,nfolds=10) {
  n <- nrow(x)
  p <- ncol(x)
  # center variables
  z <- scale(x,center = TRUE, scale = FALSE)
  yz <- scale(y,center = TRUE, scale = FALSE)
  
  svdx <- svd(z/sqrt(n),nu=0)
  A <- svdx$v%*%diag(svdx$d^(-1-phi)) 
  # A p x min (p,n)
  # b min(p,n) x 1
  
  # only cv b
  taus <- c(0)
  folds <- sample(cut(1:n,breaks = nfolds,labels = FALSE))
  bs <- matrix(c(0),min(n,p),nfolds)
  
  #create interesting threshold values tau
  for (k in 1:nfolds) {
    fold_ind <- which(folds==k)
    bs[,k] <- svdx$d^(-1+phi)*crossprod(svdx$v,crossprod(z[-fold_ind,],yz[-fold_ind])/(n-length(fold_ind)))
    taus <- c(taus,abs(bs[,k]))
  }
  taus <- quantile(taus,probs = 0:n/n)
  
  # calculate cv error for each tau and pick lowest
  cv_err <- numeric(length(taus))
  for (l in seq_along(taus)) {
    tau <- taus[l]
    cv_err[l] <- mean(sapply(1:nfolds,function(k){
      fold_ind <- which(folds==k)
      betak <- A%*%ifelse(abs(bs[,k])<tau,0,sign(bs[,k])*(abs(bs[,k])-tau))
      mean((yz[-fold_ind]-z[-fold_ind,]%*%betak)^2)
    }))
  }
  besttau <- taus[which.min(cv_err)]
  
  b <- svdx$d^(-1+phi)*crossprod(svdx$v,crossprod(z,yz)/n)
  beta <- A%*%ifelse(abs(b)<besttau,0,sign(b)*(abs(b)-besttau))
  intercept <- attributes(yz)$'scaled:center' - sum(attributes(z)$'scaled:center'*beta)
  yhat <- intercept + xtest%*%beta
  yhat_tr <- intercept + x%*%beta

  return(list(yhat=yhat, yhat_tr=yhat_tr, intercept=intercept,  beta=beta, tau=besttau, A=A, b=b))
}

# # 12 RF?, use ranger and tuneRanger (tunemtryfast)

# # 13 SPAR
mySPAR <- function(x,y,xtest,opt_par=c("best","1se"),nummods=c(20),nlambda=20) {
  if (nlambda==1 & length(nummods)==1) {
    spar_res <- spar(x,y,nummods=nummods,nlambda = nlambda)
    val_sum <- spar_res$val_res
    coef <- coef(spar_res)
  } else {
    spar_res <- spar.cv(x,y,nummods=nummods,nlambda = nlambda)
    val_sum <- spar_res$val_sum
    coef <- coef(spar_res,opt_par = opt_par)
  }
  
  yhat <- predict(spar_res,xnew = xtest,coef = coef)
  yhat_tr <- predict(spar_res,xnew = x,coef = coef)
  
  return(list(yhat=yhat, yhat_tr=yhat_tr,beta=coef$beta,intercept=coef$intercept,betas=spar_res$betas,val_sum=val_sum))
}

# # 14 RP wrapper
myRP <- function(x,y,xtest,RP=c("Gaussian","Sparse","CWSparse"),nummod=1) {
  
  p <- ncol(x)
  n <- nrow(x)
  
  RP <- match.arg(RP)
  
  ycenter <- mean(y)
  yscale <- sd(y)
  xcenter <- apply(x,2,mean)
  xscale <- apply(x,2,sd)
  xscale[xscale==0] <- 1
  
  z <- scale(x,center = xcenter,scale = xscale)
  yz <- scale(y,center = ycenter,scale = yscale)
  
  ms <- sample(seq(floor(2*log(p)),ceiling(n/2)),nummod,replace=TRUE)
  betas <- matrix(c(0),p,nummod)
  
  for (i in 1:nummod) {
    m <- ms[i]
    if (RP=="Gaussian") {
      Phi <- matrix(rnorm(p*m),m,p)
    } else if (RP=="Sparse") {
      Phi <- generate_RPmatrix(m,p,psi = 1/3)
    } else {
      Phi <- generate_RPmatrixCW(m,p)
    }
    zz <- tcrossprod(z,Phi)
    betas[,i] <- as.numeric(crossprod(Phi,solve(crossprod(zz),crossprod(zz,yz))))
  }
  
  beta <- yscale*rowMeans(betas)/xscale
  intercept <- ycenter - sum(xcenter*beta)
  
  yhat <- intercept + xtest%*%beta
  yhat_tr <- intercept + x%*%beta
  
  return(list(yhat=yhat, yhat_tr=yhat_tr,beta=beta,intercept=intercept))
}

# # 15 HOLP Screen
myHOLPScr <- function(x,y,xtest) {
  p <- ncol(x)
  n <- nrow(x)
  z <- robustHD::standardize(x)
  yz <- robustHD::standardize(y)
  
  if (p < n) {
    use_ind <- 1:p 
  } else {
    eig <- eigen(tcrossprod(z),symmetric = TRUE)
    if (sum(eig$values>1e-8) >= (n-1)) {
      myinv <- tcrossprod(eig$vectors[,eig$values>1e-8]%*%diag(1/sqrt(eig$values[eig$values>1e-8])))
      solve_res <- myinv%*%yz
    } else {
      solve_res <- solve(tcrossprod(z)+(sqrt(p)+sqrt(n))*diag(n),yz)
    }
    betaz <- crossprod(z,solve_res)
    use_ind <- which(abs(betaz)>=sort(abs(betaz),decreasing = TRUE)[n])
  }

  elnet_cv <- cv.glmnet(x[,use_ind],y,alpha=1,folds=10) # standardize = TRUE by default
  # refit with optimal lambda
  elnet <- glmnet(x[,use_ind],y,alpha=1, lambda = elnet_cv$lambda.1se)
  
  yhat <- predict(elnet,s = elnet_cv$lambda.1se,newx = xtest[,use_ind])
  yhat_tr <- predict(elnet,s = elnet_cv$lambda.1se,newx = x[,use_ind])
  
  mybeta <- numeric(p)
  mybeta[use_ind] <- elnet$beta
  return(list(yhat=yhat, yhat_tr=yhat_tr,lambda=elnet_cv$lambda.1se,beta=mybeta,intercept=elnet$a0))
}


# # for testing:
# source("./data_generation.R")
# p <- 2000
# ntest <- 100
# a <- 50
# n <- 200
# j <- 4
# 
# cov_setting <- c("ind","comsym","ar1","group","factor","extreme")[j]
# set.seed(1604879)
# dat <- generate_data_linreg(n,p,cov_setting,ntest,snr=10,a=a)
# x <- dat$x
# xtest <- dat$xtest
# y <- dat$y
# ytest <- dat$ytest
# 
# res <- myPLS(x,y,xtest)
# res <- myAdLASSO(x,y,xtest)
# res <- mySIS(x,y,xtest)
# res <- myTARP(x,y,xtest)
# 
# mean(which(res$beta!=0) %in% dat$ind)
# mean(dat$ind %in% which(res$beta!=0))
# 
# yhat <- res$yhat
# mean((ytest-yhat)^2) / mean((ytest-mean(y))^2)

