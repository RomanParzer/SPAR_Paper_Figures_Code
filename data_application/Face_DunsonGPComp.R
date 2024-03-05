
# pacman::p_load(foreach,parallel,tidyr,dplyr,R.matlab)
# source("../functions/methods.R")
# facedata <- readMat("../data/face_data.mat")
# mydata <- data.frame(t(facedata$images+matrix(rnorm(698*64^2,0,0.03),4096,698)),y=facedata$poses[1,]) # maybe transform y
# mydata$y <- scale(mydata$y)
# 
# nrep <- 50
# res <- matrix(c(0),nrep,2)
# colnames(res) <- c("MSPE","se")
# 
# unlink("../saved_results/log.txt")
# # select nuber of cores
# my.cluster <- parallel::makeCluster(4, type = "PSOCK", outfile = "../saved_results/log.txt")
# doParallel::registerDoParallel(cl = my.cluster)
# 
# clusterExport(my.cluster,c('mydata'), envir = environment())
# clusterEvalQ(my.cluster, {  
#   source("../functions/methods.R")
# })
# 
# parres <- foreach(i=1:nrep) %dopar% {
#     parresi <- c(0,0)
#    
#     set.seed((1234+i)^2)
#     testind <- sample(1:698,50,replace=FALSE)
#     x <- as.matrix(mydata[-testind,-4097])
#     xtest <- as.matrix(mydata[testind,-4097])
#     y <- mydata$y[-testind]
#     ytest <- mydata$y[testind]
#     
#     
#     callres <- tryCatch( mySPAR(x,y,xtest,nummods=c(10,20,30,50,100),opt_par = "best"),
#                          error=function(error_message) {
#                            message(paste0("Error in ",names(methods)[k],": ",error_message))
#                            return(NULL)
#                          })
#     
#     parresi[1] <- mean((callres$yhat-ytest)^2) 
#     bootrs <- bootstrap::bootstrap((callres$yhat-ytest)^2,100,mean)
#     parresi[2] <- sd(bootrs$thetastar)
#     
#     cat(sprintf('Finished rep %d / %d at %s.\n',i,nrep,Sys.time()))
#     warnings()
#     parresi
#   }
# 
# for (i in 1:nrep) {
#   res[i,] <- parres[[i]]
# }
# 
# saveRDS(res, file = "../saved_results/Face_Dunson.rds")
# parallel::stopCluster(cl = my.cluster)

res_duns <- readRDS("../saved_results/Face_Dunson.rds")
colMeans(res_duns)

