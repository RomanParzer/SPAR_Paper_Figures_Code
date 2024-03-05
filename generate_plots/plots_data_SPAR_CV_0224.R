## read in file and plot  -----------------------------------------------------------------------
## ------------------------------------------------------------------------------------------------------------

pacman::p_load(dplyr, ggplot2, tidyr, ggrepel)

resobj <- readRDS("../saved_results/data_SPAR_24_nset3_reps100_nmeth16.rds")

res <- resobj$res
methods <- dimnames(res)[[3]]
dataset_sizes <- resobj$dataset_sizes
dataset_sizes$p[3] <- 3890 # delete "constant" pixels with 0 variance, 64x64 = 4096

show_methods <- methods[c(3,6:8,11:15)]
sparse_methods <- methods[c(6,7,8,14,15)]


mydf_all <- data.frame(pivot_longer(as.data.frame(res[,1,,1]),1:(dim(res)[3]),names_to="Method",values_to="rMSPE"),
                       rMSPE_tr=pivot_longer(as.data.frame(res[,2,,1]),1:(dim(res)[3]),names_to="Method",values_to="rMSPE_tr")$rMSPE_tr,
                       Time=pivot_longer(as.data.frame(res[,3,,1]),1:(dim(res)[3]),names_to="Method",values_to="Time")$Time,
                       numAct=pivot_longer(as.data.frame(res[,4,,1]),1:(dim(res)[3]),names_to="Method",values_to="numAct")$numAct,
                       dataset_sizes[1,])

for (k in 2:nrow(dataset_sizes)) {
  mydf_all <- rbind(mydf_all,
                    data.frame(pivot_longer(as.data.frame(res[,1,,k]),1:(dim(res)[3]),names_to="Method",values_to="rMSPE"),
                               rMSPE_tr=pivot_longer(as.data.frame(res[,2,,k]),1:(dim(res)[3]),names_to="Method",values_to="rMSPE_tr")$rMSPE_tr,
                               Time=pivot_longer(as.data.frame(res[,3,,k]),1:(dim(res)[3]),names_to="Method",values_to="Time")$Time,
                               numAct=pivot_longer(as.data.frame(res[,4,,k]),1:(dim(res)[3]),names_to="Method",values_to="numAct")$numAct,
                               dataset_sizes[k,]))
}

mydf_all$Method <- factor(mydf_all$Method,levels = methods[c(1,2,5:16,3,4)])
mydf_all <- mutate(mydf_all,"isSparse"=ifelse(Method%in%c("AdLASSO","ElNet","SIS"),TRUE,FALSE))

med_df <- mydf_all %>% group_by(Method,dataset) %>% summarize(MedrMSPE=median(rMSPE),MednumAct=median(numAct)) %>% filter(dataset=="face")
med_df$MednumAct[9] <- 3890

## face application

mydf_all %>% filter(Method %in% show_methods,dataset=="face",) %>%
  ggplot(aes(x=Method,y=rMSPE,fill=Method,linetype=isSparse)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  geom_segment(aes(x = c("AdLASSO"), y = c(0.025), xend = c("AdLASSO"), yend = c(0.03)),
               arrow = arrow(length = unit(0.5, "cm")),col=scales::hue_pal()(10)[1],linetype=2) +
  geom_segment(aes(x = c("SIS"), y = c(0.025), xend = c("SIS"), yend = c(0.03)),
               arrow = arrow(length = unit(0.5, "cm")),col=scales::hue_pal()(10)[3],linetype=2) +
  geom_segment(aes(x = c("RP_CW"), y = c(0.025), xend = c("RP_CW"), yend = c(0.03)),
               arrow = arrow(length = unit(0.5, "cm")),col=scales::hue_pal()(10)[6],linetype=2) +
  annotate("text",x="AdLASSO",y=0.0225,label=paste0("median\n ",round(med_df$MedrMSPE[med_df$Method=="AdLASSO"],3))) +
  annotate("text",x="SIS",y=0.0225,label=paste0("median\n ",round(med_df$MedrMSPE[med_df$Method=="SIS"],3))) +
  annotate("text",x="RP_CW",y=0.0225,label=paste0("median\n ",round(med_df$MedrMSPE[med_df$Method=="RP_CW"],3))) +
  coord_cartesian(ylim=c(0,0.03)) +
  theme(legend.position = "none")
# ggsave(paste0("../plots/SPAR_CV_rMSPE_face.pdf"), height = 4, width = 8)

# rateye application

mydf_all %>% filter(Method %in% c(show_methods,"Ridge"),dataset!="face",Method!="PLS") %>%
  ggplot(aes(x=Method,y=rMSPE,fill=dataset,linetype=isSparse)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  coord_cartesian(ylim=c(0,1.5)) +
  scale_linetype(guide="none")
# ggsave(paste0("../plots/SPAR_CV_rMSPE_rateye.pdf"), height = 4, width = 8)


med_df_all <- mydf_all %>% group_by(Method,dataset) %>% summarize(MedrMSPE=median(rMSPE),MednumAct=median(numAct))
kab_numAct_all <- as.matrix(pivot_wider(med_df_all[med_df_all$Method%in%c(show_methods,"Ridge"),-3],names_from = dataset,values_from = MednumAct))
# row.names(kab_rMSPE) <- tab_rMSPE$method
require(knitr)
# copy to latex
kable(kab_numAct_all[,c(1,3,4,2)],format = "latex",booktabs=TRUE)

# face application close

require(R.matlab)
facedata <- readMat("../data/face_data.mat")
myfacedata <- data.frame(t(facedata$images),y=facedata$poses[1,])

# select observation i = 1, ... , 698
i <- 1
plot1 <- ggplot(data.frame(X=rep(1:64,each=64),Y=rep(64:1,64),Z=facedata$images[,i]), aes(X, Y, fill= Z)) +
  geom_tile() +
  theme_void() +
  ggtitle(paste0("y=",round(facedata$poses[1,i],1))) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))
i <- 28
plot2 <- ggplot(data.frame(X=rep(1:64,each=64),Y=rep(64:1,64),Z=facedata$images[,i]), aes(X, Y, fill= Z)) +
  geom_tile() +
  theme_void() +
  ggtitle(paste0("y=",round(facedata$poses[1,i],1))) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))
i <- 106
plot3 <- ggplot(data.frame(X=rep(1:64,each=64),Y=rep(64:1,64),Z=facedata$images[,i]), aes(X, Y, fill= Z)) +
  geom_tile() +
  theme_void() +
  ggtitle(paste0("y=",round(facedata$poses[1,i],1))) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))
require(gridExtra)
grid.arrange(plot1, plot2,plot3, ncol=3)
# ggsave("../plots/faces.pdf",arrangeGrob(plot1, plot2,plot3,ncol=3),height = 4,width=12)

# remove two test instances
x <- myfacedata[-c(9,12),-4097]
y <- myfacedata[-c(9,12),4097]
# summary(apply(x,2,sd))
const_col_ind <- which(apply(x,2,sd)<0.01)
if (length(const_col_ind)>0) {
  x <- x[,-const_col_ind]
}
xtest <- myfacedata[c(9,12),-c(const_col_ind,4097)]
ytest <- myfacedata[c(9,12),4097]

require(SPAR)
set.seed(1234)
face_res <- spar.cv(x,y,nummods = c(10,20,30,50,100))
face_coef <- coef(face_res,opt_par = "best")
abs_coef <- numeric(4096)
abs_coef[-const_col_ind] <- abs(face_coef$beta)

coef <- numeric(4096)
coef[-const_col_ind] <- face_coef$beta
# coef <- log(abs(coef)+1)* sign(coef)
 
plot1 <- ggplot(data.frame(X=rep(1:64,each=64),Y=rep(64:1,64),`pos.coefs`=ifelse(coef>0,coef,0)), aes(X, Y, fill= `pos.coefs`)) +
  geom_tile() +
  theme_void() +
  scale_fill_gradient2()
plot1

plot2 <- ggplot(data.frame(X=rep(1:64,each=64),Y=rep(64:1,64),`neg.coefs`=ifelse(coef<0,coef,0)), aes(X, Y, fill= `neg.coefs`)) +
  geom_tile() +
  theme_void() +
  scale_fill_gradient2()
plot2

i <- 9
plot3 <- ggplot(data.frame(X=rep(1:64,each=64),Y=rep(64:1,64),Z=facedata$images[,i]), aes(X, Y, fill= Z)) +
  geom_tile() +
  theme_void() +
  ggtitle(paste0("y=",round(facedata$poses[1,i],1))) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))
plot3
shap_vals <- numeric(4096)
shap_vals[-const_col_ind] <- facedata$images[-const_col_ind,i]*face_coef$beta
# shap_vals[whichmax] <- max(shap_vals[-whichmax])
plot4 <- ggplot(data.frame(X=rep(1:64,each=64),Y=rep(64:1,64),effect=shap_vals), aes(X, Y, fill= effect)) +
  geom_tile() +
  theme_void() +
  scale_fill_gradient2() +
  ggtitle(bquote(hat(y) == .(round(predict(face_res,xnew=as.matrix(xtest),coef = face_coef)[1],1)))) +
  theme(plot.title = element_text(hjust = 0.5)) 
plot4
require(gridExtra)
grid.arrange(plot1, plot2,plot3,plot4, ncol=2)
# ggsave("../plots/faces_vis_result.pdf",arrangeGrob(plot1, plot2,plot3,plot4,ncol=2),height = 8,width=8)

