## read in file and plot  -----------------------------------------------------------------------
## ------------------------------------------------------------------------------------------------------------

pacman::p_load(dplyr, ggplot2, tidyr, ggrepel,knitr,kableExtra)

resobj <- readRDS("../saved_results/sims_SPAR_24_nset45_reps100_nmeth17.rds")
res <- resobj$res

methods <- dimnames(res)[[3]]
settings <- attributes(res)$settings

show_methods <- methods[c(1,3,6:8,11:15)]
sparse_methods <- methods[c(6,7,8,14,15)]


mydf_all <- data.frame(pivot_longer(data.frame(res[,1,,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="rMSPE"),
                       rMSPE_tr=pivot_longer(data.frame(res[,2,,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="rMSPE_tr")$rMSPE_tr,
                       Time=pivot_longer(data.frame(res[,3,,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="Time")$Time,
                       Precision=pivot_longer(data.frame(res[,4,,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="Precision")$Precision,
                       Recall=pivot_longer(data.frame(res[,5,,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="Recall")$Recall,
                       activeEE=pivot_longer(data.frame(res[,6,,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="activeEE")$activeEE,
                       passiveEE=pivot_longer(data.frame(res[,7,,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="passiveEE")$passiveEE,
                       numAct=pivot_longer(data.frame(res[,8,,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="numAct")$numAct,
                       AUC=pivot_longer(data.frame(res[,9,,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="AUC")$AUC,
                       pAUC=pivot_longer(data.frame(res[,10,,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="pAUC")$pAUC,
                       settings[1,],
                       setting=1)

for (k in 2:nrow(settings)) {
  mydf_all <- rbind(mydf_all,
                    data.frame(pivot_longer(data.frame(res[,1,,k],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="rMSPE"),
                               rMSPE_tr=pivot_longer(data.frame(res[,2,,k],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="rMSPE_tr")$rMSPE_tr,
                               Time=pivot_longer(data.frame(res[,3,,k],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="Time")$Time,
                               Precision=pivot_longer(data.frame(res[,4,,k],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="Precision")$Precision,
                               Recall=pivot_longer(data.frame(res[,5,,k],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="Recall")$Recall,
                               activeEE=pivot_longer(data.frame(res[,6,,k],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="activeEE")$activeEE,
                               passiveEE=pivot_longer(data.frame(res[,7,,k],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="passiveEE")$passiveEE,
                               numAct=pivot_longer(data.frame(res[,8,,k],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="numAct")$numAct,
                               AUC=pivot_longer(data.frame(res[,9,,k],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="AUC")$AUC,
                               pAUC=pivot_longer(data.frame(res[,10,,k],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="pAUC")$pAUC,
                               settings[k,],
                               setting=k))
}

# str(mydf_all)
mydf_all <- mutate(mydf_all,"F1"=ifelse(Precision==0&Recall==0,0,2*Precision*Recall/(Precision+Recall)))
mydf_all$act_setting <- factor(mydf_all$act_setting,levels = c("sparse","medium","dense"))
mydf_all$Method <- stringr::str_replace_all(mydf_all$Method,"\\."," ")
mydf_all$Method <- factor(mydf_all$Method,levels = methods[c(2,4:8,17,9:16,3,1)])
mydf_all$cov_setting <- factor(mydf_all$cov_setting,levels = c("ind","comsym","ar1","group","factor","extreme"))
mydf_all <- mutate(mydf_all,"isSparse"=ifelse(Method%in%c("LASSO","AdLASSO","ElNet","SIS","HOLPScr"),TRUE,FALSE))
# rescale pAUC to (0,1)
mydf_all <- mutate(mydf_all,"pAUC"=pAUC*2*(p-a)/n)

# plot medium p2000 the 6 cov settings rMSPE 
mydf_all %>% filter(Method %in% show_methods,p==2000,n==200,snr==10) %>%
  ggplot(aes(x=Method,y=rMSPE,fill=Method,linetype=isSparse)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  facet_grid(cov_setting~act_setting, scales = "free_y") +
  coord_cartesian(ylim=c(0,1.35))+
  theme(legend.position = "none") +
  geom_hline(aes(yintercept=1),linetype=2)
# ggsave(paste0("../plots/SPAR_CV_rMSPE_cov_settings.pdf"), height = 10, width = 8)
# ggsave(paste0("../plots/SPAR_CV_rMSPE_cov_settings_HOLPScr.pdf"), height = 10, width = 8)


# pAUC
mydf_all %>% filter(Method %in% show_methods,p==2000,n==200,snr==10) %>%
  ggplot(aes(x=Method,y=pAUC,fill=Method,linetype=isSparse)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  ggh4x::facet_grid2(cov_setting~act_setting, scales = "free_y",independent = "y") +
  theme(legend.position = "none")
# ggsave(paste0("../plots/SPAR_CV_pAUC_cov_settings.pdf"), height = 10, width = 8)

# AUC
mydf_all %>% filter(Method %in% show_methods,p==2000,n==200,snr==10) %>%
  ggplot(aes(x=Method,y=AUC,fill=Method,linetype=isSparse)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  ggh4x::facet_grid2(cov_setting~act_setting, scales = "free_y",independent = "y") +
  theme(legend.position = "none")
# ggsave(paste0("../plots/SPAR_CV_AUC_cov_settings.pdf"), height = 10, width = 8)

# Prec
mydf_all %>% filter(Method %in% c(sparse_methods,"TARP"),p==2000,n==200,snr==10) %>%
  ggplot(aes(x=Method,y=Precision,fill=Method,linetype=isSparse)) +
  geom_boxplot() +
  geom_hline(aes(yintercept=a/p),linetype=2) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  facet_grid(cov_setting~act_setting, scales = "free_y") +
  theme(legend.position = "none")
# ggsave(paste0("../plots/SPAR_CV_Precision_cov_settings.pdf"), height = 10, width = 8)

# Rec
mydf_all %>% filter(Method %in% c(sparse_methods,"TARP"),p==2000,n==200,snr==10) %>%
  ggplot(aes(x=Method,y=Recall,fill=Method,linetype=isSparse)) +
  geom_boxplot() +
  geom_hline(yintercept=1,linetype=2) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  facet_grid(cov_setting~act_setting, scales = "free_y") +
  theme(legend.position = "none")
# ggsave(paste0("../plots/SPAR_CV_Recall_cov_settings.pdf"), height = 10, width = 8)


# compute ranks

n_showm <- length(show_methods)
mydf_all_rank_rMSPE_pAUC <- mydf_all %>% filter(Method %in% show_methods
                                                ) %>% group_by(rep,setting) %>% mutate(rank_pAUC=n_showm + 1 - rank(pAUC) ,rank_rMSPE=rank(rMSPE))
rank_tab_rMSPE_pAUC <- mydf_all_rank_rMSPE_pAUC %>% group_by(Method) %>% summarise(mean_rank_rMSPE = mean(rank_rMSPE), se_rank_rMSPE = sd(rank_rMSPE)/sqrt(100*3*6),
                                                                                   mean_rank_pAUC = mean(rank_pAUC), se_rank_pAUC = sd(rank_pAUC)/sqrt(100*3*6)) 
rank_tab_rMSPE_pAUC
# kable(rank_tab_rMSPE_pAUC,digits=3,format="latex",booktabs=TRUE)

n_sparm <- length(sparse_methods)
mydf_all_rank_Prec_Rec <- mydf_all %>% filter(Method %in% sparse_methods, p==2000,n==200,snr==10) %>% group_by(rep,setting) %>% mutate(rank_Prec=n_sparm + 1 -rank(Precision), rank_Rec=n_sparm + 1 -rank(Recall))
rank_tab_Prec_Rec <- mydf_all_rank_Prec_Rec %>% group_by(Method) %>% summarise(mean_rank_Prec = mean(rank_Prec), se_rank_Prec = sd(rank_Prec)/sqrt(100*3*6),
                                                                                   mean_rank_Rec = mean(rank_Rec), se_rank_Rec = sd(rank_Rec)/sqrt(100*3*6)) 
rank_tab_Prec_Rec
# kable(rank_tab_Prec_Rec,digits=3,format="latex",booktabs=TRUE)

# rank_tab <- matrix(NA,9,5)
# colnames(rank_tab) <- c("Method","rMSPE","pAUC","Precision","Recall")

rank_tab <- matrix(NA,length(show_methods),3)
colnames(rank_tab) <- c("Method","rMSPE","pAUC")

rank_tab[,1] <- as.character(rank_tab_rMSPE_pAUC[,1]$Method)
rank_tab[,2] <- apply(round(cbind(rank_tab_rMSPE_pAUC[,2],rank_tab_rMSPE_pAUC[,3]),3),1,function(row)paste0(row[1]," (",row[2],")"))
rank_tab[,3] <- apply(round(cbind(rank_tab_rMSPE_pAUC[,4],rank_tab_rMSPE_pAUC[,5]),3),1,function(row)paste0(row[1]," (",row[2],")"))
# rank_tab[c(1,2,3,7,8),4] <- apply(round(cbind(rank_tab_Prec_Rec[,2],rank_tab_Prec_Rec[,3]),3),1,function(row)paste0(row[1]," (",row[2],")"))
# rank_tab[c(1,2,3,7,8),5] <- apply(round(cbind(rank_tab_Prec_Rec[,4],rank_tab_Prec_Rec[,5]),3),1,function(row)paste0(row[1]," (",row[2],")"))
# copy console output to tex
kable(rank_tab,format="latex",booktabs=TRUE)

# plot group medium over p
tmp_df <- mydf_all %>% filter(Method %in% show_methods,Method!="RP_CW",cov_setting=="group",act_setting=="medium",n==200,snr==10) 
tmp_df$p <- factor(tmp_df$p,levels=c("10000","2000","500"))
tmp_df %>%
  pivot_longer(c(rMSPE,pAUC),names_to = "Measure",values_to = "Value") %>% # mutate("Meas_Cat"=ifelse(Measure=="rMSPE","rMSPE","F1/Prec/Rec")) %>%
  ggplot(aes(x=Method,y=Value,fill=Method,linetype=isSparse)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  facet_grid(Measure~p, scales = "free_y")  +
  labs(y=" ") +
  scale_linetype(guide="none")+
  theme(legend.position = "none") +
  geom_hline(data=data.frame(yintercept = 1,Measure="rMSPE"),aes(yintercept=yintercept),linetype=2)
# ggsave(paste0("../plots/SPAR_CV_group_medium_decP.pdf"), height = 5, width = 10)

# plot group medium over n 

mydf_all %>% filter(Method %in% show_methods,Method!="RP_CW",cov_setting=="group",act_setting=="medium",p==2000,snr==10,n<500) %>%
  pivot_longer(c(rMSPE,pAUC),names_to = "Measure",values_to = "Value") %>% #mutate("Meas_Cat"=ifelse(Measure=="rMSPE","rMSPE","F1/Prec/Rec")) %>%
  ggplot(aes(x=Method,y=Value,fill=Method,linetype=isSparse)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  facet_grid(Measure~n, scales = "free_y")  +
  # coord_cartesian(ylim=c(0,1.2)) +
  labs(y=" ")+
  scale_linetype(guide="none")+
  theme(legend.position = "none") +
  geom_hline(data=data.frame(yintercept = 1,Measure="rMSPE"),aes(yintercept=yintercept),linetype=2)
# ggsave(paste0("../plots/SPAR_CV_group_medium_incN.pdf"), height = 5, width = 10)

# plot group medium over snr 

mydf_all %>% filter(Method %in% show_methods,Method!="RP_CW",cov_setting=="group",act_setting=="medium",p==2000,snr<20,n==200) %>%
  pivot_longer(c(rMSPE,pAUC),names_to = "Measure",values_to = "Value") %>% #mutate("Meas_Cat"=ifelse(Measure=="rMSPE","rMSPE","F1/Prec/Rec")) %>%
  ggplot(aes(x=Method,y=Value,fill=Method,linetype=isSparse)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  facet_grid(Measure~snr, scales = "free_y")  +
  # coord_cartesian(ylim=c(0,1.2)) +
  labs(y=" ")+
  scale_linetype(guide="none") +
  theme(legend.position = "none") +
  geom_hline(data=data.frame(yintercept = 1,Measure="rMSPE"),aes(yintercept=yintercept),linetype=2)
# ggsave(paste0("../plots/SPAR_CV_group_medium_incSNR.pdf"), height = 5, width = 10)

# plot Ctime
myp <- c(500,2000,10^4)
mydf_time_sum <- mydf_all %>% filter(n==200,snr==10,cov_setting=="group",Method%in%show_methods) %>% group_by(Method,p) %>%
  summarise(Time=mean(Time)) %>% mutate(is_ref=FALSE)
mydf_time_sum <- rbind(mydf_time_sum,
                       data.frame(p=rep(myp,3),
                                  Time=c(log(myp)*myp/400,myp/300,log(myp)/20),
                                  Method=c(rep("O(plog(p))",length(myp)),rep("O(p)",length(myp)),rep("O(log(p))",length(myp))),
                                  is_ref=TRUE)
                       )
mydf_time_sum$Method <- factor(mydf_time_sum$Method,levels=c(methods[c(2,4:8,17,9:16,3,1)],"O(plog(p))","O(p)","O(log(p))"))
mydf_time_sum %>% 
  ggplot(aes(x=p,y=Time,col=Method,linetype=is_ref)) +
  geom_line() +
 scale_color_manual(values=c(scales::hue_pal()(length(show_methods)+1)[1:length(show_methods)],"darkgrey","darkgrey","darkgrey"))+ # careful about length!
  scale_x_log10() +
  scale_y_log10() +
  geom_text_repel(data = filter(mydf_time_sum,p==10000),
                  aes(x=p,y=Time,label=Method),show.legend = FALSE) +
  theme(legend.position="none") +
  labs(y="Time in s")

# ggsave("../plots/SPAR_CV_CompTime_p.pdf", height = 6, width = 9)


## # # rank sign test


cols <- c("rMSPE", "pAUC")
signs <- c(1,-1)

nem_ranks <- data.frame()
for(i in 1:length(cols)){
  rk_methods <- show_methods
  tmp <- mydf_all  %>% filter(Method %in% rk_methods
                              # cov_setting != "extreme",
                              # cov_setting != "factor"
                              # p==2000,n==200,snr==10
                              ) %>%
    select(Method,rep,setting, "score" = cols[i]) %>%
    mutate(score = score*signs[i]) %>%
    pivot_wider(names_from = Method, values_from = score) %>%
    select(-c(rep,setting)) %>%
    tsutils::nemenyi(plottype = "none", conf.level = 0.99)
  tmp2 <- t(rbind("avg_rank" = tmp$means, "lower"= tmp$intervals[1,], "upper"= tmp$intervals[2,]))
  nem_ranks <- rbind(nem_ranks,data.frame("score" = cols[i],tmp2, Method=row.names(tmp2)))
}

# str(nem_ranks)
nem_ranks$Method <- factor(nem_ranks$Method,levels = methods[c(2,4:8,17,9:16,3,1)])
nem_ranks$score <- factor(nem_ranks$score,levels = c("rMSPE","pAUC"))
# require(forcats)
# nem_ranks %>%
#   data.frame() %>%
#   mutate(Method = factor(Method, labels = rank_methods),
#          Method = fct_relevel(Method, "Ridge","LASSO","AdLASSO","ElNet","SIS","SVM","RF","RP_CW_Ensemble","TARP","SPAR","SPAR CV"),
#          score = factor(score, labels = c("pred_error","rMSLE","pAUC")),
#          score = fct_relevel(score, "pred_error","rMSLE","pAUC")
#          )


# Step 2: Identify the method with the highest median for each facet
nem_ranks_high <- nem_ranks %>%
  group_by(score) %>%
  mutate(is_best = ifelse(lower <= min(upper, na.rm = TRUE), TRUE, FALSE)) %>%
  ungroup()

# Step 3: Join this information back to the original data
nem_ranks_high <- left_join(nem_ranks, nem_ranks_high)

p_test <- ggplot(nem_ranks_high, aes(x = Method, y = avg_rank, color = is_best)) +
  geom_point(position = position_dodge(0.5), size = 0.5) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                width = 0.3, position = position_dodge(0.5)) +
  facet_grid(score ~ .,scales = "free_y") +
  # theme_bw() +
  theme(axis.text.x=element_text(angle=30, hjust=1.1, vjust = 1.05)) +
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "darkgray")) +
  guides(color = "none") +
  labs(x = "Method", y = "Mean ranks")
p_test
# ggsave("../plots/benchmark_ranks.pdf", height = 5, width = 8)
# ggsave("../plots/benchmark_ranks_HOLPScr.pdf", height = 5, width = 8)

