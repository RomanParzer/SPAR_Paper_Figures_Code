## read in file and plot  -----------------------------------------------------------------------
## ------------------------------------------------------------------------------------------------------------

pacman::p_load(dplyr, ggplot2, tidyr, ggrepel,knitr,kableExtra)

resobj <- readRDS("../saved_results/synth_rateye_SPAR_25_nset12_reps100_nmeth21.rds")

res <- resobj$res

methods <- dimnames(res)[[3]]
settings <- attributes(res)$settings

show_methods <- methods[c(1,3,6:8,11:15,18:19)]
sparse_methods <- methods[c(6,7,8,13,14,15,19)]


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
mydf_all$Method <- factor(mydf_all$Method,levels =  methods[c(5:8,17,19,9:16,3,1,2,4,18,20:21)])
mydf_all <- mutate(mydf_all,"isSparse"=ifelse(Method%in%c("LASSO","AdLASSO","ElNet","SIS","HOLPScr","SLOPE"),TRUE,FALSE))
# rescale to (0,1)
mydf_all <- mutate(mydf_all,"pAUC"=pAUC*2*(p-a)/120)

# plot medium p2000 the 6 cov settings rMSPE 
mydf_all %>% filter(Method %in% show_methods,snr==10) %>%
  ggplot(aes(x=Method,y=rMSPE,fill=Method,linetype=isSparse)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  facet_grid(p~act_setting, scales = "free_y") +
  coord_cartesian(ylim=c(0,1.35))+
  theme(legend.position = "none") +
  geom_hline(yintercept=1,linetype=2)
# ggsave(paste0("../plots/SPAR_rMSPE_synthRateye_25.pdf"), height = 6, width = 8)

# pAUC
mydf_all %>% filter(Method %in% show_methods,snr==10) %>%
  ggplot(aes(x=Method,y=pAUC,fill=Method,linetype=isSparse)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  ggh4x::facet_grid2(p~act_setting, scales = "free_y",independent = "y") +
  theme(legend.position = "none")
# ggsave(paste0("../plots/SPAR_pAUC_synthRateye_25.pdf"), height = 6, width = 10)

# optional:
# # table of ranks for pAUC and rMSPE
# 
# # remove PLS because NA for p=22905
# n_showm <- length(show_methods)-1
# mydf_all_rank_rMSPE_pAUC <- mydf_all %>% filter(Method %in% show_methods,Method!="PLS",snr==10) %>% group_by(rep,setting) %>% mutate(rank_pAUC=n_showm + 1 - rank(pAUC) ,rank_rMSPE=rank(rMSPE))
# rank_tab_rMSPE_pAUC <- mydf_all_rank_rMSPE_pAUC %>% group_by(Method) %>% summarise(mean_rank_rMSPE = mean(rank_rMSPE), se_rank_rMSPE = sd(rank_rMSPE)/sqrt(100*3*6),
#                                                                                    mean_rank_pAUC = mean(rank_pAUC), se_rank_pAUC = sd(rank_pAUC)/sqrt(100*3*6)) 
# rank_tab_rMSPE_pAUC
# 
# rank_tab <- matrix(NA,n_showm,5)
# # colnames(rank_tab) <- c("Method","rMSPE","pAUC")
# 
# rank_tab[,1] <- as.character(rank_tab_rMSPE_pAUC[,1]$Method)
# rank_tab[,2] <- round(unlist(rank_tab_rMSPE_pAUC[,2]),3)
# rank_tab[,3] <- sapply(round(rank_tab_rMSPE_pAUC[,3],3),function(tmp)paste0("(",tmp,")"))
# rank_tab[,4] <- round(unlist(rank_tab_rMSPE_pAUC[,4]),3)
# rank_tab[,5] <- sapply(round(rank_tab_rMSPE_pAUC[,5],3),function(tmp)paste0("(",tmp,")"))
# 
# # make top 3 bold
# ind_best_rMSPE <- order(unlist(rank_tab_rMSPE_pAUC[,2]))[1:3]
# ind_best_pAUC <- order(unlist(rank_tab_rMSPE_pAUC[,4]))[1:3]
# rank_tab[ind_best_rMSPE,2] <- sapply(rank_tab[ind_best_rMSPE,2],function(tmp)paste0("\textbf{",tmp,"}"))
# rank_tab[ind_best_pAUC,4] <- sapply(rank_tab[ind_best_pAUC,4],function(tmp)paste0("\textbf{",tmp,"}"))
# 
# # copy console output to tex
# kable(rank_tab,format="latex",booktabs=TRUE) %>% 
#   add_header_above(header=c("Method"=1,"rMSPE"=2,"pAUC"=2))

