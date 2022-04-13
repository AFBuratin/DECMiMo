library(plyr)
library(AUC)
library(dplyr)

evalPVals <- function(resi, alpha = 0.05, pvalsType = "adjP", rawPvalsType = "rawP") {
  # Rarely a fit has failed, then we return 0 for sens, 1 for specificity and
  # NA for FDR, AUC and the lib/cons areas
  # if (!is.matrix(resi)) {
  #   cat("Resi is not a matrix! \n")
  #   return(c(Sensitivity = 0, Specificity = 1, FDR = 0, AUC = 0.5))
  # }
  # Some DESeq2 results (for example) had NA adjusted p-values Replace NA
  # $adjP values to highest possible value (1.0)
  # resi[is.na(resi[, pvalsType]), pvalsType] <- 1
  # Or just count them
  NA_prop <- sum(is.na(resi[, pvalsType]))/nrow(resi)
  
  # Evaluate detection performance.
  wh.pred = (resi[, pvalsType] < alpha)
  wh.pos = which(wh.pred)
  wh.neg = which(!wh.pred)
  wh.TP = grep("TP", rownames(resi))
  # Calc number of differentially abundant circRNAs
  FPs = sum(!wh.pos %in% wh.TP)
  TPs = sum(wh.pos %in% wh.TP)
  TNs = sum(!wh.neg %in% wh.TP)
  FNs = sum(wh.neg %in% wh.TP)
  # Sensitivity: True Positives divided by all positives (sum of true
  # positives and false negatives)
  Sensitivity = TPs/(TPs + FNs)
  # Specificity: True negatives divided by all negatives (sum of true
  # negatives and false positives)
  Specificity = TNs/(TNs + FPs)
  # false discovery rate: false positives divided by all detected positives
  FDR = if ((TPs + FPs) == 0) 0 else FPs/(TPs + FPs)
  
  # If no true positives, return NA's for irrelevant measures
  wh.truth = (1:nrow(resi) %in% wh.TP)
  
  # IF AUC cannot be calculated, return NA
  rocObj = try(AUC::roc(1 - resi[, pvalsType], factor(as.numeric(wh.truth))))
  return(c("NA_Proportion" = NA_prop, 
           TPs = length(wh.TP),
           DECs = length(wh.pos),
           Sensitivity = Sensitivity, 
           Specificity = Specificity, 
           FDR = FDR, 
           AUC = ifelse(class(rocObj)[1] == "try-error", NA, AUC::auc(rocObj))))
}  # END - function: evalPVals
stats_sim <- function(matrix){
  
  # matrix : count.matrices for one verification and evaluation
  
  n.circ <- nrow(matrix)
  zero <- sum(matrix==0)/length(matrix)
  summarydata <- data.frame(zero,
                            n.circ)
  
  names(summarydata) <- c("perc.zeros", "n.circular")
  
  return(summarydata)
}
dataset <- "GLMM"
simulation = 1:30
sampleSize = c(5)
FPR = c(0.1)
foldEffect <- c(1.5)
#sparsityEffect <- c(0.15, 0.30)
niter = 30

simulation_flow.glmm <- data.frame(expand.grid(simulation = simulation,
                                               dataset = dataset,
                                               sampleSize = sampleSize,
                                               FPR = FPR,
                                               foldEffect = foldEffect 
                                               #sparsityEffect = sparsityEffect
))
rownames(simulation_flow.glmm) <- 1:nrow(simulation_flow.glmm)
dim(simulation_flow.glmm)
df_creator <- function(evals_file, sim_flow_file, out_dir){
  cat("Reading evals","\n")
  evals = readRDS("/blackhole/alessia/GLMM_article/parametric_sim/ALZ_glmm_parametricsimulations_power.RDS")

  names(evals) = apply(simulation_flow.glmm, 1, function(sim) paste(colnames(simulation_flow.glmm),
                                                                                 sim,
                                                                                 sep = ":",
                                                                                 collapse = "_"))
  
  cat("Creating data.frame from evals","\n")
  
  eval_stats <- ldply(.data = evals, .fun = function(simulation){
    # simulation=evals$`simulation: 1_dataset:GLMM_sampleSize:5_FPR:0.1_foldEffect:1.5`
      evalPVals(resi = simulation$pValMat, alpha = 0.05, pvalsType = "adjP", rawPvalsType = "rawP")

  })
  colnames(eval_stats) <- c("sim",colnames(eval_stats)[-1])
  simulation_flow.glmm$method = "GLMM"
  simulation_flow_df <- apply(simulation_flow.glmm, 2, 
                              function(col) sapply(col,function(cell) rep(cell,each = 1)))
  evals_stats_df <- data.frame(eval_stats,simulation_flow_df) 
  sim_stats_df <- data.frame(eval_stats,simulation_flow_df) 
  
  evals_stats_df$method <- factor(evals_stats_df$method)
  evals_stats_df$method <- factor(evals_stats_df$method, levels = levels(evals_stats_df$method), 
                                labels = c("GLMM"))
  cat("Computing ROC from pVals","\n")
  eval_ROC <- 
    ldply(.data = evals,.fun = function(m){
      ROC <- AUC::roc(predictions = 1-m$pValMat[,"adjP"], labels = as.factor(grepl(pattern = "TP", x = rownames(m$pValMat))))
      AUC = pROC::roc(as.factor(grepl(pattern = "TP", x = rownames(m$pValMat))), 1-m$pValMat[,"adjP"])$auc
      cbind(fpr = ROC$fpr, tpr = ROC$tpr, auc = AUC)
    })

  colnames(eval_ROC) <- c("sim",colnames(eval_ROC)[-1])
  
  lengths_ROC <- 
    sum(ldply(.data = evals,.fun = function(m){
      ROC <- AUC::roc(predictions = 1-m$pValMat[,"adjP"], labels = as.factor(grepl(pattern = "TP",x = rownames(m$pValMat))))
      return(length(ROC$tpr))
    })$V1)

  simulation_flow_ROC_df <- apply(simulation_flow.glmm, 2, 
                                  function(col) unlist(mapply(col,lengths_ROC,FUN = function(cell,times) rep(x = cell,times)),
                                                       use.names = FALSE))
  evals_ROC_df <- cbind(eval_ROC, simulation_flow_ROC_df)
  cat("Summarizing ROC values","\n")
  evals_ROC_summary_df <- ddply(.data = evals_ROC_df,.variables = ~ 
                                  method + 
                                  dataset + 
                                  # distribution + 
                                  sampleSize +
                                  simulation +
                                  FPR +
                                  foldEffect, 
                                  # compensation + 
                                  # sparsityEffect,
                                .fun = function(x){
                                    support <- seq(0,1,length.out = 101)
                                    fpr_tpr <- data.frame(fpr = 0, tpr = 0)
                                    for(i in 2:length(support)){
                                      fpr_s <- support[i-1]
                                      if(sum(x$fpr>=support[i-1] & x$fpr<support[i]) > 0)
                                        tpr_s <- mean(x$tpr[x$fpr>=support[i-1] & x$fpr<support[i]])
                                      else tpr_s <- fpr_tpr[i-2,2]
                                      fpr_tpr[i-1,] <- c(fpr_s,tpr_s)
                                    }
                                    fpr_tpr[1,] <- c(0,0)
                                    fpr_tpr[length(support),] <- c(1,1)
                                    return(fpr_tpr)
                                  })
  evals_ROC_summary_mean_df <- ddply(.data = evals_ROC_summary_df,.variables = ~ 
                                       method + 
                                       dataset + 
                                       # distribution + 
                                       sampleSize +
                                       simulation +
                                       FPR +
                                       foldEffect + 
                                       fpr,
                                     # compensation + 
                                     # sparsityEffect,
                                     .fun = function(x){
                                         tpr = mean(x$tpr)
                                         se = sqrt(var(x$tpr))
                                         return(data.frame(tpr = tpr, se = se))
                                       })
  cat("Saving data","\n")
  out.dir = "/blackhole/alessia/GLMM_article/parametric_sim/"
  
  write.csv(evals_stats_df,file = paste0(out.dir,"evals_stats_glmm_ALZ_df.csv"))
  write.csv(evals_ROC_df,file = paste0(out.dir,"evals_ROC_glmm_ALZ_df.csv"))
  write.csv(evals_ROC_summary_df,file = paste0(out.dir,"evals_ROC_glmm_ALZ_summary_df.csv"))
  write.csv(evals_ROC_summary_mean_df,file = paste0(out.dir,"evals_ROC_glmm_ALZ_summary_mean_df.csv"))
}

### Example code to generate power data.frames 
### The simulation files are heavy, for this reason they are not saved in github
### However the final data.frame is available.
df_creator(evals_file="/blackhole/alessia/CircModel/parametric_sim/ALZ_detmet_evals_all_parametricsimulations_power.RDS", #from eval_function_call.R
           sim_flow_file="/blackhole/alessia/CircModel/parametric_sim/ALZ_glmm_simulation_flow.RData", #form simulator
           out_dir="/blackhole/alessia/CircModel/parametric_sim/")

evals_ROC_glmm_df <- read.csv("/blackhole/alessia/GLMM_article/parametric_sim/evals_ROC_glmm_IPF_df.csv")
evals_ROC_met_df <- read.csv("/blackhole/alessia/GLMM_article/parametric_sim/evals_ROC_glm_IPF_df.csv")

evals_ROC_glmm_df <- read.csv("/blackhole/alessia/GLMM_article/parametric_sim/evals_ROC_glmm_ALZ_df.csv")
evals_ROC_met_df <- read.csv("/blackhole/alessia/GLMM_article/parametric_sim/evals_ROC_glm_ALZ_df.csv")
evals_ROC_df = rbind.fill(evals_ROC_glmm_df, evals_ROC_met_df)
evals_ROC_df$method <- factor(evals_ROC_df$method)
evals_ROC_df$method <- factor(evals_ROC_df$method, 
                              levels = c("DESeq2",
                                         "edgeR-robust",
                                         "GLMM",
                                         "voom"),
                              ordered = T)

evals_ROC_df %>% group_by(method) %>% dplyr::summarise(mean = mean(fpr))
evals_ROC_df %>% group_by(method) %>% dplyr::summarise(mean = mean(auc))
evals_ROC_df %>% group_by(method) %>% dplyr::summarise(mean = mean(tpr))

evals_AUC_ROC_mean <- ddply(evals_ROC_df, ~ method, function(x) 1-mean(x[,"auc"]))
rank_all <- data.frame(rank(evals_AUC_ROC_mean[,"V1"]))
rank_all$method = evals_AUC_ROC_mean$method
colnames(rank_all) <- c("value","method")
saveRDS(rank_all, file = "/blackhole/alessia/GLMM_article/data/simulations_summary.RDS")

library(RColorBrewer)
cols <- c(
  # GLMM
  brewer.pal(n = 9, "BuPu")[c(5)],
  
  # DEseq
  brewer.pal(n = 9, "YlOrRd")[c(3)],
  # Edger
  brewer.pal(n = 9, name = "GnBu")[c(5)],
  # limma
  brewer.pal(n = 9, "RdPu")[c(5)]
)

methods2 <- c("GLMM",
              
              "DESeq2",
             
              "edgeR-robust",
              
              "voom")
names(cols) <- methods2

png(file = "/blackhole/alessia/GLMM_article/parametric_sim/AUC_IPFSim_glmm.png",
     width = 8, height = 8, units = "in", res = 300)
ggplot(evals_ROC_df, aes(x=method, y = auc, color = method)) + 
  geom_boxplot() + 
  ylab("AUC") +
  xlab("") +
  labs(fill = "DE methods") +
  scale_color_manual(values = cols) +
  theme_bw() + 
  theme(axis.text.x = element_text(face = "bold", color = "black", 
                                   size = 14, angle = 45, vjust=.97, hjust=1.1),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(face = "bold.italic", size = 14),
        strip.text.y = element_text(face = "bold.italic", size = 14), #angle = 75),
        strip.background = element_rect(#fill = "lightblue",
          colour = "grey", size = 1),
        title = element_text(size = 14)) +
  ggtitle("AUC across sample size - Dataset IPF\n(n.sim,30; NB distr.; FC,1.5; TPR,0.1)")
dev.off()

evals_stats_met_df <- read.csv("/blackhole/alessia/GLMM_article/parametric_sim/evals_stats_glm_IPF_df.csv")
evals_stats_glmm_df <- read.csv("/blackhole/alessia/GLMM_article/parametric_sim/evals_stats_glmm_IPF_df.csv")
evals_stats_df = rbind.fill(evals_stats_glmm_df, evals_stats_met_df)
evals_stats_df$method <- factor(evals_stats_df$method)
evals_stats_df$method <- factor(evals_stats_df$method, 
                               levels = c(
                                          "DESeq2",
                                          
                                          "edgeR-robust",
                                          "GLMM",
                                          "voom"),
                               ordered = T)
evals_stats_df %>% group_by(method) %>% dplyr::summarise(mean = mean(Sensitivity))
evals_stats_df %>% group_by(method) %>% dplyr::summarise(mean = mean(Specificity))
evals_stats_df %>% group_by(method) %>% dplyr::summarise(mean = mean(FDR))

p <- evals_stats_df %>% 
  dplyr::group_by(method) %>% 
  dplyr::summarise(sens.mean = mean(Sensitivity),
            spec.mean  = mean(Specificity)) %>% 
  ggplot(aes(y=sens.mean, x=1-spec.mean, color=method))

png(file = "/blackhole/alessia/GLMM_article/parametric_sim/sens_spec_IPF_glmm.png",
     width = 8, height = 8, units = "in", res = 300)
p + geom_point(size = 3) + 
  theme_bw() + 
  scale_color_manual(values = cols) +
  ylab("Sensitivity") +
  xlab("1 - specificity (false positive rate)") + 
  # coord_cartesian(xlim=c(-.003,.035)) + 
  geom_vline(xintercept=.05) +
  theme_bw() + 
  theme(axis.text.x = element_text(face = "bold", color = "black", 
                                   size = 14, angle = 45, vjust=.97, hjust=1.1),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(face = "bold.italic", size = 14),
        strip.text.y = element_text(face = "bold.italic", size = 14), #angle = 75),
        strip.background = element_rect(#fill = "lightblue",
          colour = "grey", size = 1),
        title = element_text(size = 14)) +
  guides(shape=guide_legend(nrow=2,byrow=TRUE)) + 
  ggtitle("Power across sample size - Dataset IPF\n(n.sim,30; NB distr.; FC,1.5; TPR,0.1)")
  # scale_x_continuous(breaks=c(0,.1))
dev.off()


