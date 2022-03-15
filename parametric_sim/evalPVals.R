library(plyr)
library(AUC)
library(dplyr)

evalPVals <- function(resi, alpha = 0.05, pvalsType = "adjP", rawPvalsType = "rawP") {
  # Rarely a fit has failed, then we return 0 for sens, 1 for specificity and
  # NA for FDR, AUC and the lib/cons areas
  if (!is.matrix(resi)) {
    cat("Resi is not a matrix! \n")
    return(c(Sensitivity = 0, Specificity = 1, FDR = 0, AUC = 0.5))
  }
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
           DECs = length(wh.TP),
           Sensitivity = Sensitivity, 
           Specificity = Specificity, 
           FDR = FDR, 
           AUC = ifelse(class(rocObj)[1] == "try-error", NA, AUC::auc(rocObj))))
}  # END - function: evalPVals


methods_Sel = c("DESeq2_poscounts",
                "limma_voom_TMM",
                "edgeR_TMM_robustDisp", "GLMM")     
out.dir = "/blackhole/alessia/GLMM_article/parametric_sim/"
library(purrr)
library(cowplot)
stats_sim <- function(matrix){
  # matrix : count.matrices for one verification and evaluation
  n.circ <- nrow(matrix)
  zero <- sum(matrix==0)/length(matrix)
  summarydata <- data.frame(zero,
                            n.circ)
  names(summarydata) <- c("perc.zeros", "n.circular")
  return(summarydata)
}

evals_file="/blackhole/alessia/GLMM_article/parametric_sim/ALZ_4detmet_evals_parametricsimulations_power.RDS" #from eval_function_call.R
evals_file="/blackhole/alessia/CircModel/parametric_sim/IPF_detmet_evals_allGLMM_parametricsimulations_S123_power.RDS" #from eval_function_call.R

sim_flow_file="/blackhole/alessia/GLMM_article/parametric_sim/ALZ_simulation_flow.RData"
sim_flow_file="/blackhole/alessia/CircModel/data/IPF_simulation_flow.RData"
evals <- readRDS(file = evals_file)
evals_noTM <- lapply(evals, function(methods){methods %>% purrr::list_modify("truemodel" = NULL)})
evals_noTM <- lapply(evals_noTM, function(methods){methods %>% purrr::list_modify("Y" = NULL)})
evals_noTM = purrr::map(evals_noTM, function(x) purrr::keep(x, .p = names(x)%in%methods_Sel))

object = evals_noTM$`simulation: 1_sampleSize: 3_TPR:0.1_foldEffect:1.5_seed:210192957`
m = "limma_voom_TMM"
groundTruths = function(ob, m){
  DEC = grep("TP", rownames(as.matrix(ob[[m]]$pValMat)))
  dt = data.frame(status = ifelse((1:nrow(as.matrix(ob[[m]]$pValMat)) %in% DEC)==TRUE, 1,0))
  dt
}

object = sb
method = "limma_voom_TMM"
assay = "status"
metricsDf <- do.call(rbind, lapply(names(object), function(method) {
  #or <- order(assays(object)[[assay]][, method])
  or <- order(as.matrix(object[[method]]$pValMat[,"adjP"]))
  
  TPR <- cumsum(groundTruths(ob = object, m = method)[or,1])/sum(groundTruths(ob = object, m = method)[,1])
  FDR <- cumsum(abs(groundTruths(ob = object, m = method)[or,1] - 
                      1))/seq_along(or)
  data.frame(method = method, TPR = TPR, FDR = FDR)
}))
library(ggplot2)
ggplot(metricsDf, aes(FDR, TPR, col = method)) + geom_line(size = 1.2, 
                                                           alpha = 0.6) + xlim(0, 1) + ylim(0, 1) + geom_abline(linetype = "dashed")


df_creator <- function(evals_file, sim_flow_file, out_dir){
  cat("Reading evals","\n")
  evals <- readRDS(file = evals_file)
  evals_noTM <- lapply(evals, function(methods){methods %>% purrr::list_modify("truemodel" = NULL)})
  evals_noTM <- lapply(evals_noTM, function(methods){methods %>% purrr::list_modify("Y" = NULL)})
  evals_noTM = purrr::map(evals_noTM, function(x) purrr::keep(x, .p = names(x)%in%methods_Sel))
  load(file = sim_flow_file)
  cat("Creating data.frame from evals","\n")
  
  # evals.glmm_df = foreach(i=1:3, .combine = rbind) %dopar% {
  #   # i=2
  #   eval.glmm <- evals.glmm[[i]]
  #   eval.glmm_noTM <- lapply(eval.glmm, function(methods){methods %>% purrr::list_modify("truemodel" = NULL)})
  #   eval.glmm_noTM <- lapply(eval.glmm_noTM, function(methods){methods %>% purrr::list_modify("Y" = NULL)})
  #   eval_stats.glmm <- ldply(.data = eval.glmm_noTM,.fun = function(m){
  #     evalPVals(resi = as.matrix(m$pValMat), alpha = 0.05, pvalsType = "rawP", rawPvalsType = "rawP")
  #     })
  #   eval_stats.glmm$method = "GLMM"
  #   if(i==1){evals_stats.glmm_df = data.frame(eval_stats.glmm[,-1], simulation_flow[i:length(eval.glmm_noTM),])}
  #   if(i==2){evals_stats.glmm_df = data.frame(eval_stats.glmm[,-1], simulation_flow[31:60,])}
  #   if(i==3){evals_stats.glmm_df = data.frame(eval_stats.glmm[,-1], simulation_flow[61:90,])}
  #                                
  #   evals_stats.glmm_df
  # }
  
  eval_stats <- ldply(.data = evals_noTM,.fun = function(methods){
    # methods = evals_noTM$`simulation:30_dataset:GLMM_sampleSize:5_FPR:0.1_foldEffect:1.5`
    ldply(.data = methods,.fun = function(m){
      # m = methods
      evalPVals(resi = as.matrix(m$pValMat), alpha = 0.05, pvalsType = "adjP", rawPvalsType = "adjP")
    })
  })
  
  colnames(eval_stats) <- c("method",colnames(eval_stats)[-1])
  nmethods <- length(unique(eval_stats$method))
  simulation_flow_df <- apply(simulation_flow[1:length(evals_noTM),], 2, 
                              function(col) sapply(col,function(cell) rep(cell,each = nmethods)))
  evals_stats_df <- data.frame(eval_stats,simulation_flow_df) 
  
  evals_stats_df$method <- factor(evals_stats_df$method)
  evals_stats_df$method <- factor(evals_stats_df$method, levels = levels(evals_stats_df$method), 
                                labels = c("DESeq2",
                                           "edgeR-robust",
                                           "voom"
                                           ))
  cat("Computing ROC from pVals","\n")
  
  eval_ROC <- ldply(.data = evals_noTM,.fun = function(methods){
    ldply(.data = methods,.fun = function(m){
      ROC <- AUC::roc(predictions = 1-as.matrix(m$pValMat)[,"adjP"], labels = as.factor(grepl(pattern = "TP", 
                                                                                              x = rownames(as.matrix(m$pValMat)))))
      AUC = pROC::roc(as.factor(grepl(pattern = "TP", x = rownames(as.matrix(m$pValMat)))), 1-as.matrix(m$pValMat)[,"adjP"])$auc
      cbind(fpr = ROC$fpr, tpr = ROC$tpr, auc = AUC)
    })
  })
  
  colnames(eval_ROC) <- c("method",colnames(eval_ROC)[-1])
  eval_ROC$method <- factor(eval_ROC$method)
  eval_ROC$method <- factor(eval_ROC$method, levels = levels(eval_ROC$method), 
                                labels = c("DESeq2",
                                           "edgeR-robust",
                                           "voom"))
  lengths_ROC <- ldply(.data = evals_noTM,.fun = function(methods){
    sum(ldply(.data = methods,.fun = function(m){
      ROC <- AUC::roc(predictions = 1-as.matrix(m$pValMat)[,"adjP"], labels = as.factor(grepl(pattern = "TP",
                                                                                              x = rownames(as.matrix(m$pValMat)))))
      return(length(ROC$tpr))
    })$V1)
  })
  simulation_flow_ROC_df <- apply(simulation_flow[1:length(evals_noTM),], 2, 
                                  function(col) unlist(mapply(col,lengths_ROC$V1,FUN = function(cell,times) rep(x = cell,times)),
                                                       use.names = FALSE))
  evals_ROC_df <- cbind(eval_ROC, simulation_flow_ROC_df)
  evals_ROC_df$dataset = factor(evals_ROC_df$dataset)
  
  cat("Summarizing ROC values","\n")
  evals_ROC_summary_df <- ddply(.data = evals_ROC_df[,-ncol(evals_ROC_df)],.variables = ~ 
                                  method + 
                                  dataset + 
                                  #distribution + 
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
                                       #distribution + 
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
  write.csv(evals_stats_df,file = paste0(out.dir,"evals_stats_glm_ALZ_df.csv"))
  write.csv(evals_ROC_df,file = paste0(out.dir,"evals_ROC_glm_ALZ_df.csv"))
  write.csv(evals_ROC_summary_df,file = paste0(out.dir,"evals_ROCF_glm_ALZ_summary_df.csv"))
  write.csv(evals_ROC_summary_mean_df,file = paste0(out.dir,"evals_ROC_glm_ALZ_summary_mean_df.csv"))
}

### Example code to generate power data.frames 
### The simulation files are heavy, for this reason they are not saved in github
### However the final data.frame is available.
df_creator(evals_file="/blackhole/alessia/CircModel/parametric_sim/ALZ_detmet_evals_all_parametricsimulations_power.RDS", #from eval_function_call.R
           sim_flow_file="/blackhole/alessia/CircModel/parametric_sim/ALZ_glmm_simulation_flow.RData", #form simulator
           out_dir="/blackhole/alessia/CircModel/parametric_sim/")

evals_ROC_df$method <- factor(evals_ROC_df$method)
levels(evals_ROC_df$method)
evals_ROC_df$method <- factor(evals_ROC_df$method, 
                              levels = levels(evals_ROC_df$method),
                              labels = c(
                                         "DESeq2",
                                         
                                         "edgeR-robust",
                                        
                                         "voom"),
                              ordered = T)
evals_ROC_df$dataset = factor(evals_ROC_df$dataset)
detection.labels <- c("CIRI", "DCC",  "findcirc", "circExplorer")

names(detection.labels) <- levels(evals_ROC_df$dataset)
levels(evals_ROC_df$method)
evals_AUC_ROC_mean <- ddply(evals_ROC_df, ~ dataset + method, function(x) 1-mean(x[,"auc"]))
rank_all <- ddply(evals_AUC_ROC_mean, ~ dataset , function(x) rank(x[,"V1"]))
mean_all <- ddply(evals_AUC_ROC_mean, ~ dataset , function(x) x[,"V1"])
colnames(rank_all)[2:ncol(rank_all)] <- levels(evals_ROC_df$method)
colnames(mean_all)[2:ncol(mean_all)] <- levels(evals_ROC_df$method)

rank_agg_sampleSize <- ddply(rank_all, ~ dataset, function(x) colMeans(x[,2:ncol(rank_all)]) )
rank_agg_tot_df_melted <- melt(rank_agg_sampleSize, id.vars = "dataset")
mean_agg_tot_df_melted <- melt(mean_all, id.vars = "dataset")

rank_agg_tot_df_melted$mean <- 1-mean_agg_tot_df_melted$value
simulations_summary <- ddply(rank_agg_tot_df_melted,.variables = ~ variable,function(x) mean(x[,"value"]))
colnames(simulations_summary) <- c("method","value")
saveRDS(simulations_summary, file = "/blackhole/alessia/CircModel/data/simulations_summary.RDS")

library(RColorBrewer)
cols <- c(
  # DEseq
  brewer.pal(n = 9, "YlOrRd")[c(5)],
  # Edger
  brewer.pal(n = 9, "GnBu")[c(7)],
  # limma
  brewer.pal(n = 9, "RdPu")[c(5)]
)

methods <- c("DESeq2",
             "edgeR",
             "voom")

names(cols) <- methods

methods2 <- c("DESeq2",
              "edgeR",
              "voom")

names(cols) <- methods2

png(file = "/blackhole/alessia/CircModel/parametric_sim/Figure/AUC_ALZSim.png",
     width = 15, height = 8, units = "in", res = 150)
ggplot(evals_ROC_df , aes(x=reorder(method, -auc), y = auc, color = method)) + 
  geom_boxplot() + 
  ylab("AUC") +
  xlab("") +
  labs(fill = "DE methods") +
  scale_color_manual(values = cols) +
  facet_wrap( ~ dataset) + #,
             # labeller = labeller(dataset = detection.labels)) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(face = "bold", color = "black", 
                                   size = 15, angle = 45, vjust=.97, hjust=1.1),
        axis.title.y = element_text(face = "bold", color = "black", 
                                    size = 17),
        axis.text.y = element_text(color = "black", 
                                   size = 17),
        strip.text.x = element_text(face = "bold.italic", size = 12),
        strip.text.y = element_text(face = "bold.italic", size = 12), #angle = 75),
        strip.background = element_rect(#fill = "lightblue",
          colour = "grey", size = 1),
        title = element_text(size = 14)) +
  ggtitle("AUC across sample size - Dataset ALZ\n(n.sim,30; ZINB distr.; FC,1.5; TPR,0.1)")
dev.off()

# evals_stats_df <- readRDS("/blackhole/alessia/CircModel/parametric_sim/evals_stats_ALZ_df.RDS")

evals_stats_df$method <- factor(evals_stats_df$method)
evals_stats_df$method <- factor(evals_stats_df$method, 
                               levels = c("circMeta",
                                          # "DESeq2",
                                          "DESeq2-GamPoi",
                                          "DESeq2-ZINB Wave",
                                          "DESeq2-ZeroInflated",
                                          "edgeR-robust",
                                          "edgeR",
                                          "edgeR-ZINB Wave",
                                          "GLMM",
                                          "voom"),
                                          # "voom-ZINB Wave"),
                               ordered = T)
p <- evals_stats_df %>% 
  dplyr::group_by(method, sampleSize) %>% 
  dplyr::summarise(sens.mean = mean(Sensitivity),
            spec.mean  = mean(Specificity)) %>% 
  ggplot(aes(y=sens.mean, x=1-spec.mean, color=method))

png(file = "/blackhole/alessia/CircModel/parametric_sim/Figure/sens_spec_ALZ.png",
width = 15, height = 8, units = "in", res = 150)
p + geom_point(size = 5) + 
  facet_wrap( ~ sampleSize) +
  # scale_shape_manual(values=1:5) +
  scale_color_manual(values = cols) +
  ylab("Sensitivity") +
  xlab("1 - specificity (false positive rate)") + 
  # coord_cartesian(xlim=c(-.003,.035)) + 
  geom_vline(xintercept=.05) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(color = "black", 
                                   size = 15, angle = 45, vjust=.97, hjust=1.1),
        axis.title.y = element_text(face = "bold", color = "black", 
                                    size = 17),
        axis.title.x = element_text(face = "bold", color = "black", 
                                    size = 17),
        axis.text.y = element_text(color = "black", 
                                   size = 15),
        strip.text.x = element_text(face = "bold.italic", size = 12),
        strip.text.y = element_text(face = "bold.italic", size = 12), #angle = 75),
        strip.background = element_rect(#fill = "lightblue",
          colour = "grey", size = 1),
        title = element_text(size = 14)) +
  guides(shape=guide_legend(nrow=2,byrow=TRUE)) + 
  ggtitle("Power across sample size - Dataset ALZ\n(n.sim,30; ZINB distr.; FC,1.5; TPR,0.1)")
  # scale_x_continuous(breaks=c(0,.1))
dev.off()

legend = get_legend(p + geom_point(size = 3) + 
                      theme_bw() + 
                      facet_wrap( ~ sampleSize) +
                      # scale_shape_manual(values=1:5) +
                      scale_color_manual(values = cols) +
                      ylab("Sensitivity") +
                      xlab("1 - specificity (false positive rate)") + 
                      # coord_cartesian(xlim=c(-.003,.035)) + 
                      geom_vline(xintercept=.05) +
                      theme_classic() +
                      theme(legend.position = "bottom",
                            axis.text.x = element_text(#face = "bold", color = "#993333",
                              size = 9, angle = 45, vjust=.9, hjust=0.8),
                            strip.text.x = element_text(face = "bold.italic", size = 9),
                            strip.text.y = element_text(face = "bold.italic", size = 9), #angle = 75),
                            strip.background = element_rect(#fill = "lightblue",
                              colour = "grey", size = 1)) +
                      guides(shape=guide_legend(nrow=2,byrow=TRUE)) + 
                      ggtitle("Power across sample size - Dataset ALZ\n(n.sim,30; ZINB distr.; FC,1.5; TPR,0.1)"))

prow <- plot_grid( p1 ,
                   p2 ,
                   align = 'vh',
                   labels = c("a", "b"),
                   hjust = -1,
                   nrow = 1
)


p <- plot_grid( prow, legend, ncol = 1, rel_heights = c(1, .2))
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
grid.arrange(arrangeGrob(p1,
                         p2,
                         nrow=1),
             legend, nrow=2, heights=c(10, 1))
