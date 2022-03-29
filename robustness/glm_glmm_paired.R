# load libraries -----

library("GenomicRanges")
library("Biobase")
library("BiocParallel")
pkgs <- c("edgeR", 
          "glmGamPoi",
          "DESeq2",
          "phyloseq", 
          "plyr", 
          "reshape2",
          "ROCR",
          # "samr",
          "apeglm",
          "ZIM",
          "zinbwave",
          "AUC",
          "genefilter",
          "MAST",
          # "scde", # It is important to use flexmix v2.3-13 and scde 1.99.1
          # "Seurat",
          "crayon",
          "aod",
          "arm",
          "fdrtool",
          "lars",
          "emdist",
          # "baySeq",
          "snow",
          # "ShrinkBayes",
          # "DEsingle",
          "VGAM",
          "Matrix",
          "maxLik",
          "MASS",
          # "samr",
          #"EBseq",
          "NewWave",
          "limma",
          "data.table",
          "dplyr")
for(i in pkgs) { library(i, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE, character.only=TRUE) }

library(devtools)
install_github("egaffo/CREART", ref = "dev", force = T)
library(CREART)

# load data and meta data from previous data formatting in GOF.Rmd step -------


## DM1 data set
randomSubsets <- read.table("/blackhole/alessia/GLMM_article/robustness/random_subsets_eval_veri.txt",strings=FALSE)

meta.data <- read.csv("/blackhole/alessia/circzi/checkCircRNAnormalizationdistribution/realdata/DM1/meta_DM1.csv")
meta.data = as.data.table(meta.data)
meta.data <- meta.data[, .(sample_id = sample,
                           condition = ifelse(disease_class=="myotonic dystrophy type 1", "DM1","Normal"))]
meta.data = meta.data[order(meta.data$sample_id),][seq(1,nrow(meta.data), by = 2),]     # for odd rows

coldata <- DataFrame(group = meta.data$condition,
                     sample = meta.data$sample,
                     row.names = meta.data$sample)
coldata$group <- factor(coldata$group)
coldata$sample <- as.character(coldata$sample)
coldata

## IPF data set
randomSubsets <- read.table("/blackhole/alessia/CircModel/power/IPF_random_subsets_eval_veri.txt",strings=FALSE)

meta.data <- read.csv("/blackhole/alessia/circzi/checkCircRNAnormalizationdistribution/realdata/IPF/analyses/meta_IPF.csv")
meta.data = meta.data[seq(1,nrow(meta.data), by = 2),]     # for odd rows

coldata <- DataFrame(condition = meta.data$condition,
                     group = ifelse(meta.data$condition=="normal", "normal", "IPF"),
                     sample_id = meta.data$sample,
                     row.names = meta.data$sample)
coldata$condition <- factor(coldata$condition)
coldata$group <- factor(coldata$group)
coldata$sample_id <- as.character(coldata$sample_id)
coldata

ccp2 = read.table("/blackhole/alessia/CircModel/data/IPF_ccp2.csv", header = T)
# ccp2 = RCurl::scp(host = "threesum", path = "/home/enrico/analysis/zicirc/IPF/data/IPF_ccp2.csv", keypasswd = "alessiasucks666", 
#            user="alessia")
load("/blackhole/alessia/CircModel/data/IPFData_list.RData")

# ccp2 = IPFData_list$ccp2
chr <- sub(":.*", "", ccp2$circ_id)
start <- as.numeric(gsub(".*:(.*)\\-.*", "\\1", ccp2$circ_id))
end <- sub(".*-", "\\1", ccp2$circ_id)
circ_id = paste0(chr, ":", start+1, "-", end)
ccp2$circ_id = circ_id
ccp2.df = ccp2[,-1]
rownames(ccp2.df) = ccp2$circ_id
se.ccp2 <- SummarizedExperiment(assays = list(counts = as.matrix(ccp2.df)),
                                colData = coldata)
eset.ccp2 <- ExpressionSet(assay(se.ccp2),
                           AnnotatedDataFrame(as.data.frame(colData(se.ccp2))))
pData(eset.ccp2)$condition <- factor(pData(eset.ccp2)$group)
levels(pData(eset.ccp2)$condition) <- c("A", "B")
e = eset.ccp2


## ALZ data set
randomSubsets <- read.table("/blackhole/alessia/GLMM_article/robustness/ALZ_random_subsets_eval_veri.txt",strings=FALSE)

meta.data <- read.csv("/blackhole/alessia/GLMM_article/data/meta_alz.csv")
meta.data

coldata <- DataFrame(condition = meta.data$condition,
                     group = meta.data$group,
                     sample = meta.data$sample,
                     row.names = meta.data$sample)
coldata$condition <- factor(coldata$condition)
coldata$group <- factor(coldata$group)
coldata$sample <- as.character(coldata$sample)
coldata

load(file = "/blackhole/alessia/GLMM_article/data/ALZData_list.RData") 
load(file = "/blackhole/alessia/CircModel/data/ALZZinb_nb_Fit_detmet_models.RData") #created in datasets_and_models.R
ALZDataFilt_list = lapply(ALZData_list, function(dat){
  # dat=ALZData_list$findcirc
  new.data = as.data.frame(dat)
  new.data$circ_id = rownames(dat)
  filt.dat = CREART::smallest_group_filter(x = as.data.table(new.data), 
                                           cond = as.data.table(coldata),
                                           rthr = 1)
  count.matrix = as.matrix(filt.dat[,-"circ_id"]) 
  rownames(count.matrix) = filt.dat$circ_id
  return(count.matrix)
})
methods_det <- c("findcirc", "dcc", "ciri", "circexplorer2_star")

data_list = ALZDataFilt_list[methods_det]
models = ALZData_models[methods_det]

library(SummarizedExperiment)
se <- lapply(data_list, function(x) SummarizedExperiment(x[,colnames(x)%in%coldata$sample], colData = coldata))
eset <- lapply(se, function(x) ExpressionSet(assay(x),
                                             AnnotatedDataFrame(as.data.frame(colData(x)))))

e <- lapply(eset, function (x) {
  pData(x)$condition <- factor(pData(x)$group)
  levels(pData(x)$condition) <- c("A", "B")
  x.new = x[which(rowSums(exprs(x))!=0),]
  return(x.new)})


## TALL data set
randomSubsets <- read.table("/blackhole/alessia/CircModel/robustness/TALL_random_subsets_eval_veri.txt",strings=FALSE)

meta.data <- read.csv("/blackhole/alessia/circzi/checkCircRNAnormalizationdistribution/robustness_glmm/TALL/meta.csv")

coldata <- DataFrame(condition = meta.data$condition,
                     group = meta.data$groups,
                     sample = gsub("-", "_", meta.data$sample_id),
                     row.names = gsub("-", "_", meta.data$sample_id))
coldata$condition <- factor(coldata$group)
coldata$group <- factor(gsub("-", "_", coldata$group))
coldata$sample <- as.character(coldata$sample)
coldata

load("/blackhole/alessia/circzi/checkCircRNAnormalizationdistribution/data/TALLData_list.RData")

ccp2 = TALLData_list$ccp2

se.ccp2 <- SummarizedExperiment(assays = list(counts = as.matrix(ccp2)[, rownames(coldata)]),
                                colData = coldata)
eset.ccp2 <- ExpressionSet(assay(se.ccp2),
                           AnnotatedDataFrame(as.data.frame(colData(se.ccp2))))
pData(eset.ccp2)$condition <- factor(pData(eset.ccp2)$group)
levels(pData(eset.ccp2)$condition) <- c("A", "B")
e = eset.ccp2
# summarized experiment for GLMM matrix -----
glmm.db <- rbindlist(lapply(TALLData_list[1:6], function(x) data.frame(x, circ_id = rownames(x))), 
                     idcol = "method", use.names = TRUE)
glmm.melt <- rbindlist(lapply(TALLData_list[1:6], function(x) reshape2::melt(x)), 
                       idcol = "method", use.names = TRUE)

count.data.melt <- as.data.table(reshape2::melt(glmm.db, id.vars = c("method", "circ_id")))
name.sep <- "."
count.data.melt$variable = gsub("^X", "", count.data.melt$variable)
count.data.melt[, sample.name.ext := paste0(variable, name.sep, method)]
count.data.merge <- merge(count.data.melt, coldata, by.x = "variable", by.y = "sample")

colData.dt <- as.data.table(count.data.merge)[, .N, by = .(variable, method, 
                                                           group,
                                                           sample.name.ext)][, N := NULL][]

colData <- data.frame(colData.dt, 
                      row.names = "sample.name.ext")
table(colData$group)

glmm.wide = dcast(glmm.melt, method+Var2~Var1, value.var = "value", fill=0)
colnames(glmm.wide) = c("MethodID", "SampleID", colnames(glmm.wide)[-c(1,2)])
head(glmm.wide[,c(1:8)])
count.matrix.glmm.dt <- dcast(data = as.data.table(count.data.merge),
                              formula = circ_id ~ sample.name.ext,
                              fill = 0, fun.aggregate = sum,
                              value.var = "value")
count.matrix.filtered.glmm.dt <- as.data.table(count.matrix.glmm.dt)#[rowSums(as.data.table(count.matrix.glmm.dt)[,-"circ_id"]>=2)>1,]
count.matrix.glmm <- as.matrix(count.matrix.filtered.glmm.dt,
                               rownames = "circ_id")[, rownames(colData)]

## DM1
load("/blackhole/alessia/CircModel/data/DM1Data_list.RData")
glmm.db <- rbindlist(lapply(DM1Data_list[1:6], function(x) data.frame(x, circ_id = rownames(x))), 
                     idcol = "method", use.names = TRUE)
glmm.melt <- rbindlist(lapply(DM1Data_list[1:6], function(x) reshape2::melt(x)), 
                       idcol = "method", use.names = TRUE)

count.data.melt <- as.data.table(reshape2::melt(glmm.db, id.vars = c("method", "circ_id")))
name.sep <- "."
count.data.melt[, sample.name.ext := paste0(variable, name.sep, method)]
count.data.merge <- merge(count.data.melt, coldata, by.x = "variable", by.y = "sample")

colData.dt <- as.data.table(count.data.merge)[, .N, by = .(variable, method, 
                                                           group,
                                                           sample.name.ext)][, N := NULL][]

colData <- data.frame(colData.dt, 
                      row.names = "sample.name.ext")

glmm.wide = dcast(glmm.melt, method+Var2~Var1, value.var = "value", fill=0)
colnames(glmm.wide) = c("MethodID", "SampleID", colnames(glmm.wide)[-c(1,2)])
head(glmm.wide[,c(1:8)])

## IPF
glmm.db <- rbindlist(lapply(IPFData_list[1:6], function(x) data.frame(x, circ_id = rownames(x))), 
                     idcol = "method", use.names = TRUE)
glmm.melt <- rbindlist(lapply(IPFData_list[1:6], function(x) reshape2::melt(x)), 
                       idcol = "method", use.names = TRUE)

count.data.melt <- as.data.table(reshape2::melt(glmm.db, id.vars = c("method", "circ_id")))
name.sep <- "."
count.data.melt[, sample.name.ext := paste0(variable, name.sep, method)]
count.data.merge <- merge(count.data.melt, coldata, by.x = "variable", by.y = "sample_id")

colData.dt <- as.data.table(count.data.merge)[, .N, by = .(variable, method, 
                                                           group,
                                                           sample.name.ext)][, N := NULL][]
colData <- data.frame(colData.dt, 
                      row.names = "sample.name.ext")

# glmm.wide = dcast(glmm.melt, method+Var2~Var1, value.var = "value", fill=0)
# colnames(glmm.wide) = c("MethodID", "SampleID", colnames(glmm.wide)[-c(1,2)])
# head(glmm.wide[,c(1:8)])

## ALZ
load("/run/user/1001/gvfs/sftp:host=172.21.32.10/blackhole/alessia/GLMM_article/data/ALZData_list.RData")
ALZDataFilt_list = lapply(ALZData_list, function(dat){
  # dat=ALZData_list$findcirc
  new.data = as.data.frame(dat)
  new.data$circ_id = rownames(dat)
  filt.dat = CREART::smallest_group_filter(x = as.data.table(new.data), 
                                           cond = as.data.table(coldata),
                                           rthr = 1)
  count.matrix = as.matrix(filt.dat[,-"circ_id"]) 
  rownames(count.matrix) = filt.dat$circ_id
  return(count.matrix)
})
lapply(ALZDataFilt_list, function(x) dim(x))
glmm.db <- rbindlist(lapply(ALZDataFilt_list, function(x) data.frame(x, circ_id = rownames(x))), 
                     idcol = "method", use.names = TRUE)
methods_sel <- c("findcirc", "dcc", "ciri", "circexplorer2_star")

glmm.db = glmm.db %>% dplyr::filter(method%in%methods_sel)
glmm.melt <- rbindlist(lapply(ALZDataFilt_list, function(x) reshape2::melt(x)), 
                       idcol = "method", use.names = TRUE)
glmm.melt = glmm.melt %>% dplyr::filter(method%in%methods_sel)

count.data.melt <- as.data.table(reshape2::melt(glmm.db, id.vars = c("method", "circ_id")))
name.sep <- "."
count.data.melt[, sample.name.ext := paste0(variable, name.sep, method)]
count.data.merge <- merge(count.data.melt, coldata, by.x = "variable", by.y = "sample")
count.matrix.glmm <- dcast(count.data.melt,circ_id~sample.name.ext, value.var = "value", fill=0)
rownames(count.matrix.glmm) = count.matrix.glmm$circ_id
dim(count.matrix.glmm)
glmm.wide = dcast(count.data.melt,method+sample.name.ext~circ_id, value.var = "value", fill=0)
colnames(glmm.wide) = c("MethodID", "SampleID", colnames(glmm.wide)[-c(1,2)])
head(glmm.wide[,c(1:8)])


colData.dt <- as.data.table(count.data.merge)[, .N, by = .(variable, method, 
                                                           group,
                                                           sample.name.ext)][, N := NULL][]
colData <- data.frame(colData.dt, 
                      row.names = "sample.name.ext")
colnames(colData) = c("sample", "method", "condition")
colData$sample_id = rownames(colData)
matrix.glmm = count.matrix.glmm[,-1][,rownames(colData)]
se.glmm <- SummarizedExperiment(assays = list(counts = matrix.glmm),
                                colData = colData)
eset.glmm <- ExpressionSet(assay(se.glmm),
                           AnnotatedDataFrame(as.data.frame(colData(se.glmm))))

nreps = 30

# count.data.melt <- as.data.table(reshape2::melt(glmm.db, id.vars = c("method", "circ_id")))
# name.sep <- "."
# count.data.melt[, sample.name.ext := paste0(variable, name.sep, method)]
# count.data.merge <- merge(count.data.melt, coldata, by.x = "variable", by.y = "sample")


## DM1
count.matrix.glmm <- CREART::get_combined_matrix("/blackhole/alessia/circzi/checkCircRNAnormalizationdistribution/realdata/DM1/", 
                                                 select_methods = unique(colData$method))
## IPF
count.matrix.glmm <- CREART::get_combined_matrix("/blackhole/alessia/circzi/checkCircRNAnormalizationdistribution/realdata/IPF/analyses/", 
                                                 select_methods = unique(colData$method))
## ALZ
count.matrix.glmm <- CREART::get_combined_matrix("/blackhole/alessia/circzi/checkCircRNAnormalizationdistribution/realdata/ALZ/", 
                                                 select_methods = methods_det)

glmm.long = melt(count.matrix.glmm)
glmm.long$method = sub(".*\\.","",glmm.long$Var2)
glmm.long$sample = sub("\\..*", "", glmm.long$Var2)
glmm.wide = dcast(glmm.long, method+sample~Var1, value.var = "value", fill=0)
colnames(glmm.wide) = c("MethodID", "SampleID", colnames(glmm.wide)[-c(1,2)])
head(glmm.wide[,c(1:8)])
se.glmm <- SummarizedExperiment(assays = list(counts = count.matrix.glmm),
                                colData = colData)
eset.glmm <- ExpressionSet(assay(se.glmm),
                           AnnotatedDataFrame(as.data.frame(colData(se.glmm))))


# load function to perform DE analysis ------
source("/blackhole/alessia/GLMM_article/robustness/DEscripts.R")

## useful function for zinb-wave weights 
computeExactWeights <- function (model, x) 
{
  mu <- getMu(model)
  pi <- getPi(model)
  theta <- getTheta(model)
  theta <- matrix(rep(theta, each = ncol(x)), ncol = nrow(x))
  nb_part <- dnbinom(t(x), size = theta, mu = mu)
  zinb_part <- pi * ( t(x) == 0 ) + (1 - pi) *  nb_part
  zinbwg <- ( (1 - pi) * nb_part ) / zinb_part
  zinbwg <- t(zinbwg)
  zinbwg[x > 0] <- 1
  zinbwg[zinbwg < 1e-15] <- 1e-15
  rownames(zinbwg) = rownames(x)
  zinbwg
}

# specify the algorithms to be compared -------

algos <- list("DESeq2"= runDESeq2,
              #"DESeq2-ZI"=runDESeq2.ZI,
              #"DESeq2-glmGamPoi"= runDESeq2_gampoi,
              # "DESeq2-apeglm"=runDESeq2.apeglm, 
              #"DESeq2-ZINBWave"= DESeq_zinbweights,
              #"edgeR"=runEdgeR,
              "edgeR-robust"=runEdgeRRobust,
              "voom" = runVoom
              #"edgeR-ZINBWave"= edgeR_zinbweights
              #"circMeta" = runPois.ztest
) #"EBSeq"=runEBSeq)


namesAlgos <- names(algos)
names(namesAlgos) <- namesAlgos

# simulation results -----------

# run benchmark of NB vs GLMM models across evaluation and verification datasets 
# (using unbalanced verification and evalutation datasets) ----

resTest <- list()
resHeldout <- list()
lfcTest <- list()
lfcHeldout <- list()
ratemodel <- list()
resTestGLMM_NB <- list()
resHeldoutGLMM_NB <- list()
lfcTestGLMM_NB <- list()
lfcHeldoutGLMM_NB <- list()
resTestGLMM_ZINB <- list()
resHeldoutGLMM_ZINB <- list()
lfcTestGLMM_ZINB <- list()
lfcHeldoutGLMM_ZINB <- list()

set.seed(12388)
library("future.apply")

plan(multisession, workers = 3)

# res <- future_lapply(X = e, future.seed = T, FUN =  function(x) { 
#   for(i in pkgs) { library(i, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE, character.only=TRUE) }
#   library(dplyr)
#   library(data.table)
#   library(plyr)

## define BIOCPARALLEL parameters
hosts <- c("grigri", "grigri", "anhal")
param <- SnowParam(workers = hosts, type = "SOCK")
param

## summary stats of simulated datasets
stats_sim <- function(matrix){
  
  # matrix : count.matrices for one verification and evaluation
  
  n.circ <- nrow(matrix)
  zero <- sum(matrix==0)/length(matrix)
  summarydata <- data.frame(zero,
                            n.circ)
  
  names(summarydata) <- c("perc.zeros", "n.circular")
  
  return(summarydata)
}

library(future.apply)
resGLM_ALZ_30rep <- future_lapply(X = e, FUN =  function(x) {
  # x=e$findcirc
  res = bplapply(1:30, function(i) {   
    
    cat(i," ")
    # i = 11
    
    testSet <- as.character(randomSubsets[i,c(1:6)])
    heldOutSet <- as.character(randomSubsets[i,-c(1:6)])
    
    eTest <- x[,testSet]
    
    summary_stat_test <- stats_sim(matrix = as.matrix(exprs(eTest)))
    
    eHeldout <- x[,heldOutSet]
    
    summary_stat_heldout <- stats_sim(as.matrix(exprs(eHeldout)))
    
    resIDx = intersect(rownames(exprs(eTest)), rownames(exprs(eHeldout)))
    
    eTest = eTest[resIDx,]
    eHeldout = eHeldout[resIDx,]
    
    resTest0 <- lapply(namesAlgos, function(n) algos[[n]](e=eTest))
    
    resHeldout0 <- lapply(namesAlgos, function(n) algos[[n]](e=eHeldout))
    
    resTest <- as.data.frame(c(lapply(resTest0, function(z) z$padj[resIDx])))
    resHeldout <- as.data.frame(c(lapply(resHeldout0, function(z) z$padj[resIDx])))
    lfcTest <- as.data.frame(c(lapply(resTest0, function(z) z$beta[resIDx])))
    lfcHeldout <- as.data.frame(c(lapply(resHeldout0, function(z) z$beta[resIDx])))
    rownames(resTest) <- resIDx
    rownames(resHeldout) <- resIDx
    rownames(lfcTest) <- resIDx
    rownames(lfcHeldout) <- resIDx
    colnames(resTest) <- namesAlgos
    colnames(resHeldout) <- namesAlgos
    colnames(lfcTest) <- namesAlgos
    colnames(lfcHeldout) <- namesAlgos
    list(resTest=resTest,resHeldout=resHeldout,lfcTest=lfcTest,lfcHeldout=lfcHeldout,
         summary_stat_heldout=summary_stat_heldout,summary_stat_test=summary_stat_test
         # resTestGLMM_ZINB=resTestGLMM_ZINBNB,resHeldoutGLMM_ZINB=resHeldoutGLMM_ZINB,lfcTestGLMM_ZINB=lfcTestGLMM_ZINB,
         # lfcHeldoutGLMM_ZINB=lfcHeldoutGLMM_ZINB
    )
  }, BPPARAM =  MulticoreParam(workers = 5))
  # return(list(res=resMethods))
  res
})

resTes <- lapply(resGLM_ALZ_30rep, function(x) lapply(x, "[[", "resTest"))
resHeldout <- lapply(resGLM_ALZ_30rep, function(x) lapply(x, "[[", "resHeldout"))
lfcTest <- lapply(resGLM_ALZ_30rep, function(x) lapply(x, "[[", "lfcTest"))
lfcHeldout<- lapply(resGLM_ALZ_30rep, function(x) lapply(x, "[[", "lfcHeldout"))
summary_stat_test <-  lapply(resGLM_ALZ_30rep, function(x) bind_rows(lapply(x, "[[", "summary_stat_test")))
summary_stat_heldout <-  lapply(resGLM_ALZ_30rep, function(x) bind_rows(lapply(x, "[[", "summary_stat_heldout")))

save(resGLM_ALZ_30rep,resTes,resHeldout,lfcTest,lfcHeldout,namesAlgos,
     #/run/user/1001/gvfs/sftp:host=172.21.32.10
     file="/blackhole/alessia/GLMM_article/robustness/ALZ_sensitivityPrecision_GLM_30rep.RData")

colData = as.data.table(colData)
resGLMM_ALZ_30rep <- bplapply(1:30, function(i) {   
  
  cat(i," ")
  # i = 11
  
  testSetGLMM <- coldata[coldata$sample%in%as.character(randomSubsets[i,1:6]),]
  heldOutSetGLMM <- coldata[coldata$sample%in%as.character(randomSubsets[i,-c(1:6)]),]
  
  ## GLMM
  ## glmm-NB
  #test data
  pheno <- colData[colData$sample%in%testSetGLMM$sample,] %>% 
    dplyr::rename(SampleID = sample) %>% 
    dplyr::rename(MethodID = method)
  rownames(pheno) = pheno$sampl_id
  circularcounts <- count.matrix.glmm[, rownames(pheno)]
  
  summary_stat_test <- stats_sim(circularcounts)
  
  colnames(circularcounts) = sub("\\..*", "", colnames(circularcounts))
  testIDx = rownames(circularcounts)
  dge <- edgeR::DGEList(circularcounts)
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  offsets <- dge$samples$norm.factors # norm.factors per samples
  pheno$condition <- as.factor(pheno$condition)
  pheno$MethodID <- as.factor(pheno$MethodID)
  pheno$ln.lib.size <- offsets
  # circularcounts <- log(sweep(circularcounts,2,apply(circularcounts,2,mean),'/'))
  circularcounts <- t(circularcounts)
  allsamples_test <- cbind(pheno,circularcounts) %>% as.data.frame()
  
  
  # verification data
  pheno <- colData[colData$sample%in%heldOutSetGLMM$sample,] %>% 
    dplyr::rename(SampleID = sample) %>% 
    dplyr::rename(MethodID = method)
  rownames(pheno) = pheno$sampl_id
  
  circularcounts <- count.matrix.glmm[, rownames(pheno)]
  
  summary_stat_heldout <- stats_sim(circularcounts)
  
  
  colnames(circularcounts) = sub("\\..*", "", colnames(circularcounts))
  heldIDx = rownames(circularcounts)
  dge <- edgeR::DGEList(circularcounts)
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  offsets <- dge$samples$norm.factors # norm.factors per samples
  pheno$condition <- as.factor(pheno$condition)
  pheno$MethodID <- as.factor(pheno$MethodID)
  pheno$ln.lib.size <- offsets
  # circularcounts <- log(sweep(circularcounts,2,apply(circularcounts,2,mean),'/'))
  circularcounts <- t(circularcounts)
  allsamples_ver <- cbind(pheno,circularcounts) %>% as.data.frame()
  
  
  resIDx = intersect(testIDx, heldIDx)
  allsamples_test = allsamples_test[, colnames(allsamples_test)%in%c(colnames(allsamples_test)[1:4], resIDx)]
  allsamples_ver = allsamples_ver[, colnames(allsamples_ver)%in%c(colnames(allsamples_ver)[1:4], resIDx)]
  
  ## test GLMM-NB
  fit_GLMM_NB_test <- lapply(5:ncol(allsamples_test),
                             function(x){glmmTMB::glmmTMB(allsamples_test[,x] ~ condition + (1 | SampleID),
                                                          data=allsamples_test,
                                                          family=glmmTMB::nbinom2,
                                                          ziformula= ~0)})
  summaries <- lapply(fit_GLMM_NB_test, summary)
  pvalues <- as.numeric(unlist(lapply(summaries, function(x){stats::coef(x)$cond[2,4]})))
  pvalues[is.na(pvalues)] = 1
  padj = p.adjust(pvalues, method = "BH")
  # signif <- ifelse(pvalues < pval, 1, 0)
  # rateTMB <- mean(signif)
  lfc <- as.numeric(unlist(lapply(summaries, function(x){x$coefficients$cond["conditionnormal",1]}))) #conditionnormal
  resTestGLMM_NB = data.frame(GLMM = pvalues,
                              GLMMadj = padj,
                              row.names = colnames(allsamples_test)[5:ncol(allsamples_test)])
  lfcTestGLMM_NB <- data.frame(GLMM = lfc, row.names = colnames(allsamples_test)[5:ncol(allsamples_test)])
  
  ## verification
  fit_GLMM_NB_ver <- lapply(5:ncol(allsamples_ver),
                            function(x){glmmTMB::glmmTMB(allsamples_ver[,x] ~ condition + (1 | SampleID),
                                                         data=allsamples_ver,
                                                         family=glmmTMB::nbinom2,
                                                         ziformula= ~0)})
  summaries <- lapply(fit_GLMM_NB_ver, summary)
  pvalues <- as.numeric(unlist(lapply(summaries, function(x){stats::coef(x)$cond[2,4]})))
  pvalues[is.na(pvalues)] = 1
  padj = p.adjust(pvalues, method = "BH")
  # signif <- ifelse(pvalues < pval, 1, 0)
  # rateTMB <- mean(signif)
  lfc <- as.numeric(unlist(lapply(summaries, function(x){x$coefficients$cond["conditionnormal",1]})))
  resHeldoutGLMM_NB = data.frame(GLMM = pvalues,
                                 GLMMadj = padj,
                                 row.names = colnames(allsamples_ver)[5:ncol(allsamples_ver)])
  lfcHeldoutGLMM_NB <- data.frame(GLMM = lfc, row.names = colnames(allsamples_ver)[5:ncol(allsamples_ver)])
  
  
  resHeldoutGLMM_NB = data.frame(resHeldoutGLMM_NB[resIDx, ], row.names = resIDx)
  resTestGLMM_NB = data.frame(resTestGLMM_NB[resIDx, ], row.names = resIDx)
  lfcHeldoutGLMM_NB = data.frame(lfcHeldoutGLMM_NB[resIDx, ], row.names = resIDx)
  lfcTestGLMM_NB = data.frame(lfcTestGLMM_NB[resIDx, ], row.names = resIDx)
  
  
  list(resTestGLMM=resTestGLMM_NB,resHeldoutGLMM=resHeldoutGLMM_NB,
       lfcTestGLMM=lfcTestGLMM_NB,lfcHeldoutGLMM=lfcHeldoutGLMM_NB,
       summary_stat_heldout_GLMM=summary_stat_heldout,summary_stat_test_GLMM=summary_stat_test
       # resTestGLMM_ZINB=resTestGLMM_ZINBNB,resHeldoutGLMM_ZINB=resHeldoutGLMM_ZINB,lfcTestGLMM_ZINB=lfcTestGLMM_ZINB,
       # lfcHeldoutGLMM_ZINB=lfcHeldoutGLMM_ZINB
  )
}, BPPARAM =  MulticoreParam(workers = 5))

resTestGLMM <- lapply(resGLMM_ALZ_30rep, "[[", "resTestGLMM") #lapply(res, function(x) lapply(x, "[[", "resTest"))
resHeldoutGLMM <- lapply(resGLMM_ALZ_30rep, "[[", "resHeldoutGLMM") #lapply(res, function(x) lapply(x, "[[", "resHeldout"))
lfcTestGLMM <- lapply(resGLMM_ALZ_30rep, "[[", "lfcTestGLMM") #lapply(res, function(x) lapply(x, "[[", "lfcTest"))
lfcHeldoutGLMM<- lapply(resGLMM_ALZ_30rep, "[[", "lfcHeldoutGLMM") #lapply(res, function(x) lapply(x, "[[", "lfcHeldout"))
summary_stat_test_GLMM <- bind_rows(lapply(resGLMM_ALZ_30rep, "[[", "summary_stat_test_GLMM"))
summary_stat_heldout_GLMM <- bind_rows(lapply(resGLMM_ALZ_30rep, "[[", "summary_stat_heldout_GLMM"))

save(resGLMM_ALZ_30rep, resTestGLMM, resHeldoutGLMM,lfcTestGLMM,lfcHeldoutGLMM,
     summary_stat_test_GLMM, summary_stat_heldout_GLMM,
     file="/blackhole/alessia/GLMM_article/robustness/ALZ_sensitivityPrecision_GLMM_30rep.RData")

load("/blackhole/alessia/GLMM_article/robustness/ALZ_sensitivityPrecision_GLM_30rep.RData")
     # "/run/user/1001/gvfs/sftp:host=172.21.32.10/blackhole/alessia/GLMM_article/robustness/ALZ_sensitivityPrecision_GLM_30rep.RData")

library(tidyverse)
res_GLM_min = lapply(list(resTes, resHeldout), function(res){
  DT = melt(rbindlist(lapply(res, function(x){
    #x$circ_id = rownames(x)
    rbindlist(lapply(x, setDT, keep.rownames = TRUE),  idcol = c("Rep","rn"))
    }), idcol = "det.method"), id.vars = c("Rep", "det.method", "rn")) 
  minPadj = aggregate(value ~ variable + rn + Rep, DT, function(x) min(x))
  return(minPadj)
})
  
res_GLM_min_rep_Test = list()
for(i in 1:30){
  # i=1
  res_GLM_min_rep_Test[[i]] = res_GLM_min$resTest %>% dplyr::filter(Rep%in%i) %>% dcast(rn~variable,value.var="value")
  }
res_GLM_min_rep_Heldout = list()
for(i in 1:30){
  # i=1
  res_GLM_min_rep_Heldout[[i]] = res_GLM_min$resHeldout %>% dplyr::filter(Rep%in%i) %>% dcast(rn~variable,value.var="value")
}
res_GLM_GLMM_rep_Test = list()
for(r in 1:30){
  #r=1
  res_GLM_GLMM_rep_Test[[r]] = merge(res_GLM_min_rep_Test[[r]], resTestGLMM[[r]], by="rn",all=T)
}
res_GLM_GLMM_rep_Heldout = list()
for(r in 1:30){
  #r=1
  res_GLM_GLMM_rep_Heldout[[r]] = merge(res_GLM_min_rep_Heldout[[r]], resHeldoutGLMM[[r]], by="rn",all=T)
}
#save(res_GLM_GLMM_rep_Test, res_GLM_GLMM_rep_Heldout, file = "/blackhole/alessia/GLMM_article/robustness/robustness_30rep_GLM_GLMM.RData")
names(res_GLM_min) = c("resTest", "resHeldout")
#write.csv(minPadj, "/blackhole/alessia/GLMM_article/robustness/minPadj_30Rep_ALZ.csv")
res_GLM_min = read.csv("/blackhole/alessia/GLMM_article/robustness/minPadj_30Rep_ALZ.csv")
load("/blackhole/alessia/GLMM_article/robustness/ALZ_sensitivityPrecision_GLMM_30rep.RData")

res_GLMM = lapply(list(resTestGLMM, resHeldoutGLMM), function(res){
  resGLMM = melt(rbindlist(lapply(res, setDT, keep.rownames = TRUE),  idcol = c("Rep","rn")), id.vars = c("Rep", "rn")) 
  return(resGLMM)
})
names(res_GLMM) = c("resTest", "resHeldout")

res.concordance_30rep = rbind.fill(rbindlist(res_GLMM, idcol = "result"), rbindlist(res_GLM_min, idcol = "result"))
write.csv(res.concordance_30rep, "/blackhole/alessia/GLMM_article/robustness/GLM_GLMM_Padj_30Rep_ALZ.csv")


lfc_GLM_min = lapply(list(lfcTest, lfcHeldout), function(res){
  DT = melt(rbindlist(lapply(res, function(x){
    #x$circ_id = rownames(x)
    rbindlist(lapply(x, setDT, keep.rownames = TRUE),  idcol = c("Rep","rn"))
  }), idcol = "det.method"), id.vars = c("Rep", "det.method", "rn")) 
  minPadj = aggregate(value ~ variable + rn + Rep, DT, function(x) min(x))
  return(minPadj)
})
names(lfc_GLM_min) = c("resTest", "resHeldout")
lfc_GLM_min_rep_Test = list()
for(i in 1:30){
  # i=1
  lfc_GLM_min_rep_Test[[i]] = lfc_GLM_min$resTest %>% dplyr::filter(Rep%in%i) %>% dcast(rn~variable,value.var="value")
}
lfc_GLM_min_rep_Heldout = list()
for(i in 1:30){
  # i=1
  lfc_GLM_min_rep_Heldout[[i]] = lfc_GLM_min$resHeldout %>% dplyr::filter(Rep%in%i) %>% dcast(rn~variable,value.var="value")
}
lfc_GLM_GLMM_rep_Test = list()
for(r in 1:30){
  #r=1
  lfc_GLM_GLMM_rep_Test[[r]] = merge(lfc_GLM_min_rep_Test[[r]], lfcTestGLMM[[r]], by.x="rn", by.y="row.names", all=T)
}
lfc_GLM_GLMM_rep_Test = lapply(lfc_GLM_GLMM_rep_Test, function(x) {
  names(x)=c("rn","DESeq2", "edgeR-robust","voom","GLMM")
  return(x)
})

lfc_GLM_GLMM_rep_Heldout = list()
for(r in 1:30){
  #r=1
  lfc_GLM_GLMM_rep_Heldout[[r]] = merge(lfc_GLM_min_rep_Heldout[[r]], lfcHeldoutGLMM[[r]], by.x="rn", by.y="row.names", all=T)
}
lfc_GLM_GLMM_rep_Heldout = lapply(lfc_GLM_GLMM_rep_Heldout, function(x) {
  names(x)=c("rn","DESeq2", "edgeR-robust","voom","GLMM")
return(x)
})
#save(lfc_GLM_GLMM_rep_Test, lfc_GLM_GLMM_rep_Heldout, file = "/blackhole/alessia/GLMM_article/robustness/robustness_30rep_GLM_GLMM_lfc.RData")

# type I error control -------------

randomSubsets <- read.table("/blackhole/alessia/CircModel/type_I_error_control/random_shuffle_ALZ.txt",strings=FALSE)
randomSubsets <- read.table("/blackhole/alessia/CircModel/type_I_error_control/random_shuffle_TALL.txt",strings=FALSE)

nreps <- 30
set.seed(12388)
library("future.apply")

# plan(multisession, workers = 2)
## define BIOCPARALLEL parameters
hosts <- c(rep("grigri", 8), "threesum", "anhal")
param <- SnowParam(workers = hosts, type = "SOCK")
param

ccp2 = as.data.frame(ccp2)
# ccp2$circ_id = rownames(ccp2)
resGLM_ALZ_30rep_tIe <- future_lapply(X = e, FUN =  function(x) {
  # x=e$findcirc
  res = bplapply(1:30, function(i) {   
    
    cat(i," ")
    # i = 1
    
    testSet <- as.character(randomSubsets[i,c(1:6)])

    eTest <- x[,testSet]
    
    summary_stat_test <- stats_sim(matrix = as.matrix(exprs(eTest)))
    
    resTest0 <- lapply(namesAlgos, function(n) algos[[n]](e=eTest))
    
    resTest <- as.data.frame(c(lapply(resTest0, function(z) z$padj[resIDx])))
    resTest_raw = lapply(resTest0, function(x) {
      x$pvals[is.na(x$pvals)] = 1
      return(x)}
      )
    signif1 <- lapply(resTest_raw, function(x) ifelse(x$pvals <= 0.1, 1, 0))
    sigific05 <- lapply(resTest_raw, function(x) ifelse(x$pvals <= 0.05, 1, 0))
    ratemodel <- list(res.1=c(lapply(signif1, function(x) mean(x))),
                      res.05=c(lapply(sigific05, function(x) mean(x))))
    ratemodel.dt = rbindlist(ratemodel)
    ratemodel.dt$alpha = c(0.1, 0.05)
    
    list(resTest=resTest,
         rate = ratemodel.dt
    )
  }, BPPARAM =  MulticoreParam(workers = 5))
  res
})

resTes <- lapply(resGLM_ALZ_30rep_tIe, function(x) lapply(x, "[[", "resTest"))
rate <- lapply(resGLM_ALZ_30rep_tIe, function(x) lapply(x, "[[", "rate"))

library(tidyverse)

DT = melt(rbindlist(lapply(resTes, function(x){
  #x$circ_id = rownames(x)
  rbindlist(lapply(x, setDT, keep.rownames = TRUE),  idcol = c("Rep","rn"))
}), idcol = "det.method"), id.vars = c("Rep", "det.method", "rn")) 
minPadj = aggregate(value ~ variable + rn + Rep, DT, function(x) min(x))
res_GLM_min = minPadj

res_GLM_min_rep_Test = list()
for(i in 1:30){
  # i=1
  res_GLM_min_rep_Test[[i]] = res_GLM_min %>% dplyr::filter(Rep%in%i) %>% dcast(rn~variable,value.var="value")
}
signif1 <- lapply(res_GLM_min_rep_Test, function(x) apply(x[,-1],2,function(c) ifelse(c <= 0.1, 1, 0)))
sigific05 <- lapply(res_GLM_min_rep_Test, function(x) apply(x[,-1],2,function(c) ifelse(c <= 0.05, 1, 0)))
ratemodel <- list(res.1=lapply(signif1, function(x) apply(x,2,function(c) mean(c))),
                  res.05=lapply(sigific05, function(x) apply(x,2,function(c) mean(c))))

ratemodel.dt = rbindlist(list(res.1 = rbindlist(lapply(ratemodel$res.1, function(x) {
  dt = setDT(as.data.frame(x), keep.rownames = T)
  return(melt(dt))}), idcol = "rep"),
  res.05 = rbindlist(lapply(ratemodel$res.05, function(x) {
    
        dt = setDT(as.data.frame(x), keep.rownames = T)
        return(melt(dt))}), idcol = "rep")),idcol="alpha") %>% dplyr::select(alpha, rep, rn, value) %>% dplyr::rename("method"="rn",
                                                                                                                      "rate"="value")

save(resGLM_ALZ_30rep_tIe,resTes,rate,namesAlgos,ratemodel.dt,
     #/run/user/1001/gvfs/sftp:host=172.21.32.10
     file="/blackhole/alessia/GLMM_article/type_I_error_control/ALZ_GLM_30rep.RData")

resGLMM_ALZ_30rep_tIe <- bplapply(1:30, function(i) {   
  
  cat(i," ")
  # i = 11
  
  testSetGLMM <- coldata[coldata$sample%in%as.character(randomSubsets[i,1:6]),]
  pheno <- colData[colData$sample%in%testSetGLMM$sample,] %>% 
    dplyr::rename(SampleID = sample) %>% 
    dplyr::rename(MethodID = method)
  circularcounts <- count.matrix.glmm[, rownames(pheno)]
  
  summary_stat_test <- stats_sim(circularcounts)
  
  colnames(circularcounts) = sub("\\..*", "", colnames(circularcounts))
  testIDx = rownames(circularcounts)
  dge <- edgeR::DGEList(circularcounts)
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  offsets <- dge$samples$norm.factors # norm.factors per samples
  pheno$condition <- as.factor(pheno$condition)
  pheno$MethodID <- as.factor(pheno$MethodID)
  pheno$ln.lib.size <- offsets
  # circularcounts <- log(sweep(circularcounts,2,apply(circularcounts,2,mean),'/'))
  circularcounts <- t(circularcounts)
  allsamples_test <- cbind(pheno,circularcounts) %>% as.data.frame()

    ## test GLMM-NB
  fit_GLMM_NB_test <- lapply(5:ncol(allsamples_test),
                             function(x){glmmTMB::glmmTMB(allsamples_test[,x] ~ condition + (1 | SampleID),
                                                          data=allsamples_test,
                                                          family=glmmTMB::nbinom2,
                                                          ziformula= ~0)})
  summaries <- lapply(fit_GLMM_NB_test, summary)
  pvalues <- as.numeric(unlist(lapply(summaries, function(x){stats::coef(x)$cond[2,4]})))
  pvalues[is.na(pvalues)] = 1
  padj = p.adjust(pvalues, method = "BH")
  signif1 <- ifelse(pvalues <= 0.1, 1, 0)
  signif05 <- ifelse(pvalues <= 0.05, 1, 0)
  
  rateTMB <- data.frame(alpha=c("res.1", "res.05"),
                        method = c("GLMM","GLMM"),
                        rate = c(mean(signif1), mean(signif05)))
  
  resTestGLMM_NB = data.frame(GLMM = pvalues,
                              GLMMadj = padj,
                              row.names = colnames(allsamples_test)[5:ncol(allsamples_test)])
  
  resTestGLMM_NB = data.frame(resTestGLMM_NB[testIDx, ], row.names = testIDx)

  
  list(resTestGLMM = resTestGLMM_NB,
       rate = rateTMB,
       summary_stat_test_GLMM = summary_stat_test
  )
}, BPPARAM =  MulticoreParam(workers = 5))

resTestGLMM_tIe <- lapply(resGLMM_ALZ_30rep_tIe, "[[", "resTestGLMM") #lapply(res, function(x) lapply(x, "[[", "resTest"))
rateGLMM_tIe <- lapply(resGLMM_ALZ_30rep_tIe, "[[", "rate") #lapply(res, function(x) lapply(x, "[[", "resHeldout"))
summary_stat_test_GLMM_tIe <- bind_rows(lapply(resGLMM_ALZ_30rep_tIe, "[[", "summary_stat_test_GLMM"))

save(resGLMM_ALZ_30rep_tIe, resTestGLMM_tIe, rateGLMM_tIe,summary_stat_test_GLMM_tIe,
     file="/blackhole/alessia/GLMM_article/type_I_error_control/ALZ_typeIerror_GLMM_30rep.RData")
