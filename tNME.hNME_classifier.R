#####################################################
#### Select Tumor signature gene with SMC cohort ####
#####################################################
DE_TN <- readRDS("SMC_RNA_seq_TN222.rds")
TN_group <- ifelse(str_detect(colnames(DE_TN), "T"), "Tumor", "Normal" )
tmp <- DGEList(counts=DE_TN, group = factor(TN_group))
keep <- filterByExpr(tmp, group=TN_group)
tmp<-tmp[keep,,keep.lib.sizes=F]
tmp <- calcNormFactors(tmp) 

tmp2 <- tmp
design <- model.matrix(~TN_group)
tmp2<-estimateGLMCommonDisp(tmp2,design)
tmp2<-estimateGLMTrendedDisp(tmp2,design)
tmp2<-estimateGLMTagwiseDisp(tmp2,design)
t_fit<-glmQLFit(tmp2,design, robust = T)
tweb<-glmQLFTest(t_fit)

## Select Tumor signature gene
res<-tweb$table
res$fdr <- p.adjust(res$PValue, method="fdr")
cpm = 3
T_up <- test_t[ test_t$fdr < 0.05 & test_t$logCPM > cpm & test_t$logFC > 1,] %>% dplyr::arrange(desc(logFC))
T_top_20 <- rownames(T_up)[1:20]

## Score normal samples and divide them into groups
tmp_log <- cpm(tmp, log=T)
tmp_log <- tmp_log[,colnames(tmp_log) %in% colnames(tmp_log)[str_detect(colnames(tmp_log),"N")]]
T_top_list <- list(T_top_20)
names(T_top_list) <- c("T_top_20")
gsva_normal <- gsva(tmp_log, T_top_list, mx.diff=F, verbose=F, parallel.sz=1, kcdf="Poisson", method="ssgsea", ssgsea.norm=FALSE)
gsva_normal <- t(gsva_normal) %>% as.data.frame()
gsva_normal$group <- ifelse(gsva_normal$T_top_20 > mean(gsva_normal$T_top_20), "tNME", "hNME")

#################################################################
#### Validate grouping with Multi-center cohort (except SMC) ####
#################################################################
YS <- readRDS("Others_RNA_seq_N162.rds")
DE_TN <- YS[rownames(YS) %in% rownames(tmp),] 

YS_log <- cpm(DE_TN, log=T)
YS_log <- YS_log[rownames(tmp),]
YS_log <- as.matrix(YS_log)

YS_gsva <- gsva(YS_log, T_top_list, mx.diff=F, verbose=T, parallel.sz=1, kcdf="Poisson", method="ssgsea",ssgsea.norm=FALSE)
YS_gsva <- t(YS_gsva) %>% as.data.frame()
YS_gsva$group <- ifelse(YS_gsva$T_top_20 > mean(gsva_normal$T_top_20), "tNME", "hNME")
