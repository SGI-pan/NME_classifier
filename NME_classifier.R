# processing with study cohort
DE_TN <- readRDS("Total_raw.rds")
TN_group <- ifelse(str_detect(colnames(DE_TN), "T"), "Tumor", "Normal" )
tmp <- DGEList(counts=DE_TN, group = factor(TN_group))
keep <- filterByExpr(tmp, group=TN_group)
tmp<-tmp[keep,,keep.lib.sizes=F]
tmp <- calcNormFactors(tmp) 

tmp2 <- tmp
design <- model.matrix(~TN_group)
tmp2<-estimateDisp(tmp2,design)
t_fit<-glmQLFit(tmp2,design, robust = T)
tweb<-glmQLFTest(t_fit)

## Select Tumor signature gene
res<-tweb$table
res$fdr <- p.adjust(res$PValue, method="fdr")
T_up <- test_t[ test_t$fdr < 0.05 & test_t$logCPM > 3 & test_t$logFC > 1,] %>% dplyr::arrange(desc(logFC))
TSG <- rownames(T_up)[1:28]

## Scoring NME samples with TSG gene set and dividing them into groups (TSM or HM)
tmp_log <- cpm(tmp, log=T)
tmp_log <- tmp_log[,colnames(tmp_log) %in% colnames(tmp_log)[str_detect(colnames(tmp_log),"N")]]
T_top_list <- list(TSG)
names(T_top_list) <- c("TSG")
gsva_normal <- gsva(tmp_log, T_top_list, mx.diff=F, verbose=F, parallel.sz=1, kcdf="Poisson", method="ssgsea", ssgsea.norm=FALSE)
gsva_normal <- t(gsva_normal) %>% as.data.frame()
gsva_normal$group <- ifelse(gsva_normal$TSG > mean(gsva_normal$TSG), "TSM", "HM")
