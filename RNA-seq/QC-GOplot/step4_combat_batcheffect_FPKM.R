## remove batch effect before clustering analysis
#mat <- fread("results/monkey-RNASeq-cqnnormalized-FPKM-log2-rmd.txt", header = T)
#rownames(mat) <- mat$id
#mat$id <- NULL
mat <- read.xlsx("results/monkey-RNASeq-standward-FPKM-log2-rmd-exonbygene.xlsx", rowNames = T, colNames = T)

se <- readRDS("data/monkey-exonbygene-summaried-expriment-rmd.rds")

pheno = data.frame(colData(se)[,1:3])
rownames(pheno) = pheno$id
pheno$batch = c(rep(c(2,3,3),4), rep(1,6))
batch = pheno$batch
modcombat = model.matrix(~1, data=pheno)
combat_edata = ComBat(dat=mat, batch=batch, mod=modcombat)

write.xlsx(combat_edata, "results/monkey-RNASeq-standward-combat-removebatch-FPKM-log2-rmd.xlsx", 
            row.names = T, col.names = T)
