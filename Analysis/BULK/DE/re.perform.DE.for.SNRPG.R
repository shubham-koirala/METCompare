source('client-side/code/DEBoost.R')

################################################################################################################
# Prepare liver data, here we only want male samples
################################################################################################################
load('server-side//RData//Liver.RData')
Female.sample                 <- sample.meta.df$sample.id[sample.meta.df$gender == 'Female'] %>% as.character()
REF.log2.read.count.matrix    <- log2.read.count.matrix[,Female.sample]
REF.log2.tpm.matrix           <- log2.tpm.matrix[,Female.sample]


load('server-side/RData//Breast Invasive Carcinoma.RData')
load('~/Project/Cancer2CellLine/server-side/RData/MET500.RData')
load('client-side/output/Select.pure.sample.breast.cancer.R.output/Select.pure.sample.breast.cancer.RData')

PRI.log2.read.count.matrix <- log2.read.count.matrix[,c(pure.PRI.breast.cancer.LumB.sample)]
PRI.log2.tpm.matrix        <- log2.tpm.matrix[,c(pure.PRI.breast.cancer.LumB.sample)]

#### Prepare metastatic cancer data
pure.MET.breast.cancer.LumB.sample <- setdiff(pure.MET.breast.cancer.LumB.sample,'SRR4307452')
MET.log2.read.count.matrix         <- MET500.log2.read.count.matrix[,c(pure.MET.breast.cancer.LumB.sample)]
MET.log2.tpm.matrix                <- MET500.log2.tpm.matrix[,c(pure.MET.breast.cancer.LumB.sample)]

##### DE analysis 
BRCA.LumB.DE.rs.without.outlier <- perform.DE.analysis.between.primary.and.metastatic.cancer(
  PRI.log2.tpm.matrix = PRI.log2.tpm.matrix, PRI.log2.read.count.matrix = PRI.log2.read.count.matrix,
  MET.log2.tpm.matrix = MET.log2.tpm.matrix, MET.log2.read.count.matrix = MET.log2.read.count.matrix,
  REF.log2.tpm.matrix = REF.log2.tpm.matrix, REF.log2.read.count.matrix = REF.log2.read.count.matrix,
  TCGA.best.cell.line = TCGA.best.cell.line.LumB, MET500.best.cell.line = MET500.best.cell.line.LumB
  
)


SNRPG <- 'ENSG00000143977'

BRCA.LumB.DE.rs.without.outlier$deseq2.M.vs.P.res[SNRPG,]


