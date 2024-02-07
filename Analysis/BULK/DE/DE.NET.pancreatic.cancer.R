require(AnnotationDbi)
library(org.Hs.eg.db)
source('client-side/code/DEBoost.R')

################################################################################################################
# Prepare liver data, here we  want both male and female samples
################################################################################################################
load('server-side//RData//Liver.RData')
REF.log2.tpm.matrix          <- log2.tpm.matrix
REF.log2.read.count.matrix   <- log2.read.count.matrix


load('client-side/output/Select.pure.sample.NET.pancreatic.cancer.R.output/Select.pure.sample.NET.pancreatic.cancer.RData')
load('server-side/RData/GEP.NET.RData')
################################################################################################################
# Prepare primary cancer data
################################################################################################################
PRI.log2.read.count.matrix <- GEP.NET.log2.read.count.matrix[,pure.PRI.NET.pancreatic.cancer.sample]
PRI.log2.tpm.matrix       <- GEP.NET.log2.tpm.matrix[,pure.PRI.NET.pancreatic.cancer.sample]


################################################################################################################
# Prepare metastatic cancer data
################################################################################################################
MET.log2.read.count.matrix <- GEP.NET.log2.read.count.matrix[,pure.MET.NET.pancreatic.cancer.sample]
MET.log2.tpm.matrix       <- GEP.NET.log2.tpm.matrix[,pure.MET.NET.pancreatic.cancer.sample]



NET.PAAD.DE.rs <- perform.DE.analysis.between.primary.and.metastatic.cancer(
  PRI.log2.tpm.matrix = PRI.log2.tpm.matrix, PRI.log2.read.count.matrix = PRI.log2.read.count.matrix,
  MET.log2.tpm.matrix = MET.log2.tpm.matrix, MET.log2.read.count.matrix = MET.log2.read.count.matrix,
  REF.log2.tpm.matrix = REF.log2.tpm.matrix, REF.log2.read.count.matrix = REF.log2.read.count.matrix,
  TCGA.best.cell.line = TCGA.best.cell.line, MET500.best.cell.line = MET500.best.cell.line
  
)

save(file='client-side/output/DE.NET.pancreatic.cancer.R.output/DE.NET.pancreatic.cancer.RData',list=c('NET.PAAD.DE.rs'))





