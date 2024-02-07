source('client-side/code/DEBoost.R')

################################################################################################################
# Prepare liver data, here we only want male samples
################################################################################################################
load('server-side//RData//Liver.RData')
Male.sample                   <- sample.meta.df$sample.id[sample.meta.df$gender == 'Male'] %>% as.character()
REF.log2.read.count.matrix    <- log2.read.count.matrix[,Male.sample]
REF.log2.tpm.matrix           <- log2.tpm.matrix[,Male.sample]


load('client-side/output/Select.pure.sample.prostate.cancer.R.output/Select.pure.sample.prostate.cancer.RData')

################################################################################################################
# Prepare primary cancer data
################################################################################################################
load('server-side/RData//Prostate Adenocarcinoma.RData')
PRI.log2.read.count.matrix <- log2.read.count.matrix[,pure.PRI.prostate.cancer.sample]
PRI.log2.tpm.matrix        <- log2.tpm.matrix[,pure.PRI.prostate.cancer.sample]


################################################################################################################
# Prepare metastatic cancer data
################################################################################################################
load('~/Project/Cancer2CellLine/server-side/RData/MET500.RData')
MET.log2.read.count.matrix <- MET500.log2.read.count.matrix[,pure.MET.prostate.cancer.sample]
MET.log2.tpm.matrix       <- MET500.log2.tpm.matrix[,pure.MET.prostate.cancer.sample]



PRAD.DE.rs <- perform.DE.analysis.between.primary.and.metastatic.cancer(
  PRI.log2.tpm.matrix = PRI.log2.tpm.matrix, PRI.log2.read.count.matrix = PRI.log2.read.count.matrix,
  MET.log2.tpm.matrix = MET.log2.tpm.matrix, MET.log2.read.count.matrix = MET.log2.read.count.matrix,
  REF.log2.tpm.matrix = REF.log2.tpm.matrix, REF.log2.read.count.matrix = REF.log2.read.count.matrix,
  TCGA.best.cell.line = TCGA.best.cell.line, MET500.best.cell.line = MET500.best.cell.line
)


save(file='client-side/output/DE.prostate.cancer.R.output/DE.prostate.cancer.RData',list=c('PRAD.DE.rs'))


