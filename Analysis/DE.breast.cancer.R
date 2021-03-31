source('client-side/code/DEBoost.R')

################################################################################################################
# Prepare liver data, here we only want female samples
################################################################################################################
load('server-side//RData//Liver.RData')
Female.sample                 <- sample.meta.df$sample.id[sample.meta.df$gender == 'Female'] %>% as.character()
REF.log2.read.count.matrix    <- log2.read.count.matrix[,Female.sample]
REF.log2.tpm.matrix           <- log2.tpm.matrix[,Female.sample]


load('server-side/RData//Breast Invasive Carcinoma.RData')
load('~/Project/Cancer2CellLine/server-side/RData/MET500.RData')
load('client-side/output/Select.pure.sample.breast.cancer.R.output/Select.pure.sample.breast.cancer.RData')


### Basal-like subtype #########

###### Prepare primary cancer data
PRI.log2.read.count.matrix <- log2.read.count.matrix[,pure.PRI.breast.cancer.Basal.sample]
PRI.log2.tpm.matrix       <- log2.tpm.matrix[,pure.PRI.breast.cancer.Basal.sample]

#### Prepare metastatic cancer data
MET.log2.read.count.matrix <- MET500.log2.read.count.matrix[,pure.MET.breast.cancer.Basal.sample]
MET.log2.tpm.matrix       <- MET500.log2.tpm.matrix[,pure.MET.breast.cancer.Basal.sample]

##### DE analysis 
BRCA.Basal.DE.rs <- perform.DE.analysis.between.primary.and.metastatic.cancer(
  PRI.log2.tpm.matrix = PRI.log2.tpm.matrix, PRI.log2.read.count.matrix = PRI.log2.read.count.matrix,
  MET.log2.tpm.matrix = MET.log2.tpm.matrix, MET.log2.read.count.matrix = MET.log2.read.count.matrix,
  REF.log2.tpm.matrix = REF.log2.tpm.matrix, REF.log2.read.count.matrix = REF.log2.read.count.matrix,
  TCGA.best.cell.line = TCGA.best.cell.line.Basal, MET500.best.cell.line = MET500.best.cell.line.Basal
  
)


### LumB subtype #########

###### Prepare primary cancer data
PRI.log2.read.count.matrix <- log2.read.count.matrix[,c(pure.PRI.breast.cancer.LumB.sample)]
PRI.log2.tpm.matrix        <- log2.tpm.matrix[,c(pure.PRI.breast.cancer.LumB.sample)]

#### Prepare metastatic cancer data
MET.log2.read.count.matrix <- MET500.log2.read.count.matrix[,c(pure.MET.breast.cancer.LumB.sample)]
MET.log2.tpm.matrix       <- MET500.log2.tpm.matrix[,c(pure.MET.breast.cancer.LumB.sample)]

##### DE analysis 
BRCA.LumB.DE.rs <- perform.DE.analysis.between.primary.and.metastatic.cancer(
  PRI.log2.tpm.matrix = PRI.log2.tpm.matrix, PRI.log2.read.count.matrix = PRI.log2.read.count.matrix,
  MET.log2.tpm.matrix = MET.log2.tpm.matrix, MET.log2.read.count.matrix = MET.log2.read.count.matrix,
  REF.log2.tpm.matrix = REF.log2.tpm.matrix, REF.log2.read.count.matrix = REF.log2.read.count.matrix,
  TCGA.best.cell.line = TCGA.best.cell.line.LumB, MET500.best.cell.line = MET500.best.cell.line.LumB
  
)


### Her2 subtype #########

###### Prepare primary cancer data
PRI.log2.read.count.matrix <- log2.read.count.matrix[,pure.PRI.breast.cancer.Her2.sample]
PRI.log2.tpm.matrix        <- log2.tpm.matrix[,pure.PRI.breast.cancer.Her2.sample]

#### Prepare metastatic cancer data
MET.log2.read.count.matrix <- MET500.log2.read.count.matrix[,pure.MET.breast.cancer.Her2.sample]
MET.log2.tpm.matrix       <- MET500.log2.tpm.matrix[,pure.MET.breast.cancer.Her2.sample]

##### DE analysis 
BRCA.Her2.DE.rs <- perform.DE.analysis.between.primary.and.metastatic.cancer(
  PRI.log2.tpm.matrix = PRI.log2.tpm.matrix, PRI.log2.read.count.matrix = PRI.log2.read.count.matrix,
  MET.log2.tpm.matrix = MET.log2.tpm.matrix, MET.log2.read.count.matrix = MET.log2.read.count.matrix,
  REF.log2.tpm.matrix = REF.log2.tpm.matrix, REF.log2.read.count.matrix = REF.log2.read.count.matrix,
  TCGA.best.cell.line = TCGA.best.cell.line.Her2, MET500.best.cell.line = MET500.best.cell.line.Her2
  
)


save(file = 'client-side/output/DE.breast.cancer.R.output/DE.breast.cancer.RData',list=c('BRCA.Basal.DE.rs','BRCA.Her2.DE.rs','BRCA.LumB.DE.rs'))



