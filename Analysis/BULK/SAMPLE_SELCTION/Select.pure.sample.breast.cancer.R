
require(plyr)
require(dplyr)
require(genefu)
require(Rtsne)
require(ggplot2)
require(dplyr)
source('client-side/code/util.R')
load('server-side/RData/Breast Invasive Carcinoma.RData')


brca.tcga.meta <- read.delim("~/Project/BreastCancerMetaPotenial/client-side/Data/cBioPortal/brca_tcga_pan_can_atlas_2018/data_clinical_sample.txt", stringsAsFactors=FALSE,comment.char = '#',header = TRUE)
duc.sample     <- brca.tcga.meta$SAMPLE_ID[brca.tcga.meta$CANCER_TYPE_DETAILED == 'Breast Invasive Ductal Carcinoma'] %>% as.character()
duc.sample     <- intersect(duc.sample,colnames(log2.read.count.matrix))
TCGA.breast.cancer.polyA.sample           <- duc.sample
TCGA.breast.cancer.log2.read.count.matrix <- log2.read.count.matrix[,duc.sample]
TCGA.breast.cancer.log2.fpkm.matrix       <- log2.fpkm.matrix[,duc.sample]


#### Perform pam50 subtyping for TCGA polyA samples###########
pam50.gene.df               <- read.delim("~/Project/Cancer2CellLine/client-side/meta.data/pam50_entrez_gene_id_and_ensemble_gene_id.txt", stringsAsFactors=FALSE,comment.char = '#') # well, I borrow some information from the metastatic breast cancer evaluation project 
pam50.gene                  <- pam50.gene.df$ensemble.gene.id %>% as.character
colnames(pam50.gene.df)[1]  <- 'EntrezGene.ID'
colnames(pam50.gene.df)[2]  <- 'probe' #damn it, genefu package doesnot mentioning this!
pam50.gene.df$EntrezGene.ID <- as.character(pam50.gene.df$EntrezGene.ID)
pam50.gene.expr             <- TCGA.breast.cancer.log2.fpkm.matrix[pam50.gene.df$probe %>% as.character,TCGA.breast.cancer.polyA.sample] %>% t 
annot.matrix                <- pam50.gene.df[,1:2] %>% as.matrix
rownames(annot.matrix)      <- annot.matrix[,'probe']
pam50.subtype.rs            <- intrinsic.cluster.predict(sbt.model = pam50.robust,data = pam50.gene.expr,annot = annot.matrix,mapping=annot.matrix ,do.mapping = TRUE )


TCGA.breast.cancer.polyA.LumB.sample    <- TCGA.breast.cancer.polyA.sample[pam50.subtype.rs$subtype == 'LumB']
TCGA.breast.cancer.polyA.Basal.sample   <- TCGA.breast.cancer.polyA.sample[pam50.subtype.rs$subtype == 'Basal']
TCGA.breast.cancer.polyA.Her2.sample    <- TCGA.breast.cancer.polyA.sample[pam50.subtype.rs$subtype == 'Her2']
TCGA.breast.cancer.polyA.LumA.sample    <- TCGA.breast.cancer.polyA.sample[pam50.subtype.rs$subtype == 'LumA']
TCGA.breast.cancer.polyA.Normal.sample  <- TCGA.breast.cancer.polyA.sample[pam50.subtype.rs$subtype == 'Normal']




############## TCGA.sample - CCLE.cell.line correlation analysis ################
load('~/Project/Cancer2CellLine/server-side/RData/CCLE.RData')
CCLE.median                 <- apply(CCLE.log2.rpkm.matrix,1,median)
CCLE.expressed.gene         <- names(CCLE.median)[CCLE.median > 1]
tmp                         <- CCLE.log2.rpkm.matrix[CCLE.expressed.gene,]
tmp.rank                    <- apply(tmp,2,rank)
rank.mean                   <- apply(tmp.rank,1,mean)
rank.sd                     <- apply(tmp.rank,1,sd)
plot(x=rank.mean,y=rank.sd)
lowess(x=rank.mean,y=rank.sd) %>% lines(lwd=5,col='red')
CCLE.rna.seq.marker.gene.1000   <- names(sort(rank.sd,decreasing =TRUE))[1:1000]

#### TCGA #######
cor.cut.off                           <- 0.3
rank.cut.off                          <- 1010

rs.TCGA                               <- pick.out.cell.line(expr.of.samples = TCGA.breast.cancer.log2.fpkm.matrix[,TCGA.breast.cancer.polyA.Basal.sample],expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
expr.cor                              <- rs.TCGA$correlation.matrix[,rs.TCGA$best.cell.line]
cell.line.rank                        <- apply(rs.TCGA$correlation.matrix, 1, function(x) rank(x)[rs.TCGA$best.cell.line] )
pure.PRI.breast.cancer.Basal.sample   <- names(expr.cor)[(expr.cor > cor.cut.off)  & (cell.line.rank > rank.cut.off) ]
TCGA.best.cell.line.Basal             <- rs.TCGA$best.cell.line


rs.TCGA                               <- pick.out.cell.line(expr.of.samples = TCGA.breast.cancer.log2.fpkm.matrix[,TCGA.breast.cancer.polyA.Her2.sample],expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
expr.cor                              <- rs.TCGA$correlation.matrix[,rs.TCGA$best.cell.line]
cell.line.rank                        <- apply(rs.TCGA$correlation.matrix, 1, function(x) rank(x)[rs.TCGA$best.cell.line] )
pure.PRI.breast.cancer.Her2.sample    <- names(expr.cor)[(expr.cor > cor.cut.off)  & (cell.line.rank > rank.cut.off) ]
TCGA.best.cell.line.Her2              <- rs.TCGA$best.cell.line


rs.TCGA                               <- pick.out.cell.line(expr.of.samples = TCGA.breast.cancer.log2.fpkm.matrix[,TCGA.breast.cancer.polyA.LumB.sample],expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
expr.cor                              <- rs.TCGA$correlation.matrix[,rs.TCGA$best.cell.line]
cell.line.rank                        <- apply(rs.TCGA$correlation.matrix, 1, function(x) rank(x)[rs.TCGA$best.cell.line] )
pure.PRI.breast.cancer.LumB.sample    <- names(expr.cor)[(expr.cor > cor.cut.off)  & (cell.line.rank > rank.cut.off) ]
TCGA.best.cell.line.LumB               <- rs.TCGA$best.cell.line


rs.TCGA                               <- pick.out.cell.line(expr.of.samples = TCGA.breast.cancer.log2.fpkm.matrix[,TCGA.breast.cancer.polyA.LumA.sample],expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
expr.cor                              <- rs.TCGA$correlation.matrix[,rs.TCGA$best.cell.line]
cell.line.rank                        <- apply(rs.TCGA$correlation.matrix, 1, function(x) rank(x)[rs.TCGA$best.cell.line] )
pure.PRI.breast.cancer.LumA.sample    <- names(expr.cor)[(expr.cor > cor.cut.off)  & (cell.line.rank > rank.cut.off) ]
TCGA.best.cell.line.LumA              <- rs.TCGA$best.cell.line




######## MET500 ##########

load('~/Project/Cancer2CellLine/client-side/output//MET500.breast.cancer.meta.R.output//MET500.breast.cancer.meta.RData')
load('~/Project/Cancer2CellLine/server-side/RData/MET500.RData')

rs.MET500                             <- pick.out.cell.line(expr.of.samples = MET500.log2.fpkm.matrix[,MET500.breast.cancer.polyA.Basal.sample],expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
expr.cor                              <- rs.MET500$correlation.matrix[,rs.MET500$best.cell.line]
cell.line.rank                        <- apply(rs.MET500$correlation.matrix, 1, function(x) rank(x)[rs.MET500$best.cell.line] )
pure.MET.breast.cancer.Basal.sample   <- names(expr.cor)[(expr.cor > cor.cut.off)  & (cell.line.rank > rank.cut.off) ]
biopsy.site                           <- MET500.sample.meta[pure.MET.breast.cancer.Basal.sample,'biopsy.site'] %>% as.character()
pure.MET.breast.cancer.Basal.sample   <- pure.MET.breast.cancer.Basal.sample[biopsy.site == 'LIVER']
MET500.best.cell.line.Basal           <- rs.MET500$best.cell.line




rs.MET500                             <- pick.out.cell.line(expr.of.samples = MET500.log2.fpkm.matrix[,MET500.breast.cancer.polyA.LumB.sample],expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
expr.cor                              <- rs.MET500$correlation.matrix[,rs.MET500$best.cell.line]
cell.line.rank                        <- apply(rs.MET500$correlation.matrix, 1, function(x) rank(x)[rs.MET500$best.cell.line] )
pure.MET.breast.cancer.LumB.sample    <- names(expr.cor)[(expr.cor > cor.cut.off)  & (cell.line.rank > rank.cut.off) ]
biopsy.site                           <- MET500.sample.meta[pure.MET.breast.cancer.LumB.sample,'biopsy.site'] %>% as.character()
pure.MET.breast.cancer.LumB.sample    <- pure.MET.breast.cancer.LumB.sample[biopsy.site == 'LIVER']
MET500.best.cell.line.LumB            <- rs.MET500$best.cell.line


rs.MET500                             <- pick.out.cell.line(expr.of.samples = MET500.log2.fpkm.matrix[,MET500.breast.cancer.polyA.Her2.sample],expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
expr.cor                              <- rs.MET500$correlation.matrix[,rs.MET500$best.cell.line]
cell.line.rank                        <- apply(rs.MET500$correlation.matrix, 1, function(x) rank(x)[rs.MET500$best.cell.line] )
pure.MET.breast.cancer.Her2.sample    <- names(expr.cor)[(expr.cor > cor.cut.off)  & (cell.line.rank > rank.cut.off) ]
biopsy.site                           <- MET500.sample.meta[pure.MET.breast.cancer.Her2.sample,'biopsy.site'] %>% as.character()
pure.MET.breast.cancer.Her2.sample    <- pure.MET.breast.cancer.Her2.sample[biopsy.site == 'LIVER']
MET500.best.cell.line.Her2            <- rs.MET500$best.cell.line


rs.MET500                             <- pick.out.cell.line(expr.of.samples = MET500.log2.fpkm.matrix[,MET500.breast.cancer.polyA.LumA.sample],expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
expr.cor                              <- rs.MET500$correlation.matrix[,rs.MET500$best.cell.line]
cell.line.rank                        <- apply(rs.MET500$correlation.matrix, 1, function(x) rank(x)[rs.MET500$best.cell.line] )
pure.MET.breast.cancer.LumA.sample    <- names(expr.cor)[(expr.cor > cor.cut.off)  & (cell.line.rank > rank.cut.off) ]
biopsy.site                           <- MET500.sample.meta[pure.MET.breast.cancer.LumA.sample,'biopsy.site'] %>% as.character()
pure.MET.breast.cancer.LumA.sample    <- pure.MET.breast.cancer.LumA.sample[biopsy.site == 'LIVER']
MET500.best.cell.line.LumA            <- rs.MET500$best.cell.line


# attentation! SRR4238348 and SRR4305661 are duplicated samples. It is NOT our fault, maybe the author upload the wrong sequence files
pure.MET.breast.cancer.LumB.sample <- setdiff(pure.MET.breast.cancer.LumB.sample,'SRR4238348')

save(file = 'client-side/output/Select.pure.sample.breast.cancer.R.output/Select.pure.sample.breast.cancer.RData',
     list = c(
              'pure.PRI.breast.cancer.LumB.sample',  'pure.PRI.breast.cancer.LumA.sample',
              'pure.PRI.breast.cancer.Her2.sample',  'pure.PRI.breast.cancer.Basal.sample',
              'pure.MET.breast.cancer.LumB.sample',  'pure.MET.breast.cancer.LumA.sample',
              'pure.MET.breast.cancer.Her2.sample',  'pure.MET.breast.cancer.Basal.sample',
              'TCGA.best.cell.line.Basal','TCGA.best.cell.line.LumB',
              'TCGA.best.cell.line.Her2' ,'TCGA.best.cell.line.LumA',
              'MET500.best.cell.line.Basal' ,'MET500.best.cell.line.LumB',
              'MET500.best.cell.line.Her2'  ,'MET500.best.cell.line.LumA'
              
              )
)


