require(plyr)
require(dplyr)
require(genefu)
require(Rtsne)
require(ggplot2)
require(dplyr)
source('client-side/code/util.R')


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




load('server-side/RData/COLORECTAL_SRP029880.RData')
met.sample <- COLORECTAL_SRP029880_Metadata$Run[COLORECTAL_SRP029880_Metadata$tissue == 'metastatic colorectal cancer to the liver'] %>% as.character()
pri.sample <- COLORECTAL_SRP029880_Metadata$Run[COLORECTAL_SRP029880_Metadata$tissue == 'primary colorectal cancer'] %>% as.character()






cor.cut.off <- 0.3
rank.cut.off <- 1010

TCGA.colorectal.cancer.log2.fpkm.matrix  <- COLORECTAL_SRP029880_log2.fpkm.matrix[,pri.sample]
rs.TCGA                                  <- pick.out.cell.line(expr.of.samples = TCGA.colorectal.cancer.log2.fpkm.matrix,expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
expr.cor                                 <- rs.TCGA$correlation.matrix[,rs.TCGA$best.cell.line]
cell.line.rank                           <- apply(rs.TCGA$correlation.matrix, 1, function(x) rank(x)[rs.TCGA$best.cell.line] )
pure.PRI.colorectal.cancer.sample        <- names(expr.cor)[(expr.cor > cor.cut.off)  & (cell.line.rank > rank.cut.off) ]




MET500.colorectal.cancer.log2.fpkm.matrix  <- COLORECTAL_SRP029880_log2.fpkm.matrix[,met.sample]
rs.MET500                                  <- pick.out.cell.line(expr.of.samples = MET500.colorectal.cancer.log2.fpkm.matrix,expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
expr.cor                                   <- rs.MET500$correlation.matrix[,rs.MET500$best.cell.line]
cell.line.rank                             <- apply(rs.MET500$correlation.matrix, 1, function(x) rank(x)[rs.MET500$best.cell.line] )
pure.MET.colorectal.cancer.sample          <- names(expr.cor)[(expr.cor > cor.cut.off)  & (cell.line.rank > rank.cut.off) ]

MET500.best.cell.line <- rs.MET500$best.cell.line
TCGA.best.cell.line   <- rs.TCGA$best.cell.line


save(file = 'client-side/output/Select.pure.sample.colorectal.cancer.R.output/Select.pure.sample.colorectal.cancer.RData',list=c('pure.PRI.colorectal.cancer.sample','pure.MET.colorectal.cancer.sample','MET500.best.cell.line','TCGA.best.cell.line'))




