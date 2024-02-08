require(dplyr)

load('client-side/Data/GSE97693/COLORECTAL_GSE97693_single.RData')
flag <- GSE97693_FPKM$X %in% c('1-Mar','2-Mar')
GSE97693_FPKM <- GSE97693_FPKM[!flag,]
rownames(GSE97693_FPKM) <- GSE97693_FPKM$X
GSE97693_FPKM$X <- NULL
GSE97693_FPKM.matrix <- apply(GSE97693_FPKM,2,as.numeric)
rownames(GSE97693_FPKM.matrix) <- rownames(GSE97693_FPKM)
flag   <- GSE97693_FPKM_METADATA$Sample_class %in% c(' Post-treatment Liver Metastasis','Liver Metastasis')
gsm.id <- GSE97693_FPKM_METADATA$GSM[flag] %>% as.character()


load('client-side/output/DE.colorectal.cancer.R.output/DE.colorectal.cancer.RData')
gene_with_protein_product <- read.delim("~/Project/BreastCancerMetaPotenial/client-side/Data/HGNC/gene_with_protein_product.txt", stringsAsFactors=FALSE)
mapping.df                <- gene_with_protein_product[,c('ensembl_gene_id','symbol')]
mapping.df                <- mapping.df[mapping.df$ensembl_gene_id != '',]
rownames(mapping.df)      <- mapping.df$ensembl_gene_id

g.vec <- mapping.df[COAD.DE.rs$tumor.intrinsic.DE.gene.rs$up.outlier.gene,'symbol']


GSE97693_FPKM.matrix['CD5L',gsm.id] %>% sort %>% boxplot
GSE97693_FPKM.matrix['SLC13A5',gsm.id] %>% sort %>% boxplot
GSE97693_FPKM.matrix['HABP2',gsm.id] %>% sort %>% boxplot
GSE97693_FPKM.matrix['F11',gsm.id] %>% sort %>% boxplot
GSE97693_FPKM.matrix['IGF2',gsm.id] %>% sort %>% boxplot
GSE97693_FPKM.matrix['MARCO',gsm.id] %>% sort %>% boxplot
GSE97693_FPKM.matrix['CPN1',gsm.id] %>% sort %>% boxplot
GSE97693_FPKM.matrix['INS-IGF2',gsm.id] %>% sort %>% boxplot


COAD.ectopic.liver.gene.expr.matrix <- GSE97693_FPKM.matrix[c('CD5L','SLC13A5','HABP2','F11','IGF2','MARCO','CPN1','INS-IGF2'),gsm.id]


GSE97693_FPKM.matrix['CPN1',gsm.id] %>% sort(decreasing = TRUE) %>% head
GSM2697027.expr.vec <- GSE97693_FPKM.matrix[,'GSM2697027'] # GSM2697027 has highest CPN1 expression



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

g            <- intersect(rownames(CCLE.log2.rpkm.matrix),mapping.df$ensembl_gene_id)
CC           <- CCLE.log2.rpkm.matrix[g,]
rownames(CC) <- mapping.df[g,'symbol']

m         <- mapping.df[CCLE.rna.seq.marker.gene.1000,'symbol']
m         <- m[is.na(m) == FALSE]
m         <- intersect(names(GSM2697027.expr.vec),m)
GSM2697027.cor.value <- cor(CC[m,],GSM2697027.expr.vec[m],method='spearman')

save(file='client-side/output/COAD.scRNAseq.analysis.R.output/COAD.scRNAseq.analysis.RData',list=c('GSM2697027.cor.value','COAD.ectopic.liver.gene.expr.matrix'))





