source('client-side/code/util.R')
require(DESeq2)
require(org.Hs.eg.db)
require(plyr)
require(dplyr)


get.expressed.gene <- function(expr.matrix,cut.off=1){
  m.expr <- apply(expr.matrix,1,median)
  rownames(expr.matrix)[m.expr >= cut.off ]
}

############ DE between HES4-high vs HES4-low samples, TCGA samples ###########
load('server-side/RData/Breast Invasive Carcinoma.RData')
load('client-side/output/Select.pure.sample.breast.cancer.R.output/Select.pure.sample.breast.cancer.RData')

TCGA.breast.cancer.log2.tpm.matrix        <- log2.tpm.matrix
TCGA.breast.cancer.log2.read.count.matrix <- log2.read.count.matrix

HES4       <- 'ENSG00000188290'
NOTCH1     <- 'ENSG00000148400'
NOTCH3     <- 'ENSG00000074181'
HES4.expr  <- TCGA.breast.cancer.log2.tpm.matrix[HES4,pure.PRI.breast.cancer.Basal.sample] %>% sort
hist(HES4.expr,breaks=30)
quantile(HES4.expr,probs = seq(0, 1, 0.05))  %>% plot



l.sample      <- names(HES4.expr)[HES4.expr <  median(HES4.expr)]
h.sample      <- names(HES4.expr)[HES4.expr >=  median(HES4.expr)]
l.expr.matrix <- TCGA.breast.cancer.log2.tpm.matrix[,l.sample]
h.expr.matrix <- TCGA.breast.cancer.log2.tpm.matrix[,h.sample]

g1                     <- get.expressed.gene(l.expr.matrix)
g2                     <- get.expressed.gene(h.expr.matrix)
expressed.gene         <- c(g1,g2) %>% unique

tmp                  <- cbind(TCGA.breast.cancer.log2.read.count.matrix[expressed.gene,h.sample],TCGA.breast.cancer.log2.read.count.matrix[expressed.gene,l.sample])
read.count.matrix    <- 2^tmp - 1


df             <- data.frame(condition=c(rep(x='HES4.high',times=length(h.sample)), rep(x='HES4.low',times=length(l.sample))))
df$condition   <- factor(df$condition,levels = c('HES4.low','HES4.high'))

dds            <- DESeqDataSetFromMatrix(countData = round(read.count.matrix),
                                         colData = df,
                                         design = ~ condition )
dds            <- DESeq(dds)
res            <- results(dds,contrast = c('condition','HES4.high','HES4.low')) %>% as.data.frame
res            <- res[order(res$pvalue),]
res            <- res[complete.cases(res),]     
res            <- as.data.frame(res)


HES4.h.up.gene    <- rownames(res)[res$log2FoldChange >  1 & res$padj < 0.05]
HES4.h.up.gene.df <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = HES4.h.up.gene)
HES4.h.dn.gene    <- rownames(res)[res$log2FoldChange < -1 & res$padj < 0.05]
HES4.h.dn.gene.df <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = HES4.h.dn.gene)

TCGA.HES4.high.vs.low.res <- res 



############ DE between HES4-high vs HES4-low samples, SRP157974 samples ###########
load('server-side/RData/SRP157974_PrimaryTumor.RData')
require(Rtsne)
require(genefu)
TCGA.breast.cancer.polyA.sample <- sample.meta.df$sample.id[sample.meta.df$primary.disease.or.tissue == 'Breast Invasive Carcinoma']

#### Perform pam50 subtyping for TCGA polyA samples###########
pam50.gene.df               <- read.delim("~/Project/Cancer2CellLine/client-side/meta.data/pam50_entrez_gene_id_and_ensemble_gene_id.txt", stringsAsFactors=FALSE,comment.char = '#') # well, I borrow some information from the metastatic breast cancer evaluation project 
pam50.gene                  <- pam50.gene.df$ensemble.gene.id %>% as.character
colnames(pam50.gene.df)[1]  <- 'EntrezGene.ID'
colnames(pam50.gene.df)[2]  <- 'probe' #damn it, genefu package doesnot mentioning this!
pam50.gene.df$EntrezGene.ID <- as.character(pam50.gene.df$EntrezGene.ID)
pam50.gene.expr             <- cbind(TCGA.breast.cancer.log2.tpm.matrix[pam50.gene.df$probe %>% as.character,TCGA.breast.cancer.polyA.sample],
                                     SRP157974_PrimaryTumor_log2.tpm.matrix[pam50.gene.df$probe %>% as.character,]
) %>% t
annot.matrix                <- pam50.gene.df[,1:2] %>% as.matrix
rownames(annot.matrix)      <- annot.matrix[,'probe']
pam50.subtype.rs            <- intrinsic.cluster.predict(sbt.model = pam50.robust,data = pam50.gene.expr,annot = annot.matrix,mapping=annot.matrix ,do.mapping = TRUE )


dist.obj       <- as.dist(1- cor(pam50.gene.expr %>% t,method='spearman'))
set.seed(8) # I want to reproduce the tsne results, 8 is just a arbitrary numnber, it DOES NOT change the conclusion
tsne.rs        <- Rtsne(dist.obj,perplexity = 15)

draw.df <- data.frame(dim1=tsne.rs$Y[,1],
                      dim2=tsne.rs$Y[,2],
                      subtype       = pam50.subtype.rs$subtype[rownames(pam50.gene.expr)]
)
ggplot(draw.df,aes(x=dim1,y=dim2,color=subtype)) + geom_point(size=6) +  scale_shape_manual(values=c(8,15:18))

basal.sample <- names(pam50.subtype.rs$subtype)[pam50.subtype.rs$subtype == 'Basal']
basal.sample <- basal.sample[grepl(x=basal.sample,pattern = 'TCGA') == FALSE]



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

rs <- pick.out.cell.line(expr.of.samples = SRP157974_PrimaryTumor_log2.fpkm.matrix[,basal.sample], expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
r  <- apply(rs$correlation.matrix,1,function(x) rank(x)['HDQP1_BREAST'])
plot(sort(r))
basal.sample <- names(r)[r >= 1000]


cDNA.sample  <- SRP157974_PrimaryTumor_Metadata$Run[SRP157974_PrimaryTumor_Metadata$LibrarySelection =='cDNA'] %>% as.character()
basal.sample <- intersect(basal.sample,cDNA.sample)

SRP157974_PrimaryTumor_log2.tpm.matrix       <- SRP157974_PrimaryTumor_log2.tpm.matrix[,basal.sample]
SRP157974_PrimaryTumor_log2.read.count.matrix <- SRP157974_PrimaryTumor_log2.read.count.matrix[,basal.sample]


HES4.expr <- SRP157974_PrimaryTumor_log2.tpm.matrix[HES4,] %>% sort
hist(HES4.expr,breaks=40)
quantile(HES4.expr,probs = seq(0, 1, 0.05))  %>% plot


l.sample      <- names(HES4.expr)[HES4.expr  <  median(HES4.expr)]
h.sample      <- names(HES4.expr)[HES4.expr >=  median(HES4.expr)]
l.expr.matrix <- SRP157974_PrimaryTumor_log2.tpm.matrix[,l.sample]
h.expr.matrix <- SRP157974_PrimaryTumor_log2.tpm.matrix[,h.sample]

g1                     <- get.expressed.gene(l.expr.matrix)
g2                     <- get.expressed.gene(h.expr.matrix)
expressed.gene         <- c(g1,g2) %>% unique

tmp    <- cbind(SRP157974_PrimaryTumor_log2.read.count.matrix[expressed.gene,h.sample],SRP157974_PrimaryTumor_log2.read.count.matrix[expressed.gene,l.sample])
read.count.matrix    <- 2^tmp - 1
df             <- data.frame(condition=c(rep(x='HES4.high',times=length(h.sample)), rep(x='HES4.low',times=length(l.sample))))
df$condition   <- factor(df$condition,levels = c('HES4.low','HES4.high'))
dds            <- DESeqDataSetFromMatrix(countData = round(read.count.matrix),
                                         colData = df,
                                         design = ~ condition  )
dds            <- DESeq(dds)
res            <- results(dds,contrast = c('condition','HES4.high','HES4.low')) %>% as.data.frame
res            <- res[order(res$pvalue),]
res            <- res[complete.cases(res),]     
res            <- as.data.frame(res)
SRP157974.HES4.high.vs.low.res <- res 
SRP157974.HES4.high.vs.low.res[NOTCH3,]

save(file='client-side/output/NOTCH3-HES4.R.output/NOTCH3-HES4.RData',list =c ('SRP157974.HES4.high.vs.low.res','TCGA.HES4.high.vs.low.res'))



