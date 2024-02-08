require(clusterProfiler)
library(org.Hs.eg.db)
require(dplyr)
require(pheatmap)

gene_with_protein_product <- read.delim("~/Project/BreastCancerMetaPotenial/client-side/Data/HGNC/gene_with_protein_product.txt", stringsAsFactors=FALSE)
mapping.df                <- gene_with_protein_product[,c('ensembl_gene_id','symbol')]
mapping.df                <- mapping.df[mapping.df$ensembl_gene_id != '',]
rownames(mapping.df)      <- mapping.df$ensembl_gene_id


load('client-side/output/DE.breast.cancer.R.output/DE.breast.cancer.RData')
load('client-side/output/DE.colorectal.cancer.R.output/DE.colorectal.cancer.RData')
load('client-side/output/DE.NET.pancreatic.cancer.R.output/DE.NET.pancreatic.cancer.RData')
load('client-side/output/DE.NET.si.cancer.R.output/DE.NET.si.cancer.RData')
load('client-side/output/DE.prostate.cancer.R.output/DE.prostate.cancer.RData')
load('client-side/output/DE.prostate.cancer.SRP253428.R.output/DE.prostate.cancer.SRP253428.RData')

DE.rs.list        <- list(BRCA.Basal.DE.rs, BRCA.Her2.DE.rs, BRCA.LumB.DE.rs, PRAD.DE.rs, COAD.DE.rs, NET.PAAD.DE.rs, NET.SI.DE.rs)
names(DE.rs.list) <- c(   'BRCA.Basal',     'BRCA.Her2',     'BRCA.LumB',     'PRAD',     'COAD',     'PNET',     'SINET')


perform.GO.analysis <- function(DE.rs) {
    up.gene.annotation.df  <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = intersect(DE.rs$tumor.intrinsic.DE.gene.rs$up.gene,mapping.df$ensembl_gene_id))
    up.gene.annotation.df  <- up.gene.annotation.df[complete.cases(up.gene.annotation.df),]
    GO.rs.1.up             <- enrichGO(gene=up.gene.annotation.df$SYMBOL %>% unique ,keyType='SYMBOL',OrgDb=org.Hs.eg.db,ont='BP',qvalueCutoff = 0.2,pvalueCutoff=0.001,maxGSSize = 200)@result
    GO.rs.1.up             <- GO.rs.1.up[ GO.rs.1.up$Count >= 5 , ]
    GO.rs.1.up             <- GO.rs.1.up[order(GO.rs.1.up$pvalue),]
    
    dn.gene.annotation.df  <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = intersect(DE.rs$tumor.intrinsic.DE.gene.rs$dn.gene,mapping.df$ensembl_gene_id))
    dn.gene.annotation.df  <- dn.gene.annotation.df[complete.cases(dn.gene.annotation.df),]
    GO.rs.1.dn             <- enrichGO(gene=dn.gene.annotation.df$SYMBOL %>% unique,keyType='SYMBOL',OrgDb=org.Hs.eg.db,ont='BP',qvalueCutoff = 0.2,pvalueCutoff=0.001,maxGSSize = 200)@result
    GO.rs.1.dn             <- GO.rs.1.dn[ GO.rs.1.dn$Count >= 5 , ]
    GO.rs.1.dn             <- GO.rs.1.dn[order(GO.rs.1.dn$pvalue),]
    
    list(up = GO.rs.1.up,dn = GO.rs.1.dn)
}


GO.rs.list <- foreach(DE.rs = DE.rs.list) %do% {
    perform.GO.analysis(DE.rs)
}
names(GO.rs.list) <- names(DE.rs.list)

save(file='client-side/output/DE.gene.overview.R.output/DE.gene.overview.RData',list=c('GO.rs.list'))





