require(ggplot2)
require(CePa)
require(dplyr)
require(ComplexHeatmap)
require(RColorBrewer)
require(circlize)
require(segmented)
source('client-side/code/Manuscript/ggplot.style.R')

load('client-side/output/DE.breast.cancer.R.output/DE.breast.cancer.RData')
load('client-side/output/DE.prostate.cancer.R.output/DE.prostate.cancer.RData')
load('client-side/output/DE.colorectal.cancer.R.output/DE.colorectal.cancer.RData')
load('client-side/output/DE.NET.pancreatic.cancer.R.output/DE.NET.pancreatic.cancer.RData')
load('client-side/output/DE.NET.si.cancer.R.output/DE.NET.si.cancer.RData')

DE.rs.list        <- list(BRCA.LumB.DE.rs, BRCA.Basal.DE.rs, BRCA.Her2.DE.rs, PRAD.DE.rs,COAD.DE.rs,NET.PAAD.DE.rs,NET.SI.DE.rs)
names(DE.rs.list) <- c('BRCA.LumB','BRCA.Basal', 'BRCA.Her2','PRAD', 'COAD', 'PNET', 'SINET')

for(cancer.type in names(DE.rs.list)) {
  DE.rs         <- DE.rs.list[[cancer.type]]
  DE.rs.R.vs.P  <- DE.rs$deseq2.R.vs.P.res
  DE.rs.M.vs.P  <- DE.rs$deseq2.M.vs.P.res
  c.gene        <- intersect(rownames(DE.rs.R.vs.P),rownames(DE.rs.M.vs.P))
  x             <- DE.rs.R.vs.P[c.gene,'log2FoldChange']
  y             <- DE.rs.M.vs.P[c.gene,'log2FoldChange']
  lin.mod       <- lm(y~x)
  psi           <- -1
  while(psi < 0){ # well, for BRCA.Her2, this is necessary! psi = 6.595793
    lin.mod       <- lm(y~x)
    segmented.mod <- segmented(lin.mod, seg.Z = ~x, psi=2)    
    tmp           <- summary(segmented.mod)
    psi           <- tmp$psi[1,'Est.']
  }
  df            <- data.frame(x=DE.rs.R.vs.P[c.gene,'log2FoldChange'],y=DE.rs.M.vs.P[c.gene,'log2FoldChange'],fitted = fitted(segmented.mod),residual = segmented.mod$residuals)
  rownames(df)  <- c.gene
  
  
  #########################################################
  ## (a) : log2FC (MET.vs.PRI)') log2FC (LIVER.vs.PRI)')
  ######################################################## 

  
ggplot(df[df$x > psi,],aes(x=x,y=y)) + geom_point(size=2.5) + ggplot.style + 
    ylab('log2FC (MET.vs.PRI)') + xlab('log2FC (LIVER.vs.PRI)') + 
    xlim(c(-15,22))   + 
    geom_vline(xintercept =psi,linetype=2,size=2) + 
    geom_line(aes(x=x,y=fitted),col='red',lwd=3.5) +  geom_point(data = df[df$x < psi,], aes(x=x,y=y),color='grey',size=2.5)
  file.name <- sprintf('~/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section5/%s.log2FC.MP.and.normal.grey.pdf',cancer.type)
  ggsave(filename = file.name ,width = 20,height=20)
  
  
 h <- df[df$x > psi,]
  h$r <- scale(h$residual)
  ggplot(h,aes(x=x,y=r)) + geom_point(size=2.5) + ggplot.style + 
    ylab('log2FC (MET.vs.PRI)') + xlab('log2FC (LIVER.vs.PRI)') + 
    xlim(c(0,22))  +  
    geom_point(data = h[h$r > 3,],aes(x=x,y=r),color='red',size = 8.5 ) +
    geom_hline(yintercept =3,linetype=2,size=2) 
  file.name <- sprintf('~/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section5/%s.log2FC.MP.and.normal.residual.pdf',cancer.type)
  ggsave(filename = file.name ,width = 20,height=20)
  
  print(sprintf('cancer.type = %s, psi = %f',cancer.type,psi))
  

#############################################
## (b) : boxplot of CPN1's max expression 
#############################################
load('client-side/output/COAD.scRNAseq.analysis.R.output/COAD.scRNAseq.analysis.RData')
max.expr <- apply(COAD.ectopic.liver.gene.expr.matrix,1,max)
df <- data.frame(max.expr = max.expr) 
ggplot(df,aes(x='COAD',y=(max.expr))) + geom_boxplot(outlier.shape=NA,lwd=2)  + ggplot.style + geom_jitter(aes(x='COAD',y=max.expr),size=8.5) 

file.name <- sprintf('~/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section5/max.value.box.plot.pdf')
ggsave(filename = file.name ,width = 20,height=20)


###################################################
## (d) : PNET SINET log2FC.MP.and.normal.residual
###################################################

if(cancer.type  %in% c('PNET','SINET')){
  file.name <- sprintf('~/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section5/%s.log2FC.MP.and.normal.residual.zoomin.pdf',cancer.type)
  ggplot(h[h$r > 3,],aes(x=x,y=r)) + geom_point(size = 6,color='red') + ggplot.style
  ggsave(filename = file.name ,width = 20,height=20)         
}

}

