require(ggplot2)
require(org.Hs.eg.db)
require(AnnotationDbi)
source('client-side/code/Manuscript/ggplot.style.R')
load('server-side/RData/Breast Invasive Carcinoma.RData')
load('client-side/output/TCGA.breast.cancer.meta.R.output/TCGA.breast.cancer.meta.RData')
load('~/Project/Cancer2CellLine/client-side/output//MET500.breast.cancer.meta.R.output//MET500.breast.cancer.meta.RData')
load('~/Project/Cancer2CellLine/server-side/RData/MET500.RData')
load('client-side/output/estimate.hepatocytes.abundance.R.output/estimate.hepatocytes.abundance.RData')


TCGA.breast.cancer.log2.fpkm.matrix       <- log2.fpkm.matrix
MET500.liver.sample                       <- rownames(MET500.sample.meta)[MET500.sample.meta$biopsy.site == 'LIVER'] 
PRRX1                                     <- 'ENSG00000116132'

################################################################################   
#----------------------------- Figure 3 ---------------------------------------#
################################################################################  



#############################################
### Fig 3b, PRRX1 expression 
#############################################
MET500.sample  <- MET500.breast.cancer.polyA.Basal.sample
MET500.sample  <- intersect(MET500.sample,MET500.liver.sample)
MET500.sample  <- MET500.sample[hepatocyte.abundance.vec[MET500.sample] < 0.2]
TCGA.sample    <- pure.TCGA.breast.cancer.polyA.Basal.sample
df1            <- data.frame(expr=MET500.log2.fpkm.matrix[PRRX1,MET500.sample],condition='MET500')
df2            <- data.frame(expr=TCGA.breast.cancer.log2.fpkm.matrix[PRRX1,TCGA.sample],condition='TCGA')
df             <- rbind(df1,df2)
df$condition   <- factor(df$condition,levels = c('TCGA','MET500'))
wilcox.test(df1$expr,df2$expr)
df.basal         <- df
df.basal$subtype <- 'Basal-like'


MET500.sample  <- MET500.breast.cancer.polyA.Her2.sample
MET500.sample  <- intersect(MET500.sample,MET500.liver.sample)
MET500.sample  <- MET500.sample[hepatocyte.abundance.vec[MET500.sample] < 0.2]
TCGA.sample    <- pure.TCGA.breast.cancer.polyA.Her2.sample
df1            <- data.frame(expr=MET500.log2.fpkm.matrix[PRRX1,MET500.sample],condition='MET500')
df2            <- data.frame(expr=TCGA.breast.cancer.log2.fpkm.matrix[PRRX1,TCGA.sample],condition='TCGA')
df             <- rbind(df1,df2)
df$condition   <- factor(df$condition,levels = c('TCGA','MET500'))
wilcox.test(df1$expr,df2$expr)
df.her2         <- df
df.her2$subtype <- 'Her2-enriched'



MET500.sample  <- MET500.breast.cancer.polyA.LumB.sample
MET500.sample  <- intersect(MET500.sample,MET500.liver.sample)
MET500.sample  <- MET500.sample[hepatocyte.abundance.vec[MET500.sample] < 0.2]
TCGA.sample    <- pure.TCGA.breast.cancer.polyA.LumB.sample
df1            <- data.frame(expr=MET500.log2.fpkm.matrix[PRRX1,MET500.sample],condition='MET500')
df2            <- data.frame(expr=TCGA.breast.cancer.log2.fpkm.matrix[PRRX1,TCGA.sample],condition='TCGA')
df             <- rbind(df1,df2)
df$condition   <- factor(df$condition,levels = c('TCGA','MET500'))
wilcox.test(df1$expr,df2$expr)
df.lumb         <- df
df.lumb$subtype <- 'LuminalB'

draw.df         <- rbind(df.basal,df.her2,df.lumb)
draw.df$subtype <- factor(draw.df$subtype,levels = c('Basal-like','LuminalB','Her2-enriched'))
ggplot(draw.df,aes(x=subtype,y=expr,fill=condition)) +
  geom_boxplot(lwd=1.5,outlier.shape = NA) + 
  geom_point(position=position_jitterdodge(),size=3.0) + 
  scale_fill_manual(values = c('MET500'='red','TCGA'='grey')) + 
  theme_bw(base_size = 55) + theme(axis.title = element_text( size=25, face="bold"),
                                   axis.text.y  = element_text( size=55, face="bold"),
                                   plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
                                   axis.line.x = element_line(colour = "black",size = 3),
                                   axis.line.y = element_line(colour = "black",size = 3),
                                   axis.text.x = element_text(angle = 45, hjust = 1,size=10, face="bold"),
                                   legend.position= 'none')            + xlab('') 
ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section3/PRRX1.boxplot.pdf',width = 20,height=20)








#################################################
# Fig 3a, co-expr between the ten genes and PRRX1 
#################################################
load('client-side/output/PRRX1.R.output/PRRX1.RData')

draw.df <- rbind(Basal.cor.df,LumB.cor.df,Her2.cor.df)
draw.df$subtype <- factor(draw.df$subtype,levels = c('Basal-like','LuminalB','Her2-enriched'))
ggplot(draw.df,aes(x=subtype,y=cor.value,fill=gene.type)) + geom_violin(lwd=3)  + 
  scale_fill_manual(values = c('other.gene'='grey','ECM.gene'='blue')) + 
  theme_bw(base_size = 55) + theme(axis.title = element_text( size=25, face="bold"),
                                   axis.text.y  = element_text( size=55, face="bold"),
                                   plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
                                   axis.line.x = element_line(colour = "black",size = 3),
                                   axis.line.y = element_line(colour = "black",size = 3),
                                   axis.text.x = element_text(angle = 45, hjust = 1,size=10, face="bold"),
                                   legend.position= 'none')            + xlab('')  + ylim (-1.2,1.2)

ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section3/co.expr.violin.plot.pdf',width = 20,height=20)







################################################################################   
#----------------------------- Figure S5 ---------------------------------------#
################################################################################  


####################################################################    
### Fig 54a: boxplot of PRRX1 in dataset GSE58078 
####################################################################    
load('server-side/RData/BRACA_SRP043470.RData')
pri.expr         <- BRACA_SRP043470_log2.fpkm.matrix[PRRX1,c('SRR1427482','SRR1427483','SRR1427484')]
liver.met.expr   <- BRACA_SRP043470_log2.fpkm.matrix[PRRX1,c('SRR1427487','SRR1427488','SRR1427489')]

df1            <- data.frame(expr=pri.expr,condition='PRI')
df2            <- data.frame(expr=liver.met.expr,condition='MET')
df             <- rbind(df1,df2)
df$condition   <- factor(df$condition,levels = c('PRI','MET'))
ggplot(df) + geom_boxplot(aes(x=condition,y=expr),outlier.shape=NA,lwd=2)  + geom_jitter(aes(x=condition,y=expr),size=5.5) + xlab('') + ggplot.style
ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section3/PRRX1.boxplot.SRP043470.pdf',width = 20,height=20)
wilcox.test(df1$expr,df2$expr,paired=TRUE)
t.test(df1$expr,df2$expr,paired=TRUE)

# p-value two high, maybe due to sample size
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4289456/


###############################
### Fig S4b: boxplot of PRRX1 in brain metastsis and TCGA  
###############################
MET500.brain.sample <- rownames(MET500.sample.meta)[MET500.sample.meta$biopsy.site == 'BRAIN'] 
MET500.sample       <- c(MET500.breast.cancer.polyA.Basal.sample)
MET500.sample       <- intersect(MET500.sample,MET500.brain.sample)

TCGA.sample    <- pure.TCGA.breast.cancer.polyA.Basal.sample
df1            <- data.frame(expr=TCGA.breast.cancer.log2.fpkm.matrix[PRRX1,TCGA.sample],condition='TCGA')
df2            <- data.frame(expr=MET500.log2.fpkm.matrix[PRRX1,MET500.sample],          condition='MET500')
draw.df <- rbind(df1,df2)
ggplot(draw.df) + geom_boxplot(aes(x=condition,y=expr),outlier.shape=NA,lwd=1.2) + ggplot.style + geom_jitter(aes(x=condition,y=expr),size=4.5) + xlab('') + ylim(0,9)
ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section3/PRRX1.brain.metastasis.pdf',width = 20,height=20)
wilcox.test(df1$expr,df2$expr)









