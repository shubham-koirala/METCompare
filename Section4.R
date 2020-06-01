require(ggplot2)
source('client-side/code/Manuscript/ggplot.style.R')

#################################################################
#### Fig 4a Compare NOTCH3 CNV between primary and metestatic cancer
#################################################################
load('~/Project/Cancer2CellLine/client-side/output//MET500.breast.cancer.meta.R.output//MET500.breast.cancer.meta.RData')
load('~/Project/Cancer2CellLine/server-side/RData/MET500.RData')
load('client-side/output/TCGA.breast.cancer.meta.R.output/TCGA.breast.cancer.meta.RData')
load('client-side/output/estimate.hepatocytes.abundance.R.output/estimate.hepatocytes.abundance.RData')
load('client-side/output/organize.TCGA.and.MET500.breast.cancer.cnv.data.R.output/organize.TCGA.and.MET500.breast.cancer.cnv.data.RData')

MET500.liver.sample  <- rownames(MET500.sample.meta)[MET500.sample.meta$biopsy.site == 'LIVER'] 
MET500.sample        <- MET500.breast.cancer.polyA.Basal.sample
MET500.sample        <- intersect(MET500.sample,MET500.liver.sample)
MET500.sample        <- MET500.sample[hepatocyte.abundance.vec[MET500.sample] < 0.2]
TCGA.sample          <- pure.TCGA.breast.cancer.polyA.Basal.sample

MET500.subject.id    <- MET500.sample.meta[MET500.sample %>% as.character(),'MET500.id']
TCGA.subject.id      <- intersect(TCGA.sample,colnames(TCGA.breast.cancer.cnv.matrix))

df1                  <- data.frame(cnv= MET500.breast.cancer.cnv.matrix['NOTCH3',MET500.subject.id], site= rep('MET500',length(MET500.subject.id)))
df2                  <- data.frame(cnv= TCGA.breast.cancer.cnv.matrix  ['NOTCH3',TCGA.subject.id],   site= rep('TCGA',length(TCGA.subject.id)))
draw.df              <- rbind(df1,df2)
draw.df$site         <- factor(draw.df$site,levels = c('TCGA','MET500'))

ggplot(draw.df) + geom_boxplot(aes(x=site,y=cnv),outlier.shape=NA,lwd=3) + ggplot.style + geom_jitter(aes(x=site,y=cnv),size=5.5) + xlab('')
wilcox.test(MET500.breast.cancer.cnv.matrix['NOTCH3',MET500.subject.id],TCGA.breast.cancer.cnv.matrix['NOTCH3',TCGA.subject.id])
ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section4/NOTCH3.cnv.basal.pdf',width = 20,height=20)


#################################################################
#### Fig 4c  DE analysis between HES4.high vs HES4.low
#################################################################
load('client-side/output/NOTCH3-HES4.R.output/NOTCH3-HES4.RData')
NOTCH3     <- 'ENSG00000074181'
NOTCH3.df <- TCGA.HES4.high.vs.low.res[rownames(TCGA.HES4.high.vs.low.res) == NOTCH3,]
ggplot() + geom_point(data=TCGA.HES4.high.vs.low.res,aes(x=log2FoldChange,y= -1 * log10(padj)),size=2.5) + ggplot.style + ylim(0,15) + geom_point(data=NOTCH3.df,aes(x=log2FoldChange,y= -1 * log10(padj)),size=7.5,color='red')
ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section4/TCGA.HES4.high.vs.low.de.pdf',width = 20,height=20)



#################################################################
#### Fig 4d  KM plot of HES4 in basal-like breast cancer
#################################################################

load('client-side/output/DE.gene.clinics.R.output/DE.gene.clinics.RData')
load('client-side/output//TCGA.breast.cancer.meta.R.output//TCGA.breast.cancer.meta.RData')
load('server-side/RData//Breast Invasive Carcinoma.RData')
load('~/Project/InSilicoCRISPR/client-side/output/organize.TCGA.clinical.data.R.output/organize.TCGA.clinical.data.RData')
#load('client-side/output/tumor.purity.based.on.cell.line.R.output//tumor.purity.based.on.cell.line.RData')
require(survival)
require(survminer)


TCGA.breast.cancer.log2.fpkm.matrix       <- log2.fpkm.matrix
clinical.data                             <- TCGA.clinical.data.list$BRCA


draw.KM.plot <- function(g) {
  expr           <- TCGA.breast.cancer.log2.fpkm.matrix[g,subtype.sample]
  patient.name   <- gsub(expr %>% names,pattern = '\\-01',replacement = '')
  df <- data.frame(expr=expr,
                   status=clinical.data[patient.name,'vital_status'], 
                   d1=clinical.data[patient.name,'days_to_death'] ,
                   d2=clinical.data[patient.name,'days_to_last_followup'],
                   age= clinical.data[patient.name,'years_to_birth'],
                   pathologic_stage= clinical.data[patient.name,'pathologic_stage']
  )
  df[is.na(df$d1),'d1'] <- 0
  df[is.na(df$d2),'d2'] <- 0
  df$time <- (df$d1 + df$d2) / 30
  df      <- df[complete.cases(df),]
  q                          <- quantile(df$expr)
  df$group                   <- 'median'
  df[df$expr > q[3],'group'] <- 'high'
  df[df$expr <= q[3],'group']<- 'low'
  df                         <- df[df$group!= 'median',]
  df$group                   <- factor(df$group,levels = c('low','high'))
  fit                        <- survfit(Surv(time, status) ~ group , data=df)
  km.style <- theme_bw(base_size = 55) + theme(axis.title = element_text( size=30, face="bold"),
                                               axis.text  = element_text( size=70, face="bold"),
                                               plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
                                               axis.line.x = element_line(colour = "black",size = 3),
                                               axis.line.y = element_line(colour = "black",size = 3),
                                               plot.title = element_text(size =40, face = "bold")
  ) 
  ggsurvplot(fit=fit,data=df,legend='none',ggtheme = km.style,xlab='',ylab='',palette = c("blue", "red"),xlim=c(0,350),title=g,size=5,break.x.by=100) 
}


HES4  <- 'ENSG00000188290'
subtype.sample           <- intersect(pure.TCGA.breast.cancer.polyA.Basal.sample,colnames(TCGA.breast.cancer.log2.fpkm.matrix))
HES4.KM.plot <- draw.KM.plot(HES4)
basal.survival.rs.up[HES4,]
pdf(file ='~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section4/HES4.KM.plot.pdf',width=20,height=15)
print(HES4.KM.plot[[1]])
dev.off()



#################################################################
#### Fig S6a  cnv value of NOTCH genes in MET500
#################################################################
load('client-side/output/coompare.CNV.between.TCGA.and.MET500.Basal.R.output/coompare.CNV.between.TCGA.and.MET500.Basal.RData')
NOTCH.gene.MET500.vs.TCGA.cnv.comparision.df <- NOTCH.gene.MET500.vs.TCGA.cnv.comparision.df[order(NOTCH.gene.MET500.vs.TCGA.cnv.comparision.df$fdr),]
sig.NOTCH.gene <- rownames(NOTCH.gene.MET500.vs.TCGA.cnv.comparision.df)[NOTCH.gene.MET500.vs.TCGA.cnv.comparision.df$fdr < 0.05]

draw.df <- foreach(g= sig.NOTCH.gene,.combine = 'rbind') %do% {
  data.frame(cnv= MET500.breast.cancer.cnv.matrix[g,MET500.subject.id], gene=g)
  
}
draw.df$gene <- factor(draw.df$gene,levels=sig.NOTCH.gene %>% as.character())
ggplot(draw.df)  + geom_point(aes(x=gene,y=cnv),size=5.5) + coord_flip() + ggplot.style + geom_abline(intercept = 0,slope = 0,linetype=2,size=2)
ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section4/significant.NOTCH.gene.cnv.in.MET500.pdf',width = 20,height=20)

#################################################################
#### Fig S6b  CNV-expression correlation
#################################################################
NOTCH3    <- 'ENSG00000074181'
df.TCGA   <- data.frame (expr= TCGA.breast.cancer.log2.fpkm.matrix[NOTCH3,TCGA.subject.id],cnv=TCGA.breast.cancer.cnv.matrix  ['NOTCH3',TCGA.subject.id])
df.MET500 <- data.frame (expr= MET500.log2.fpkm.matrix[NOTCH3,MET500.sample],              cnv=MET500.breast.cancer.cnv.matrix['NOTCH3',MET500.subject.id])
ggplot(rbind(df.MET500,df.TCGA)) + geom_point(aes(x=cnv,y=expr),size=5.5) + ggplot.style
ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section4/NOTCH3.cnv.expr.correlation.pdf',width = 20,height=20)



#################################################################
#### Fig S6c  boxplot of NOTCH3 cnv in Her2 and LuminalB subtype
#################################################################
MET500.liver.sample  <- rownames(MET500.sample.meta)[MET500.sample.meta$biopsy.site == 'LIVER'] 
MET500.sample        <- MET500.breast.cancer.polyA.LumB.sample
MET500.sample        <- intersect(MET500.sample,MET500.liver.sample)
MET500.sample        <- MET500.sample[hepatocyte.abundance.vec[MET500.sample] < 0.2]
TCGA.sample          <- pure.TCGA.breast.cancer.polyA.LumB.sample
MET500.subject.id    <- MET500.sample.meta[MET500.sample %>% as.character(),'MET500.id']
TCGA.subject.id      <- intersect(TCGA.sample,colnames(TCGA.breast.cancer.cnv.matrix))
df1                  <- data.frame(cnv= MET500.breast.cancer.cnv.matrix['NOTCH3',MET500.subject.id], site= rep('MET500',length(MET500.subject.id)))
df2                  <- data.frame(cnv= TCGA.breast.cancer.cnv.matrix  ['NOTCH3',TCGA.subject.id],   site= rep('TCGA',length(TCGA.subject.id)))
draw.df              <- rbind(df1,df2)
draw.df$site         <- factor(draw.df$site,levels = c('TCGA','MET500'))
ggplot(draw.df) + geom_boxplot(aes(x=site,y=cnv),outlier.shape=NA,lwd=2) + ggplot.style + geom_jitter(aes(x=site,y=cnv),size=5.5) + xlab('')
wilcox.test(MET500.breast.cancer.cnv.matrix['NOTCH3',MET500.subject.id],TCGA.breast.cancer.cnv.matrix['NOTCH3',TCGA.subject.id])
ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section4/NOTCH3.cnv.lumb.pdf',width = 20,height=20)


MET500.liver.sample  <- rownames(MET500.sample.meta)[MET500.sample.meta$biopsy.site == 'LIVER'] 
MET500.sample        <- MET500.breast.cancer.polyA.Her2.sample
MET500.sample        <- intersect(MET500.sample,MET500.liver.sample)
MET500.sample        <- MET500.sample[hepatocyte.abundance.vec[MET500.sample] < 0.2]
TCGA.sample          <- pure.TCGA.breast.cancer.polyA.Her2.sample
MET500.subject.id    <- MET500.sample.meta[MET500.sample %>% as.character(),'MET500.id']
TCGA.subject.id      <- intersect(TCGA.sample,colnames(TCGA.breast.cancer.cnv.matrix))
df1                  <- data.frame(cnv= MET500.breast.cancer.cnv.matrix['NOTCH3',MET500.subject.id], site= rep('MET500',length(MET500.subject.id)))
df2                  <- data.frame(cnv= TCGA.breast.cancer.cnv.matrix  ['NOTCH3',TCGA.subject.id],   site= rep('TCGA',length(TCGA.subject.id)))
draw.df              <- rbind(df1,df2)
draw.df$site         <- factor(draw.df$site,levels = c('TCGA','MET500'))
ggplot(draw.df) + geom_boxplot(aes(x=site,y=cnv),outlier.shape=NA,lwd=2) + ggplot.style + geom_jitter(aes(x=site,y=cnv),size=5.5) + xlab('')
wilcox.test(MET500.breast.cancer.cnv.matrix['NOTCH3',MET500.subject.id],TCGA.breast.cancer.cnv.matrix['NOTCH3',TCGA.subject.id])
ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section4/NOTCH3.cnv.her2.pdf',width = 20,height=20)




########################################################################
### Fig S6d valcano plot of SRP157974
############################################################################
NOTCH3.df <- SRP157974.HES4.high.vs.low.res[rownames(SRP157974.HES4.high.vs.low.res) == NOTCH3,]
ggplot() + geom_point(data=SRP157974.HES4.high.vs.low.res,aes(x=log2FoldChange,y= -1 * log10(padj)),size=2.5) + ggplot.style + ylim(0,15) + geom_point(data=NOTCH3.df,aes(x=log2FoldChange,y= -1 * log10(padj)),size=7.5,color='red')
ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section4/SRP157974.HES4.high.vs.low.de.pdf',width = 20,height=20)






