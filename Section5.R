require(ggplot2)
require(dplyr)
source('client-side/code/Manuscript/ggplot.style.R')

load('client-side/output/pharmacological.interaction.R.output/pharmacological.interaction.RData')


####################################################################################
### Fig 5a, the pipeline
#####################################################################################


####################################################################################
### Fig 5b, scatterplot of RNAi and CRISPR residual score, BASAL cell line 
#####################################################################################

c.gene            <- intersect(rownames(BASAL.rs$CRISPR.result.df), rownames(BASAL.rs$RNAi.result.df))
draw.df           <- data.frame(x=BASAL.rs$CRISPR.result.df[c.gene,'sr'],y= BASAL.rs$RNAi.result.df[c.gene,'sr'])
rownames(draw.df) <- c.gene
ggplot(draw.df,aes(x=x,y=y)) + geom_point(size=8.5) + ggplot.style + ylab('RNAi') + xlab('CRISPR') + xlim(-10,10) + ylim(-10,10) 
ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section5/BASAL.screen.rs.pdf',width = 20,height=20)



####################################################################################
### Fig 5c, scatterplot of RNAi and CRISPR residual score, LYMPHOID 
#####################################################################################

c.gene            <- intersect(rownames(LYMPHOID.rs$CRISPR.result.df), rownames(LYMPHOID.rs$RNAi.result.df))
draw.df           <- data.frame(x=LYMPHOID.rs$CRISPR.result.df[c.gene,'sr'],y= LYMPHOID.rs$RNAi.result.df[c.gene,'sr'])
rownames(draw.df) <- c.gene
ggplot(draw.df,aes(x=x,y=y)) + geom_point(size=8.5) + ggplot.style + ylab('RNAi') + xlab('CRISPR') + xlim(-12,10) + ylim(-12,10) 
ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section5/LYMPHOID.screen.rs.pdf',width = 20,height=20)



####################################################################################
### Fig 5d, scatterplot of GDSC score 
#####################################################################################

ggplot(GDSC.IC50.result.df,aes(x=pan.cell.line.score,y=amp.cell.line.score)) + geom_point(size=8.5) + ggplot.style + ylab('amp.cell.line') + xlab('pan.cell.line') + xlim(-9,9) + ylim(-9,9) +  geom_abline(intercept = 0, slope = 1,linetype="dashed",size=2,colour='red') 
ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section5/GDSC.screen.rs.pdf',width = 20,height=20)




##############################################################################################
### Fig S6a,b: scattler plot of median screen scroe in the two groups, BASAL breast cancer cell lines  
###############################################################################################

ggplot(BASAL.rs$CRISPR.result.df,aes(x=pan.cell.line.score,y=amp.cell.line.score))  + geom_point(size=3.5) + ggplot.style + ylab('amp.cell.line') + xlab('pan.cell.line') + xlim(-6,6) + ylim(-6,6) +  geom_abline(intercept = 0, slope = 1,linetype="dashed",size=2,colour='red')
ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section5/BASAL.CRISPR.pdf',width = 20,height=20)

ggplot(BASAL.rs$RNAi.result.df,aes(x=pan.cell.line.score,y=amp.cell.line.score))  + geom_point(size=3.5) + ggplot.style + ylab('amp.cell.line') + xlab('pan.cell.line') + xlim(-2,1.5) + ylim(-2,1.5) +  geom_abline(intercept = 0, slope = 1,linetype="dashed",size=2,colour='red')
ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section5/BASAL.RNAi.pdf',width = 20,height=20)



##############################################################################################
### Fig S6c,d: scattler plot of median screen scroe in the two groups, LYMPHOID cancer cell lines  
###############################################################################################

ggplot(LYMPHOID.rs$CRISPR.result.df,aes(x=pan.cell.line.score,y=amp.cell.line.score))  + geom_point(size=3.5) + ggplot.style + ylab('amp.cell.line') + xlab('pan.cell.line') + xlim(-6,6) + ylim(-6,6) +  geom_abline(intercept = 0, slope = 1,linetype="dashed",size=2,colour='red')
ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section5/LYMPHOID.CRISPR.pdf',width = 20,height=20)


ggplot(LYMPHOID.rs$RNAi.result.df,aes(x=pan.cell.line.score,y=amp.cell.line.score))  + geom_point(size=3.5) + ggplot.style + ylab('amp.cell.line') + xlab('pan.cell.line') + xlim(-2,1.5) + ylim(-2,1.5) +  geom_abline(intercept = 0, slope = 1,linetype="dashed",size=2,colour='red')
ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section5/LYMPHOID.RNAi.pdf',width = 20,height=20)





