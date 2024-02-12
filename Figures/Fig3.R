##################################################################################################################################
#--------------------------- Figure 3 ---------------------------------------#
##################################################################################################################################

#################################################################
# (a): Table to show DE gene clusters
#################################################################


#################################################################
# (b): screen shot of UCSC genome browser to show the genomic amplification on chr19 (BRCA.Basal)
#################################################################


#################################################################
# (c): Survival analysis 19p13.12-amp vs normal
#################################################################
require(ggplot2)
source('client-side/code/Manuscript/ggplot.style.R')

library(cgdsr)
library(survival)
library(survminer)

mycgds  <-  CGDS('http://www.cbioportal.org/',verbose=TRUE)
symbols <- c("CASP14", "ILVBL", "NOTCH3", "BRD4", "AKAP8L", "WIZ") 


mycancerstudy  <-  "brca_metabric" 
mycaselist     <- "brca_metabric_all"

myclinicaldata  <-  getClinicalData(mycgds,mycaselist)
myclinicaldata  <-  myclinicaldata[myclinicaldata$HER2_STATUS == "Negative" & myclinicaldata$ER_STATUS == "Negative" &  myclinicaldata$PR_STATUS == "Negative",]
myclinicaldata  <-  myclinicaldata[is.na(myclinicaldata$TUMOR_STAGE) == FALSE,]
myclinicaldata  <-  myclinicaldata[myclinicaldata$TUMOR_STAGE <= 2,]


mygeneticprofile <- 'brca_metabric_cna'
mrna             <- getProfileData(mycgds,c(as.character(symbols)) ,mygeneticprofile,mycaselist)
mrna             <- mrna[complete.cases(mrna),]
mrna.matrix      <- as.matrix(mrna)


expr_avg = apply(mrna.matrix, 1, function(x) median(x, na.rm = T))
expr_clinic <-  merge(x= data.frame(barcode = rownames(myclinicaldata), myclinicaldata), 
                      y= data.frame(barcode = rownames(mrna.matrix),    infect_expr = expr_avg  ), 
                      by = "barcode" 
)

expr_clinic  <- expr_clinic[,c('infect_expr','OS_MONTHS','OS_STATUS')]
expr_clinic  <- expr_clinic[complete.cases(expr_clinic),]


expr_clinic$group <- ifelse( expr_clinic$infect_expr >= 1, "GAIN", "OTHER")


expr_clinic$OS_STATUS_BIN = 1
expr_clinic$OS_STATUS_BIN[expr_clinic$OS_STATUS == "LIVING"] = 0

expr_clinic$group <- factor(expr_clinic$group,levels =  c('OTHER','GAIN'))

my.fit1 <- survfit(Surv(OS_MONTHS, OS_STATUS_BIN) ~ group, data =  expr_clinic) 
km.style <- theme_bw(base_size = 55) + theme(axis.title = element_text( size=30, face="bold"),
                                             axis.text  = element_text( size=55, face="bold"),
                                             plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
                                             axis.line.x = element_line(colour = "black",size = 3),
                                             axis.line.y = element_line(colour = "black",size = 3),
                                             plot.title = element_text(size =40, face = "bold")
) 
p <- ggsurvplot(my.fit1,data =expr_clinic,  pval = F,legend='none',ggtheme =km.style, xlab='',ylab='',palette = c("blue", "red"),size=4)
ggsave(filename = '~/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section4/19p13.12.KM.plot.pdf',width=20,height=15)


#################################################################
# (d):  KM plot of HES4 in basal-like breast cancer
#################################################################

load('server-side/RData//Breast Invasive Carcinoma.RData')
load('client-side/output/Select.pure.sample.breast.cancer.R.output//Select.pure.sample.breast.cancer.RData')
load('client-side/output/DE.breast.cancer.R.output/DE.breast.cancer.RData')


TCGA.CDR.data <- read_excel("client-side/Data/TCGA-CDR/TCGA-CDR-SupplementalTableS1.xlsx", sheet = "TCGA-CDR") %>% as.data.frame
get.patient.id <- function(x) { 
  tmp <- strsplit(x,split = '-') %>% unlist
  paste(tmp[1:3],collapse = '-')
}

PRI.log2.tpm.matrix    <- log2.tpm.matrix[,pure.PRI.breast.cancer.Basal.sample]

g                   <- 'ENSG00000188290'
v                   <- PRI.log2.tpm.matrix[g,] 
high.group          <- sapply(names(v)[ v >= median(v)],get.patient.id)
low.group           <- sapply(names(v)[ v <  median(v)],get.patient.id)
high.group.df       <- TCGA.CDR.data[TCGA.CDR.data$bcr_patient_barcode %in% high.group,]
low.group.df        <- TCGA.CDR.data[TCGA.CDR.data$bcr_patient_barcode %in% low.group,]
high.group.df$group <- 'high'
low.group.df$group  <- 'low'
c.df                <- rbind(high.group.df,low.group.df)
c.df$group          <- factor(c.df$group,levels = c('low','high'))
c.df <- c.df[,c('PFI.time','PFI','group')]
c.df$PFI.time <- c.df$PFI.time / 30 # survival time in months
fit                 <- survfit(Surv(PFI.time,PFI) ~ group , data=c.df)
km.style <- theme_bw(base_size = 55) + theme(axis.title = element_text( size=30, face="bold"),
                                             axis.text  = element_text( size=70, face="bold"),
                                             plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
                                             axis.line.x = element_line(colour = "black",size = 3),
                                             axis.line.y = element_line(colour = "black",size = 3),
                                             plot.title = element_text(size =40, face = "bold")
) 
p <- ggsurvplot(fit=fit,data=c.df,legend='none',ggtheme =km.style, xlab='',ylab='',palette = c("blue", "red"),title= 'HES4',size=5) 
ggsave(filename = '~/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section4/HES4.KM.plot.pdf',width=20,height=15)
sdf                 <- survdiff(Surv(PFI.time,PFI) ~  group, data = c.df)
p.val               <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)



