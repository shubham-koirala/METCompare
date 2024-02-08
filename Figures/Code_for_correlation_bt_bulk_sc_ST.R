#bulk metastatic cancer
setwd("Ke_metastatic_work/")
load("TCGA_pancran_data/tcga_panc_symbol.RData") ### This file can be downloaded from https://Octad.org
pheno= read.csv("TCGA_pancran_data/Liver_enriched_cells_octad/phenoDF_withage.csv") ### This file can be downloaded from https://Octad.org

#load("/Users/chenbi12/Documents/msu/tumor_cell_line/data/metastatic/MET500.RData")
colon_liver_metastatsis = pheno[pheno$biopsy.site == "LIVER" & pheno$cancer == "Colon Adenocarcinoma", ]

colon_liver_met_counts = TCGA_panc_symbol_exp[, colon_liver_met$run.id]

mapping = read.csv("~/Desktop/code/gencode.v23.annotation.gene.probeMap.csv") ##This annotation file can be downloaded from https://www.gencodegenes.org/human/release_23lift37.html
mapping = mapping[, c("id", "gene")]

mapping = merge(mapping, colon_liver_met_counts, by.x = "id", by.y = 0)
mapping = mapping[, -1]

colon_liver_metastatsis_counts_merged = aggregate(. ~ gene, mapping, mean)
rownames(colon_liver_metastatsis_counts_merged)= colon_liver_metastatsis_counts_merged[,1]
colon_liver_metastatsis_counts_merged= colon_liver_metastatsis_counts_merged[,-1]

#scRNA-seq
#### The data was obtained from GEO: GSE110009
d = read.csv("Ke_metastatic_work/GSE110009_metastatic_Colon_TPM.txt", sep ="\t")
d = d[!duplicated(d$X), ]
rownames(d) = d$X
d = d[, -1]

library(Seurat)
GSE110009= CreateSeuratObject(counts = d)
for (i in 1:nrow(GSE110009@meta.data)){
  GSE110009@meta.data$sample[i]= strsplit(rownames(GSE110009@meta.data)[i], "_")[[1]][2]
}

GSE110009= NormalizeData(GSE110009)

GSE110009@meta.data$patients_sample= paste(GSE110009@meta.data$orig.ident, GSE110009@meta.data$sample, sep = "_")
Idents(GSE110009)= GSE110009$patients_sample
GSE110009_bulk=AverageExpression(GSE110009, group.by = "patients_sample")
GSE110009_bulk= as.data.frame(GSE110009_bulk$RNA)
GSE110009_bulk_LM= GSE110009_bulk[,c(30, 31, 37, 38, 46,47)]

#test=cbind.data.frame(r = GSE110009_bulk$RNA + 1, c = b@j + 1, x = b@x)

df2 <- as.data.frame(matrix(NA, ncol = ncol(GSE110009_bulk_LM), nrow = nrow(GSE110009_bulk_LM)))
colnames(df2) <- colnames(GSE110009_bulk_LM)
rownames(df2) <- rownames(GSE110009_bulk_LM)
df2= df2[,-c(2,4,6)]

for (i in 1: ncol(GSE110009_bulk_LM)){
  df= GSE110009_bulk_LM[,c(i, i+1)]
  df2=cbind(df2,apply(df, 1, function(x) sum(x)))
  
}

GSE110009_bulk= df2[,c(4,6,8)]
colnames(GSE110009_bulk)=c("P7_LM", "P8_LM", "P9_LM")



bulk_sc = merge(colon_liver_metastatsis_counts_merged, GSE110009_bulk, by = 0)
co_bulk_sc = cor(bulk_sc[, -1], method = "spearman")

plot(bulk_sc$MO_1001.poly.SI_4110.D0CG3ABXX, bulk_sc$P7_LM)

#scRAN-seq 2
library(Seurat)
### The data was obtained from GEO: GSE178318
expression_matrix <- ReadMtx(
  mtx = "Ke_metastatic_work//GSE178318_matrix.mtx.gz", features = "Ke_metastatic_work//GSE178318_genes.tsv.gz",
  cells = "Ke_metastatic_work//GSE178318_barcodes.tsv.gz"
)
seurat_object = CreateSeuratObject(counts = expression_matrix)
cell_name = sapply(colnames(seurat_object@assays$RNA@features), function(x)
  paste(unlist(strsplit(x,  "_"))[-1], collapse = "_")
)
seurat_object=NormalizeData(seurat_object)

for (i in 1:nrow(seurat_object@meta.data)){
  seurat_object@meta.data$sample[i]= strsplit(rownames(seurat_object@meta.data)[i], "_")[[1]][3]
}
seurat_LM= subset(seurat_object, subset=sample=="LM")

for (i in 1:nrow(seurat_LM@meta.data)){
  seurat_LM@meta.data$patients[i]= strsplit(rownames(seurat_LM@meta.data)[i], "_")[[1]][2]
}

Idents(seurat_LM)= seurat_LM$patients
seurat_LM_bulk=AverageExpression(seurat_LM, group.by = "patients")
seurat_LM_bulk= as.data.frame(seurat_LM_bulk$RNA)

bulk_sc_sc2 = merge(bulk_sc, seurat_LM_bulk, by.x = "Row.names", by.y=0)
co_bulk_sc_sc2 = cor(bulk_sc_sc2[, -1], method = "spearman")

#for ST analysis
##The ST data was obtained from PMID: 34919432
load("Ke_metastatic_work/Rama_ST.RData")
d_st = as.matrix(Rama_ST_exp)
cells = rownames(Rama_ST_metadata[Rama_ST_metadata$tissue %in% c("liver2"), ]) #, "liver2", "liver3", "liver4" Colon1

d_st_p1 = log(apply(d_st[, rownames(Rama_ST_metadata[Rama_ST_metadata$tissue %in% c("liver1"), ])], 1, sum) + 1)
d_st_p2 = log(apply(d_st[, rownames(Rama_ST_metadata[Rama_ST_metadata$tissue %in% c("liver2"), ])], 1, sum) + 1)
d_st_p3 = log(apply(d_st[, rownames(Rama_ST_metadata[Rama_ST_metadata$tissue %in% c("liver3"), ])], 1, sum) + 1)
d_st_p4 = log(apply(d_st[, rownames(Rama_ST_metadata[Rama_ST_metadata$tissue %in% c("liver4"), ])], 1, sum) + 1)

bulk_sc_st = merge(bulk_sc_sc2, data.frame(gene = names(d_st_p1), d_st_p1, d_st_p2, d_st_p3, d_st_p4), by.x= "Row.names", by.y = "gene")
names(bulk_sc_st)[c(2:5)]=c("MET1", "MET2", "MET3", "MET4")
co_bulk_sc_st = data.frame(cor(bulk_sc_st_log[, -1], method = "spearman"))

library(reshape2)
co_bulk_sc_st$sample=rownames(co_bulk_sc_st)
#d= melt(liver_cell_atlas_cells_log, id.vars = "sample")
co_bulk_sc_st$sample= factor(co_bulk_sc_st$sample, levels = c("SRR4306096", "SRR4306310", "SRR4306794",    "P7_LM",   "P8_LM",   
                                                              "P9_LM",   "COL07",   "COL12",   "COL15",   "COL16",   "COL17",   "COL18",   "d_st_p1",
                                                              "d_st_p2", "d_st_p3", "d_st_p4"))


d= melt(co_bulk_sc_st, id.vars = "sample")
library(ggplot2)
ggplot(d, aes(variable, forcats::fct_rev(sample), fill = value, size = value)) +
  geom_point(shape = 21, stroke = 0) +
  #geom_hline(yintercept = seq(.5, 4.5, 1), size = .2) +
  scale_x_discrete(position = "top") +
  scale_radius(range = c(2, 15)) +
  scale_fill_gradient(low = "blue", high = "orange", breaks = c(0,  0.50,  1.00),  labels = c("0", "0.50", "1.00"), limits = c(0, 1)) +
  theme_bw() +
  theme(legend.position = "right", 
        panel.grid.major = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        axis.title.x = element_text(size = 6, angle = 45)) +
  guides(size = guide_legend(override.aes = list(fill = NA, color = "black", stroke = .25), 
                             label.position = "bottom",
                             title.position = "right", 
                             order = 1),
         fill = guide_colorbar(ticks.colour = NA, title.position = "top", order = 2)) +
  labs(size = "Correlation", fill = "Correlation", x = NULL, y = NULL)


ggsave("~/Desktop/Ke_metastatic_work/Correlation_bulk_sc_ST.pdf", width = 10, height = 8)

###########
## Scatter plot for bulk, scRNA-seq and ST data
###########

bulk= apply(bulk_sc_st[,c(1:3)], 1, function(x) mean(x))
sc1= apply(bulk_sc_st[,c(4:6)], 1, function(x) mean(x))
sc2=apply(bulk_sc_st[,c(7:12)], 1, function(x) mean(x))
ST=apply(bulk_sc_st[,c(13:16)], 1, function(x) mean(x))


df_mean= data.frame(cbind(bulk, log2(sc1+1), log2(sc2+1), log2(ST+1)))
names(df_mean)= c("Bulk.RNA.seq", "scRNA.seq.1", "scRNA.seq.2", "ST")
samples= c("Bulk.RNA.seq", "scRNA.seq.1", "scRNA.seq.2", "ST")
for (sam in samples){
  ggscatter(df_mean, x = "Bulk.RNA.seq", y = sam, add = "reg.line") +
    stat_cor(label.x = 3, label.y = 14) +
    #stat_regline_equation(label.x = 3, label.y = 12)+
    theme_bw()+theme(axis.text = element_text( size=12, face="bold"),
                     axis.title.y = element_text(size=12, face="bold"),
                     axis.title.x = element_text(size=12, face="bold"),
                     legend.position = "none",
                     plot.title = element_text(size=14, face="bold",lineheight = 0.9),
                     panel.border = element_rect(fill=NA, colour = "black", size=1),
                     axis.ticks = element_line(size = 0.5))+
    labs(title = paste("Bulk.RNA.seq.vs", sam, sep = "."))
  
  ggsave(paste("bulk", sam,"_corr.pdf", sep = "_"), width = 4, height = 4)
}

########### 
##Enrichment of genes shwoing the zero value
###########
genes_with_zero_exp_in_sc <- df_mean[(df_mean$scRNA.seq.1 < 0.01 & df_mean$Bulk.RNA.seq > 0) & (df_mean$scRNA.seq.2 <= 0.01 & df_mean$Bulk.RNA.seq > 0), ]
#genes= read.csv("~/Desktop/Ke_metastatic_work/genes_with_zero_in_sc.csv")
signature_list= list(genes_with_zero_exp_in_sc$genes)
names(signature_list)="Genes"
liver.marker= read.csv("~/Desktop/Liver_sc_Atlas/Refined_celltype_markers.csv", row.names = 1)


genes= liver.marker
cluster_celltype_score = sapply(unique(liver.marker$cluster), function(x){
  idx = genes$cluster==x
  avglogFC = genes$avg_log2FC[idx]
  names(avglogFC) = toupper(genes$gene[idx])
  score_cluster = sapply(signature_list, function(y){
    score = sum(avglogFC[y], na.rm = TRUE) / log2(length(y))
    return(score)
  })
})

colnames(cluster_celltype_score) = as.character(unique(liver.marker$cluster))
min.score=0.1
cellscore_max = apply(cluster_celltype_score, 2, max, na.rm = TRUE)
cellscore_max_celltype = apply(cluster_celltype_score, 2, function(x){
  if (max(x) < min.score){
    return("Others")
  }else{
    return(rownames(cluster_celltype_score)[which.max(x)])
  }
})


cell_score= data.frame(cluster_celltype_score)
cell_score$names= substr(rownames(cell_score), 1, nchar(rownames(cell_score)) - 6)
cell_score$reg= "up"
library(ggplot2)
library(ggthemes)
ggplot(cell_score, aes(x = cluster_celltype_score, y= reorder(names, cluster_celltype_score), color=reg)) + geom_point(aes(size=cluster_celltype_score))+
  ylab(label = NULL)+ xlab(label = "celltype_score")+
  theme_few()+theme(axis.text.x = element_text(size=10, face="bold"), 
                    axis.text.y = element_text(size=10, face="bold"),
                    axis.title =  element_text(size=10, face="bold"),
                    legend.position = "right",
                    panel.border = element_rect(fill=NA, colour = "black", size=0.8), 
                    axis.ticks = element_line(size = 0.5))

ggsave("Genes_with_low_exp_in_scRNA_seq.pdf", width = 7, height = 5)



