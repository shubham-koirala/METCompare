library(Seurat)
library(dplyr)
library(reshape2)

setwd("~/Desktop/infercnv")

#LIVER METASTASIS SINGLE CELL ATLAS DATA
#load ("/home/ubuntu/single_cell/metastasis_study_RNA/cell_atlas/primary_metastasis_combined_cell_atlas/results/TISCH_primary_metastasis_RNA_20231127.RData")

#In this directory save the Single Cell Liver Metastatsis Atlas and Infercnv_agg_avg.txt 
# in this code 'RNA_ref_query' is Single Cell Liver Metastatsis Atlas seurat object

##########################################################################################

# Selecting only malignant cells from Single Cell Liver Metastatsis Atlas based on Infercnv results

##########################################################################################

# Adding Infercnv score to Single Cell Liver Metastatsis Atlas

cnv=read.table("Infercnv_agg_avg.txt",header=T)

comm=intersect(rownames(RNA_ref_query@meta.data),cnv$cellID)

length(comm)

RNA_ref_query$cellname=rownames(RNA_ref_query@meta.data)

RNA_ref_query2=subset(RNA_ref_query, subset =  cellname %in% comm) #subsetting cells found in Infercnv file 

RNA_ref_query2$cnv_score=cnv$score[match(RNA_ref_query2$cellname,cnv$cellID)]

RNA_ref_query2$celltype15=cnv$cellType_threshold_0.015[match(RNA_ref_query2$cellname,cnv$cellID)]

RNA_ref_query2$celltype20=cnv$cellType_threshold_0.02[match(RNA_ref_query2$cellname,cnv$cellID)]

saveRDS (RNA_ref_query2, file = IC_all.RDS)

#############################################################################################

load(IC_all.RDS)

RNA_ref_query2=subset(RNA_ref_query2, subset = celltype15== "Malignant") # Only subsetting 'Malignant cells' based on Infercnv score 0.015

saveRDS (RNA_ref_query2, file = IC_malignant_threshold15.RDS)

######### Checking number of malignant cells in Primary and Metastaic Site of each dataset

load(IC_malignant_threshold15.RDS)

# RNA_ref_query2 seurat object containing only malignant cells 

cell_counts <- as.data.frame(table(RNA_ref_query2$DATASET, RNA_ref_query2$TISSUE_TYPE))

colnames(cell_counts)[1]="DATASET"

cell_counts_table <- dcast(cell_counts, DATASET ~ Var2, value.var = "Freq")

write.csv(cell_counts_table,"Malignant_Cellcount.csv",rownames=F,quote=F)

ds=cell_counts_table[cell_counts_table$PT >10 & cell_counts_table$METASTASIS > 10,c(1)]

###### Creating Primary vs Malignant Signatures from each dataset ##################

for (i in 1:length(ds)) 
   {
  
   try {
    subset <- subset(RNA_ref_query2, subset = DATASET == ds[i])
    
    ds_name=ds[i]
    
    cancer=unique(subset$DISEASE)
    
    seu_name=paste0(ds_name,"_",cancer,"_mali.RDS")
    
    saveRDS (subset, file = seu_name)
    
    Idents(subset)=subset$TISSUE_TYPE
    
    lm_markers <- FindMarkers(subset, ident.1 = "METASTASIS", ident.2 = "PT", test.use = 'wilcox')
    
    lm_markers$log2FoldChange <- lm_markers$avg_log2FC
    
    lm_markers$Symbol <- row.names(lm_markers)
    
    ds_name=ds[i]
    
    cancer=unique(subset$DISEASE)
    
    name=paste0("Signatures/",ds_name,"_",cancer,".csv")
    
    write.csv(lm_markers, name)
   }
}

##### The Signatures and pathway analysis data is provided in Primary_vs_Metastasis_DE.zip
