library(Seurat)
library(dplyr)
library(anndata)

library(reticulate)
use_python('/home/ubuntu/chenlab_deeplearning/chenlab_deeplearning_V2/anaconda3/envs/SCANVI/bin/python', required = T)
sc <- import('scanpy', convert = FALSE)
scvi <- import('scvi', convert = FALSE)
setwd('/home/ubuntu/single_cell/scarches-0.4-2.0/')
sca <- import('scarches', convert = FALSE)
torch <- import('torch', convert = FALSE)
remove_sparsity <- import('scarches.dataset.trvae.data_handling', convert = FALSE)
plt <- import('matplotlib.pyplot', convert = FALSE)
np <- import('numpy', convert = FALSE)
os <- import('os', convert = FALSE)

load("/home/ubuntu/single_cell/Liver_sc_atlas/results/combined_15_datasets_shared_genes_reference_RNA_count_20231009.RData")

RNA_ref <- RNA
RNA_ref[['study']] = ''
RNA_ref$study[grep('HCC', RNA_ref$orig.ident)] = 'GSE149614'
RNA_ref$study[grep('TLH', RNA_ref$orig.ident)] = 'GSE115469'

#Merge query datasets:
RNA <- NULL

for (i in c('GSE242889', 'GSE151530', 'GSE217235', 'GSE234241', 'GSE140228', 'GSE136103',
            'GSE202642', 'GSE125449', 'GSE162616', 'GSE182159', 'GSE169446', 'GSE236382', 'GSE124395')) {
  
  tmp <- get(paste0(i, '_RNA'))
  tmp[['study']] = i
  rownames(tmp@meta.data) <- colnames(tmp@assays$RNA@counts)
  
  if(length(RNA) == 0){
    RNA <- tmp
  }else{
    RNA <- merge(RNA, tmp)
  }
  print(i)
  
}

#--------------------------------------------------------------------------------------------------
#changeable parameters:
#HCC, combined tumor and control group:
save_prediction_query_set_cell_type_file = '/home/ubuntu/single_cell/Liver_sc_atlas/results/liver_combined_13_query_scArches-SCANVI_latent_pred_celltype_reference_2_datasets_20231010.csv'
save_query_set_cell_latent_representations_file = '/home/ubuntu/single_cell/Liver_sc_atlas/results/liver_combined_13_query_RNA_latent_representations_scArches-SCANVI_reference_2_datasets_20231010.csv'
reference_latent_representations_file = '/home/ubuntu/single_cell/Liver_sc_atlas/results/reference_combined_2_datasets_RNA_scArches-SCANVI_latent_representations_20231010.csv'
save_visualize_ref_query_sets_by_cell_type_file = '/home/ubuntu/single_cell/Liver_sc_atlas/results/visualize_reference_combined_2_datasets_query_combined_13_RNA_scArches-SCANVI_by_celltype_20231011.pdf'
save_visualize_ref_query_sets_by_patient_file = '/home/ubuntu/single_cell/Liver_sc_atlas/results/visualize_reference_combined_2_datasets_query_combined_13_RNA_scArches-SCANVI_by_patient_20231011.pdf'
save_visualize_ref_query_sets_by_study_file = '/home/ubuntu/single_cell/Liver_sc_atlas/results/visualize_reference_combined_2_datasets_query_combined_13_RNA_scArches-SCANVI_by_study_20231011.pdf'
save_visualize_query_sets_by_cell_type_file = '/home/ubuntu/single_cell/Liver_sc_atlas/results/visualize_query_combined_13_RNA_scArches-SCANVI_by_celltype_20231011.pdf'
save_visualize_query_sets_by_study_file = '/home/ubuntu/single_cell/Liver_sc_atlas/results/visualize_query_combined_13_RNA_scArches-SCANVI_by_study_20231011.pdf'
save_visualize_query_sets_by_patient_file = '/home/ubuntu/single_cell/Liver_sc_atlas/results/visualize_query_combined_13_RNA_scArches-SCANVI_by_patient_20231011.pdf'
save_visualize_ref_sets_by_cell_type_file = '/home/ubuntu/single_cell/Liver_sc_atlas/results/visualize_ref_combined_2_RNA_scArches-SCANVI_by_celltype_20231011.pdf'
save_visualize_ref_sets_by_study_file = '/home/ubuntu/single_cell/Liver_sc_atlas/results/visualize_ref_combined_2_RNA_scArches-SCANVI_by_study_20231011.pdf'
save_visualize_ref_sets_by_patient_file = '/home/ubuntu/single_cell/Liver_sc_atlas/results/visualize_ref_combined_2_RNA_scArches-SCANVI_by_patient_20231011.pdf'


reference_latent <- read.csv(reference_latent_representations_file, stringsAsFactors = F, row.names = 1)
query_latent <- read.csv(save_query_set_cell_latent_representations_file, stringsAsFactors = F, row.names = 1)

ref_query_latent <- rbind(reference_latent, query_latent)

RNA$cell_type = read.csv(save_prediction_query_set_cell_type_file, stringsAsFactors = F)$x


RNA_ref_query <- merge(RNA_ref, RNA)

rownames(ref_query_latent) = colnames(RNA_ref_query)

RNA_ref_query[['scvi']] <- CreateDimReducObject(embeddings = as.matrix(ref_query_latent), key = "scvi_", assay = DefaultAssay(RNA_ref_query))
# Find clusters, then run UMAP, and visualize
RNA_ref_query <- FindNeighbors(RNA_ref_query, dims = 1:10, reduction = 'scvi')
RNA_ref_query <- FindClusters(RNA_ref_query, resolution =1)

RNA_ref_query <- RunUMAP(RNA_ref_query, dims = 1:10, reduction = 'scvi', n.components = 2)

pdf(save_visualize_ref_query_sets_by_cell_type_file, width = 12, height = 10)
p = DimPlot(RNA_ref_query, reduction = "umap",label=T,group.by='cell_type', raster=FALSE)
print(p)
dev.off()

pdf(save_visualize_ref_query_sets_by_patient_file, width = 20, height = 10)
p = DimPlot(RNA_ref_query, reduction = "umap",label=F,group.by='orig.ident', raster=FALSE)
print(p)
dev.off()


pdf(save_visualize_ref_query_sets_by_study_file, width = 12, height = 10)
p = DimPlot(RNA_ref_query, reduction = "umap",label=F,group.by='study', raster=FALSE)
print(p)
dev.off()

#visualize query RNA alone:
RNA[['scvi']] <- CreateDimReducObject(embeddings = as.matrix(query_latent), key = "scvi_", assay = DefaultAssay(RNA))
# Find clusters, then run UMAP, and visualize
RNA <- FindNeighbors(RNA, dims = 1:10, reduction = 'scvi')
RNA <- FindClusters(RNA, resolution =1)

RNA <- RunUMAP(RNA, dims = 1:10, reduction = 'scvi', n.components = 2)

pdf(save_visualize_query_sets_by_cell_type_file, width = 12, height = 10)
p = DimPlot(RNA, reduction = "umap",label=T,group.by='cell_type', raster=FALSE)
print(p)
dev.off()

pdf(save_visualize_query_sets_by_patient_file, width = 20, height = 10)
p = DimPlot(RNA, reduction = "umap",label=F,group.by='orig.ident', raster=FALSE)
print(p)
dev.off()

pdf(save_visualize_query_sets_by_study_file, width = 12, height = 10)
p = DimPlot(RNA, reduction = "umap",label=F,group.by='study', raster=FALSE)
print(p)
dev.off()

#visualize reference RNA alone:
RNA_ref[['scvi']] <- CreateDimReducObject(embeddings = as.matrix(reference_latent), key = "scvi_", assay = DefaultAssay(RNA_ref))
# Find clusters, then run UMAP, and visualize
RNA_ref <- FindNeighbors(RNA_ref, dims = 1:10, reduction = 'scvi')
RNA_ref <- FindClusters(RNA_ref, resolution =1)

RNA_ref <- RunUMAP(RNA_ref, dims = 1:10, reduction = 'scvi', n.components = 2)

pdf(save_visualize_ref_sets_by_cell_type_file, width = 12, height = 10)
p = DimPlot(RNA_ref, reduction = "umap",label=T,group.by='cell_type', raster=FALSE)
print(p)
dev.off()

pdf(save_visualize_ref_sets_by_patient_file, width = 20, height = 10)
p = DimPlot(RNA_ref, reduction = "umap",label=F,group.by='orig.ident', raster=FALSE)
print(p)
dev.off()

pdf(save_visualize_ref_sets_by_study_file, width = 12, height = 10)
p = DimPlot(RNA_ref, reduction = "umap",label=F,group.by='study', raster=FALSE)
print(p)
dev.off()


